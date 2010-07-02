/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2010 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:
*/

#include <ctype.h>  /* isspace */

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "counter.h"

int streamID = 0;
int vlistID = 0;
int varID = 0;
int levelID = 0;
int tsID1 = 0;
int tsID2 = 59;
double *array = NULL;

#define MAX_LINE 256

int Done = 0;

int com_help(char *);
int com_list(char *);
int com_quit(char *);
int com_stat(char *);
int com_vars(char *);
//int com_stat(char *);


typedef struct {
  char  *name; /* User printable name of the function. */
  int  (*func)(char *); /* Function to call to do the job. */
  char  *doc; /* Documentation for this function. */
}
command_t;

command_t commands[] = {
  { "help", com_help, "Display this text" },
  { "?",    com_help, "Synonym for 'help'" },
  { "list", com_list, "List files in DIR" },
  { "quit", com_quit, "Quit using CDO" },
  { "stat", com_stat, "Statistic for selected field" },
  { "vars", com_vars, "list variables" },
  //  { "stat", com_stat, "Print out statistics on FILE" },
  { NULL, NULL, NULL }
};


/* Return non-zero if ARG is a valid argument for CALLER, else print
an error message and return zero. */
int valid_argument (char *caller, char *arg)
{
  if (!arg || !*arg)
    {
      fprintf (stderr, "%s: Argument required.\n", caller);
      return (0);
    }
  return (1);
}

/* Print out help for ARG, or for all of the commands if ARG is not present. */
int com_help(char *arg)
{
  int i;
  int printed = 0;

  for ( i = 0; commands[i].name; i++ )
    {
      if (!*arg || (strcmp (arg, commands[i].name) == 0))
	{
	  printf ("%s\t\t%s.\n", commands[i].name, commands[i].doc);
	  printed++;
	}
    }

  if ( !printed )
    {
      printf ("No commands match â%sâ. Possibilties are:\n", arg);
      for (i = 0; commands[i].name; i++)
	{
	  /* Print in six columns. */
	  if ( printed == 6 )
	    {
	      printed = 0;
	      printf ("\n");
	    }
	  printf("%s\t", commands[i].name);
	  printed++;
	}

      if (printed) printf ("\n");
    }

  return (0);
}


/* List the file(s) named in arg. */
int com_list(char *arg)
{
  if (!arg)
    arg = "";

  return (0);
}

/* The user wishes to quit using this program. Just set DONE non-zero. */
int com_quit(char *arg)
{
  Done = 1;

  return (0);
}


int com_stat(char *arg)
{
  int nrecs;
  int tsID;

  for ( tsID = tsID1; tsID <= tsID2; ++tsID )
    {
      nrecs = streamInqTimestep(streamID, tsID);
      if ( nrecs == 0 )
	{
	  fprintf(stderr, "Timestep %d out of range!\n", tsID+1);
	  break;
	}
      else
	{
	  int i;
	  int nmiss;
	  int gridsize;
	  double fmin = 1.e50 , fmax = -1.e50, fmean = 0;
	  counter_t counter;
	  
	  counter_start(&counter);
	  streamReadVarSlice(streamID, varID, levelID, array, &nmiss);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
	  for ( i = 0; i < gridsize; ++i )
	    {
	      if ( array[i] < fmin ) fmin = array[i];
	      if ( array[i] > fmax ) fmax = array[i];
	      fmean += array[i];
	    }
	  fmean /= gridsize;
	  counter_stop(&counter);
	  
	  fprintf(stdout, "%d %d %d %g %g %g (%gs)\n",
		  tsID+1, varID+1, levelID+1, fmin, fmean, fmax,
		  counter_cputime(counter));
	}
    }
  return (0);
}


int com_vars(char *arg)
{
  int varID;
  int nvars = vlistNvars(vlistID);

  if (!arg)
    arg = "";
  printf("com_vars: %s %d\n", arg, nvars);

  for ( varID = 0; varID < nvars; ++varID )
    {
      fprintf(stdout,"varID %d\n", varID);
    }

  return (0);
}

/* Look up NAME as the name of a command, and return a pointer to that
command. Return a NULL pointer if NAME isn't a command name. */
command_t *find_command(char *name)
{
  int i;

  for ( i = 0; commands[i].name; i++ )
    if ( strcmp(name, commands[i].name) == 0)
      return (&commands[i]);

  return ((command_t *)NULL);
}

/* Execute a command line. */
int execute_line(char *line)
{
  int i;
  command_t *command;
  char *word;

  /* Isolate the command word. */
  i = 0;
  while ( line[i] && isspace(line[i]))  i++;
  word = line + i;
  while ( line[i] && !isspace(line[i])) i++;

  if ( line[i] ) line[i++] = '\0';

  command = find_command(word);
  if ( !command )
    {
      fprintf (stderr, "%s: No such command!\n", word);
      return (-1);
    }
  /* Get argument to command, if any. */
  while ( isspace(line[i]) ) i++;

  word = line + i;
  /* Call the function. */
  return ((*(command->func)) (word));
}

/* Strip isspace from the start and end of STRING. Return a pointer into STRING. */
char *stripwhite(char *string)
{
  char *s, *t;
  for (s = string; isspace(*s); s++)
    ;
  if (*s == 0)
    return (s);
  t = s + strlen (s) - 1;
  while (t > s && isspace(*t))
    t--;
  *++t = '\0';

  return s;
}


void readcmd(const char *prompt, char *line, int size)
{
  fputs(prompt, stdout);
  fflush(stdout);

  *line = '\0';
  if ( fgets(line, size, stdin) )
    {
      char *newline = strchr(line, '\n'); /* check for trailing '\n' */
      if ( newline )
	*newline = '\0'; /* overwrite the '\n' with a terminating null */
    }
}



void *Command(void *argument)
{
  static char func[] = "Command";
  int nrecs;
  int recID, varID, levelID;
  int gridsize;
  int nmiss;
  int taxisID;
  double s_utime, s_stime;
  double e_utime, e_stime;
  double c_cputime = 0, c_usertime = 0, c_systime = 0;
  char line[MAX_LINE];
  char *s;
 
  cdoInitialize(argument);

  processStartTime(&s_utime, &s_stime);

  streamID = streamOpenRead(cdoStreamName(0));
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));
  
  vlistID = streamInqVlist(streamID);
  taxisID = vlistInqTaxis(vlistID);

  gridsize = vlistGridsizeMax(vlistID);
  array = (double *) malloc(gridsize*sizeof(double));

  /* Loop reading and executing lines until the user quits. */
  while ( !Done )
    {
      readcmd("cdo cmd> ", line, MAX_LINE);

      /* Remove leading and trailing whitespace from the line.
	 Then, if there is anything left, add it to the history list
	 and execute it. */
      s = stripwhite(line);
      if ( *s ) execute_line(s);
    }


/*
  while ( readline(stdin, line, MAX_LINE) )
    {
      linep = line;
      
      if      ( strcmp(linep, "exit") == 0 ) break;
      else if ( strcmp(linep, "vars") == 0 )
	{
	  int varID;
	  int nvars = vlistNvars(vlistID);
	  for ( varID = 0; varID < nvars; ++nvars )
	    {
	      fprintf(stdout,"varID %d\n", varID);
	    }
	}
    }
*/
  
  streamClose(streamID);

  if ( array ) free(array);
  
  cdoProcessTime(&e_utime, &e_stime);

  c_usertime = e_utime - s_utime;
  c_systime  = e_stime - s_stime;
  c_cputime  = c_usertime + c_systime;

  s_utime = e_utime;
  s_stime = e_stime;

  cdoPrint("%.2fs %.2fs %.2fs", c_usertime, c_systime, c_cputime);

  cdoFinish();

  return (0);
}
