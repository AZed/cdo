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

#include <stdio.h>
#include <string.h>
#include <ctype.h>   /* tolower */

#include "cdo.h"
#include "cdo_int.h"
#include "modules.h"
#include "util.h"

char *getProgname(char *string)
{
  char *progname;

#if defined (_WIN32)
  /*  progname = strrchr(string, '\\'); */
  progname = " cdo";
#else
  progname = strrchr(string, '/');
#endif

  if ( progname == NULL ) progname = string;
  else                    progname++;

  return (progname);
}

char *getOperator(const char *argument)
{
  char *operatorArg = NULL;
  char *blankpos;
  size_t len;

  if ( argument )
    {
      blankpos = strchr(argument, ' ');

      if ( blankpos )
	len = blankpos - argument;
      else
	len = strlen(argument);

      operatorArg = (char *) malloc(len+1);

      memcpy(operatorArg, argument, len);
      operatorArg[len] = '\0';
    }

  return (operatorArg);
}

char *operatorAlias(char *operatorName);

char *getOperatorName(const char *operatorArg)
{
  char *commapos;
  char *operatorName = NULL;
  size_t len;

  if ( operatorArg )
    {
      if ( operatorArg[0] == '-' ) operatorArg++;

      commapos = strchr(operatorArg, ',');

      if ( commapos )
	len = commapos - operatorArg;
      else
	len = strlen(operatorArg);

      operatorName = (char *) malloc(len+1);

      memcpy(operatorName, operatorArg, len);
      operatorName[len] = '\0';
    }

  /*  return (operatorName); */
  return (operatorAlias(operatorName));
}


char *makeArgument(int argc, char *argv[])
{
  char *argument = NULL;
  int iarg;
  size_t len, pos = 0, off = 0;

  if ( argv[0][0] == '-' ) off = 1;
  for ( iarg = 0; iarg < argc; iarg++ )
    {
      len = strlen(argv[iarg]) + 1 - off;
      argument = (char *) realloc(argument, pos+len);
      strcpy(&argument[pos], argv[iarg]+off);
      pos += len;
      argument[pos-1] = ' ';
      off = 0;
    }

  if ( argc )
    argument[pos-1] = '\0';

  return (argument);
}


char *getFileArg(char *argument)
{
  char *fileArg = NULL;
  char *parg;
  char *blankpos;
  size_t len;

  if ( argument )
    {
      blankpos = strchr(argument, ' ');

      if ( blankpos )
	{
	  parg = blankpos + 1;
	  len = strlen(parg);
	  fileArg = (char *) malloc(len+1);
	  strcpy(fileArg, parg);
	}
    }

  return (fileArg);
}


void input_int(char *arg, int intarr[], int maxint, int *nintfound)
{
  int nint = 0;

  intarr[nint++] = atoi(arg);

  while ( (arg = strchr(arg, ',')) && (nint < maxint) )
    intarr[nint++] = atoi(++arg);
    
  *nintfound = nint;
}


void strtolower(char *str)
{
  int i, len;

  if ( str )
    {
      len = (int) strlen(str);
      for ( i = 0; i < len; i++ )
	str[i] = tolower((int) str[i]);
    }
}


const char *seas_name_dec[4] = {"DJF", "MAM", "JJA", "SON"};
const char *seas_name_jan[4] = {"JFM", "AMJ", "JAS", "OND"};

int get_season_start(void)
{
  int season_start = START_DEC;
  char *envstr;

  envstr = getenv("CDO_SEASON_START");
  if ( envstr )
    {
      if      ( strcmp(envstr, "DEC") == 0 ) season_start = START_DEC;
      else if ( strcmp(envstr, "JAN") == 0 ) season_start = START_JAN;
      
      if ( cdoVerbose )
	{
	  if      ( season_start == START_DEC )
	    cdoPrint("Set SEASON_START to December");
	  else if ( season_start == START_JAN )
	    cdoPrint("Set SEASON_START to January");
	}
    }

  return (season_start);
}


void get_season_name(const char *seas_name[4])
{
  long i;

  if ( get_season_start() == START_DEC )
    for ( i = 0; i < 4; ++i ) seas_name[i] = seas_name_dec[i];
  else
    for ( i = 0; i < 4; ++i ) seas_name[i] = seas_name_jan[i];
}


//#include <sys/types.h>
#include <sys/stat.h>
//#include <unistd.h>

int fileExist(const char *filename)
{
  int status = 0;
  struct stat buf;

  if ( stat(filename, &buf) == 0 )
    {
      if ( buf.st_size > 0 ) status = 1;
    }

  return (status);
}


int userFileOverwrite(const char *filename)
{
  int status = 0;
  char line[1024], *pline;

  fprintf(stderr, "File %s already exist, overwrite? (yes/no): ", filename);
  readline(stdin, line, 1024);
  pline = line;
  while ( isspace((int) *pline) ) pline++;
  if ( pline[0] == 'y' && pline[1] == 'e' && pline[2] == 's' )
    status = 1;
  else if ( pline[0] == 'Y' && pline[1] == 'E' && pline[2] == 'S' )
    status = 1;

  return (status);
}

int stdin_is_tty  = 0;
int stdout_is_tty = 0;

void init_is_tty(void)
{
  struct stat statbuf;
  fstat(0, &statbuf);
  if ( S_ISCHR(statbuf.st_mode) ) stdin_is_tty = 1;  
  fstat(1, &statbuf);
  if ( S_ISCHR(statbuf.st_mode) ) stdout_is_tty = 1;  
}


int ps_lhead = FALSE;
int ps_nch   = 0;
int ps_cval  = -1;

void progressInit(void)
{
  ps_lhead = FALSE;
  ps_nch   = 0;;
  ps_cval  = -1;
}


void progressStatus(double offset, double refval, double curval)
{
  int ival;

  if ( !stdout_is_tty ) return;

  offset = offset < 0 ? 0: offset;
  offset = offset > 1 ? 1: offset;
  refval = refval < 0 ? 0: refval;
  refval = refval > 1 ? 1: refval;
  curval = curval < 0 ? 0: curval;
  curval = curval > 1 ? 1: curval;

  ival = (offset + refval*curval)*100;

  if ( ps_cval == -1 )
    {
      ps_nch = fprintf(stdout, "%s: %3d%%", processInqPrompt(), 0);
      fflush(stdout);
      ps_lhead = TRUE;
    }

  if ( ival != ps_cval )
    {
      ps_cval = ival;
      fprintf(stdout, "\b\b\b\b%3d%%", ps_cval);
      fflush(stdout);
    }

  if ( ps_cval == 100 && ps_lhead )
    {
      ps_lhead = FALSE;
      while ( ps_nch-- ) fprintf(stdout, "\b \b");
      fflush(stdout);
    }
}
