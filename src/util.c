/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 /* ftello */
#endif

#include <stdio.h>
#include <string.h>
#include <ctype.h>   /* tolower */

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "modules.h"
#include "util.h"


char *getProgname(char *string)
{
  char *progname;

#if defined(_WIN32)
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
  size_t len;

  if ( argument )
    {
      len = 1 + strlen(argument);

      operatorArg = (char*) malloc(len);

      memcpy(operatorArg, argument, len);
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

      operatorName = (char*) malloc(len+1);

      memcpy(operatorName, operatorArg, len);
      operatorName[len] = '\0';
    }

  /*  return (operatorName); */
  return (operatorAlias(operatorName));
}


argument_t *file_argument_new(const char *filename)
{
  argument_t *argument;

  argument = (argument_t*) calloc(1, sizeof(argument_t));

  argument->argc = 1;
  argument->argv = (char **) calloc(1, sizeof(char *));
  argument->argv[0] = (char *) filename;
  argument->args = (char *) filename;

  return (argument);
}


void file_argument_free(argument_t *argument)
{
  if ( argument )
    {
      if ( argument->argc )
	{
	  assert(argument->argc == 1);
	  free(argument->argv);
	}
      free(argument);
    }
}


argument_t *argument_new(size_t argc, size_t len)
{
  argument_t *argument;

  argument = (argument_t*) calloc(1, sizeof(argument_t));

  if ( argc > 0 )
    {
      argument->argc = argc;
      argument->argv = (char **) calloc(argc, sizeof(char *));
    }

  if ( len > 0 )
    argument->args = (char*) calloc(len, sizeof(char));

  return (argument);
}


void argument_free(argument_t *argument)
{
  if ( argument )
    {
      if ( argument->argc )
	{
	  int argc =  argument->argc;
	  for ( int i = 0; i < argc; ++i )
	    {
	      if ( argument->argv[i] )
		{
		  free(argument->argv[i]);
		  argument->argv[i] = NULL;
		}
	    }

	  free(argument->argv);
	  argument->argv = NULL;
	  argument->argc = 0;
	}

      if ( argument->args )
	{
	  free(argument->args);
	  argument->args = NULL;
	}

      free(argument);
    }
}


void argument_fill(argument_t *argument, int argc, char *argv[])
{
  int iarg;

  assert(argument->argc == argc);

  for ( iarg = 0; iarg < argc; ++iarg )
    argument->argv[iarg] = strdup(argv[iarg]);
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
	  fileArg = (char*) malloc(len+1);
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


void get_season_name(const char *seas_name[])
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

int fileExists(const char *filename)
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
  int status = 0, len;
  char line[1024], *pline;

  fprintf(stderr, "File %s already exists, overwrite? (yes/no): ", filename);
  readline(stdin, line, 1024);
  pline = line;
  while ( isspace((int) *pline) ) pline++;
  len = strlen(pline);
  if ( len == 3 )
    {
      if ( pline[0] == 'y' && pline[1] == 'e' && pline[2] == 's' )
	status = 1;
      else if ( pline[0] == 'Y' && pline[1] == 'E' && pline[2] == 'S' )
	status = 1;
    }
  else if ( len == 1 )
    {
      if ( pline[0] == 'y' ) status = 1;
    }

  return (status);
}


int ps_lhead = FALSE;
int ps_nch   = 0;
int ps_cval  = -1;

void progressInit(void)
{
  ps_lhead = FALSE;
  ps_nch   = 0;
  ps_cval  = -1;
}


void progressStatus(double offset, double refval, double curval)
{
  int ival;

  if ( cdoSilentMode ) return;
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


int datatype2str(int datatype, char *datatypestr)
{
  int status = 0;

  if      ( datatype == DATATYPE_PACK   ) strcpy(datatypestr, "P0");
  else if ( datatype > 0 && datatype <= 32  ) sprintf(datatypestr, "P%d", datatype);
  else if ( datatype == DATATYPE_CPX32  ) strcpy(datatypestr, "C32");
  else if ( datatype == DATATYPE_CPX64  ) strcpy(datatypestr, "C64");
  else if ( datatype == DATATYPE_FLT32  ) strcpy(datatypestr, "F32");
  else if ( datatype == DATATYPE_FLT64  ) strcpy(datatypestr, "F64");
  else if ( datatype == DATATYPE_INT8   ) strcpy(datatypestr, "I8");
  else if ( datatype == DATATYPE_INT16  ) strcpy(datatypestr, "I16");
  else if ( datatype == DATATYPE_INT32  ) strcpy(datatypestr, "I32");
  else if ( datatype == DATATYPE_UINT8  ) strcpy(datatypestr, "U8");
  else if ( datatype == DATATYPE_UINT16 ) strcpy(datatypestr, "U16");
  else if ( datatype == DATATYPE_UINT32 ) strcpy(datatypestr, "U32");
  else                                  { strcpy(datatypestr, "-1"); status = -1;}

  return (status);
}


int str2datatype(const char *datatypestr)
{
  int datatype = -1;
  size_t len;

  len = strlen(datatypestr);

  if ( len > 1 )
    {
      int ilen = atoi(datatypestr+1);
      if      ( memcmp(datatypestr, "P0",  len) == 0 ) datatype = DATATYPE_PACK;
      else if ( memcmp(datatypestr, "P",     1) == 0 &&
		ilen > 0 && ilen <= 32 )               datatype = atoi(datatypestr+1);
      else if ( memcmp(datatypestr, "C32", len) == 0 ) datatype = DATATYPE_CPX32;
      else if ( memcmp(datatypestr, "C64", len) == 0 ) datatype = DATATYPE_CPX64;
      else if ( memcmp(datatypestr, "F32", len) == 0 ) datatype = DATATYPE_FLT32;
      else if ( memcmp(datatypestr, "F64", len) == 0 ) datatype = DATATYPE_FLT64;
      else if ( memcmp(datatypestr, "I8",  len) == 0 ) datatype = DATATYPE_INT8;
      else if ( memcmp(datatypestr, "I16", len) == 0 ) datatype = DATATYPE_INT16;
      else if ( memcmp(datatypestr, "I32", len) == 0 ) datatype = DATATYPE_INT32;
      else if ( memcmp(datatypestr, "U8",  len) == 0 ) datatype = DATATYPE_UINT8;
      else if ( memcmp(datatypestr, "U16", len) == 0 ) datatype = DATATYPE_UINT16;
      else if ( memcmp(datatypestr, "U32", len) == 0 ) datatype = DATATYPE_UINT32;
      else if ( memcmp(datatypestr, "real",   len) == 0 ) datatype = DATATYPE_FLT32;
      else if ( memcmp(datatypestr, "double", len) == 0 ) datatype = DATATYPE_FLT64;
    }

  return (datatype);
}


off_t filesize(const char *filename)
{
  FILE *fp;
  off_t pos = 0;

  if ( filename[0] == '(' && filename[1] == 'p' )
    {
    }
  else
    {
      fp = fopen(filename, "r");
      if ( fp == NULL )
	{
	  fprintf(stderr, "Open failed on %s\n", filename);
	}
      else
	{
	  fseek(fp, 0L, SEEK_END);
	  pos = ftello(fp);
	}
    }
  
  return pos;
}
