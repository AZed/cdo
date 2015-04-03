/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2009 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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
#include "modules.h"
#include "dmemory.h"
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
  static char func[] = "getOperator";
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

char *getOperatorName(char *operatorArg)
{
  static char func[] = "getOperatorName";
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
  static char func[] = "makeArgument";
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
  static char func[] = "getFileArg";
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

