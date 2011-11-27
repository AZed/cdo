/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <string.h>
#include <time.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"

static char *ghistory = NULL;
static size_t ghistorysize = 0;

static char strtime[32];

void init_strtime()
{
  time_t tp;
  struct tm *ltime;

  tp = time(NULL);

  if ( tp != -1 )
    {
      ltime = localtime(&tp);
      (void) strftime(strtime, sizeof(strtime), "%a %b %d %H:%M:%S %Y: ", ltime);
    }
}


char *get_strtimeptr()
{
  if ( strlen(strtime) == 0 )
    init_strtime();

  return (strtime);
}


void cdoInqHistory(int fileID)
{
  if ( ghistory )
    {
      free(ghistory);
      ghistorysize = 0;
      ghistory = NULL;
    }

  ghistorysize = streamInqHistorySize(fileID);
  if ( ghistorysize > 0 )
    {
      size_t len;
      ghistory = (char *) malloc(ghistorysize+1);
      ghistory[ghistorysize] = 0;
      streamInqHistoryString(fileID, ghistory);
      len = strlen(ghistory);
      if ( len < ghistorysize )
	{
	  /* printf("%d %d\n", len, ghistorysize); */
	  ghistorysize = len;
	}
    }
}


void cdoDefHistory(int fileID, char *histstring)
{
  char *history = NULL;
  size_t historysize = 0;
  char *strtimeptr;
  extern int cdoDisableHistory;

  strtimeptr = get_strtimeptr();
  
  historysize = ghistorysize+strlen(strtimeptr)+strlen(histstring)+2;
  history = (char *) malloc(historysize);

  strcpy(history, strtimeptr);
  strcat(history, histstring);

  if ( cdoDisableHistory == FALSE )
    if ( ghistory )
      {
	strcat(history, "\n");
	strcat(history, ghistory);
      }

  streamDefHistory(fileID, strlen(history), history);
  
  free(history);
}
