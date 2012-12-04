/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Splittime  splithour       Split hours
      Splittime  splitday        Split days
      Splittime  splitmon        Split months
      Splittime  splitseas       Split seasons
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


#define  MAX_STREAMS 32

void *Splittime(void *argument)
{
  int SPLITHOUR, SPLITDAY, SPLITMON, SPLITSEAS;
  int operatorID;
  int operfunc, operintval;
  int nchars;
  int streamID1, streamID2;
  int varID;
  int nrecs;
  int tsID, recID, levelID;
  int vlistID1, vlistID2;
  int  streamIDs[MAX_STREAMS], tsIDs[MAX_STREAMS];
  char filesuffix[32];
  char filename[8192];
  int index = 0;
  int i;
  int taxisID1, taxisID2;
  int vdate, vtime;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  int gridID;
  int nvars, nlevel;
  int nconst;
  double *array = NULL;
  field_t **vars = NULL;
  int season_start;
  const char *seas_name[4];

  cdoInitialize(argument);

  SPLITHOUR = cdoOperatorAdd("splithour", func_time, 10000, NULL);
  SPLITDAY  = cdoOperatorAdd("splitday",  func_date,     1, NULL);
  SPLITMON  = cdoOperatorAdd("splitmon",  func_date,   100, NULL);
  SPLITSEAS = cdoOperatorAdd("splitseas", func_date,   100, NULL);

  operatorID = cdoOperatorID();
  operfunc   = cdoOperatorF1(operatorID);
  operintval = cdoOperatorF2(operatorID);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  season_start = get_season_start();
  get_season_name(seas_name);

  for ( i = 0; i < MAX_STREAMS; i++ ) streamIDs[i] = -1;
  for ( i = 0; i < MAX_STREAMS; i++ ) tsIDs[i] = 0;

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  strcpy(filename, cdoStreamName(1));
  nchars = strlen(filename);

  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), cdoDefaultFileType, vlistID1);

  //  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double *) malloc(gridsize*sizeof(double));
    }

  nvars = vlistNvars(vlistID1);
  nconst = 0;
  for ( varID = 0; varID < nvars; varID++ )
    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) nconst++;

  if ( nconst )
    {
      vars = (field_t **) malloc(nvars*sizeof(field_t *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT )
	    {
	      gridID  = vlistInqVarGrid(vlistID1, varID);
	      nlevel  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      gridsize = gridInqSize(gridID);
		  
	      vars[varID] = (field_t *) malloc(nlevel*sizeof(field_t));

	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  vars[varID][levelID].grid    = gridID;
		  vars[varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		}
	    }
	}
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( operfunc == func_date )
	{
	  index = (vdate/operintval)%100;

	  if ( operatorID == SPLITSEAS )
	    {
	      if ( index < 0 || index > 16 )
		cdoAbort("Month %d out of range!", index);

	      if ( season_start == START_DEC )
		{
		  if ( index <= 12 )
		    index = (index % 12) / 3;
		  else
		    index = index - 13;
		}
	      else
		{
		  if ( index <= 12 )
		    index = (index - 1) / 3;
		  else
		    index = index - 13;
		}
	      
	      if ( index < 0 || index > 3 )
		cdoAbort("Season %d out of range!", index+1);
	    }
	}
      else if ( operfunc == func_time )
	{
	  index = (vtime/operintval)%100;
	}

      if ( index < 0 || index >= MAX_STREAMS )
	cdoAbort("Index out of range!");

      streamID2 = streamIDs[index];
      if ( streamID2 < 0 )
	{
	  if ( operatorID == SPLITSEAS )
	    {
	      sprintf(filename+nchars, "%3s", seas_name[index]);
	      if ( filesuffix[0] )
		sprintf(filename+nchars+3, "%s", filesuffix);
	    }
	  else
	    {
	      sprintf(filename+nchars, "%02d", index);
	      if ( filesuffix[0] )
		sprintf(filename+nchars+2, "%s", filesuffix);
	    }

	  if ( cdoVerbose ) cdoPrint("create file %s", filename);

	  streamID2 = streamOpenWrite(filename, cdoFiletype());

	  streamDefVlist(streamID2, vlistID2);

	  streamIDs[index] = streamID2;
	}

      taxisCopyTimestep(taxisID2, taxisID1);
      streamDefTimestep(streamID2, tsIDs[index]);

      if ( tsID > 0 && tsIDs[index] == 0 && nconst )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT )
		{
		  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    {
		      streamDefRecord(streamID2, varID, levelID);
		      nmiss = vars[varID][levelID].nmiss;
		      streamWriteRecord(streamID2, vars[varID][levelID].ptr, nmiss);
		    }
		}
	    }
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  if ( lcopy && !(tsID == 0 && nconst) )
	    {
	      streamCopyRecord(streamID2, streamID1);
	    }
	  else
	    {
	      streamReadRecord(streamID1, array, &nmiss);
	      streamWriteRecord(streamID2, array, nmiss);

	      if ( tsID == 0 && nconst )
		{
		  if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT )
		    {
		      gridID  = vlistInqVarGrid(vlistID1, varID);
		      gridsize = gridInqSize(gridID);
		      memcpy(vars[varID][levelID].ptr, array, gridsize*sizeof(double));
		      vars[varID][levelID].nmiss = nmiss;
		    }
		}
	    }
	}

      tsIDs[index]++;
      tsID++;
    }

  streamClose(streamID1);

  for ( index = 0; index < MAX_STREAMS; index++ )
    {
      streamID2 = streamIDs[index];
      if ( streamID2 >= 0 ) streamClose(streamID2);
    }
 
  if ( array ) free(array);

  if ( nconst )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTime(vlistID2, varID) == TIME_CONSTANT )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		if ( vars[varID][levelID].ptr )
		  free(vars[varID][levelID].ptr);

	      free(vars[varID]);
	    }
	}

      if ( vars  ) free(vars);
    }

  vlistDestroy(vlistID2);

  cdoFinish();

  return (0);
}
