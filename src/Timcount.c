/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
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

      Timcount    timcount          Time counts
      Hourcount   hourcount         Hourly counts
      Daycount    daycount          Daily counts
      Moncount    moncount          Monthly counts
      Yearcount   yearcount         Yearly counts
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Timcount(void *argument)
{
  int operatorID;
  int cmplen;
  char indate1[DATE_LEN+1], indate2[DATE_LEN+1];
  int gridsize;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int nrecs, nrecords;
  int gridID, varID, levelID, recID;
  int tsID;
  int otsID;
  long nsets;
  int i;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  double missval;
  field_t **vars1 = NULL;
  field_t field;

  cdoInitialize(argument);

  cdoOperatorAdd("timcount",  0, 31, NULL);
  cdoOperatorAdd("yearcount", 0, 10, NULL);
  cdoOperatorAdd("moncount",  0,  8, NULL);
  cdoOperatorAdd("daycount",  0,  6, NULL);
  cdoOperatorAdd("hourcount", 0,  4, NULL);

  operatorID = cdoOperatorID();

  cmplen = DATE_LEN - cdoOperatorF2(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  nvars    = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
      vlistDefVarUnits(vlistID2, varID, "No.");

  if ( cdoOperatorF2(operatorID) == 16 ) vlistDefNtsteps(vlistID2, 1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  nrecords = vlistNrecs(vlistID1);
  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);

  field.ptr = (double *) malloc(gridsize*sizeof(double));

  vars1 = (field_t **) malloc(nvars*sizeof(field_t *));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval  = vlistInqVarMissval(vlistID1, varID);

      vars1[varID] = (field_t *)  malloc(nlevel*sizeof(field_t));

      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  vars1[varID][levelID].grid    = gridID;
	  vars1[varID][levelID].nmiss   = 0;
	  vars1[varID][levelID].missval = missval;
	  vars1[varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
	}
    }

  tsID    = 0;
  otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
	{
	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);

	  if ( nsets == 0 ) SET_DATE(indate2, vdate, vtime);
	  SET_DATE(indate1, vdate, vtime);

	  if ( DATE_IS_NEQ(indate1, indate2, cmplen) ) break;

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);

	      if ( tsID == 0 )
		{
		  recVarID[recID]   = varID;
		  recLevelID[recID] = levelID;
		}

	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	      if ( nsets == 0 )
		{
		  for ( i = 0; i < gridsize; i++ )
		    vars1[varID][levelID].ptr[i] = vars1[varID][levelID].missval;
		  vars1[varID][levelID].nmiss = gridsize;
		}

              streamReadRecord(streamID1, field.ptr, &field.nmiss);
              field.grid    = vars1[varID][levelID].grid;
	      field.missval = vars1[varID][levelID].missval;

              farcount(&vars1[varID][levelID], field);
	    }

	  vdate0 = vdate;
	  vtime0 = vtime;
	  nsets++;
	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      taxisDefVdate(taxisID2, vdate0);
      taxisDefVtime(taxisID2, vtime0);
      streamDefTimestep(streamID2, otsID);

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  varID   = recVarID[recID];
	  levelID = recLevelID[recID];

	  if ( otsID && vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, vars1[varID][levelID].ptr,  vars1[varID][levelID].nmiss);
	}

      if ( nrecs == 0 ) break;
      otsID++;
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  free(vars1[varID][levelID].ptr);
	}

      free(vars1[varID]);
    }

  free(vars1);

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
