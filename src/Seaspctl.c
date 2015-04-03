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

      Seaspctl   seaspctl        Seasonal percentiles
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "percentiles.h"
#include "util.h"


void *Seaspctl(void *argument)
{
  static char func[] = "Seaspctl";
  int gridsize;
  int vdate1 = 0, vtime1 = 0;
  int vdate2 = 0, vtime2 = 0;
  int vdate3 = 0, vtime3 = 0;
  int vdate4 = 0, vtime4 = 0;
  int nrecs, nrecords;
  int gridID, varID, levelID, recID;
  int tsID;
  int otsID;
  long nsets;
  int year, month, day, seas, seas0 = 0;
  int streamID1, streamID2, streamID3, streamID4;
  int vlistID1, vlistID2, vlistID3, vlistID4, taxisID1, taxisID2, taxisID3, taxisID4;
  int nmiss;
  int nvars, nlevels;
  int *recVarID, *recLevelID;
  int newseas, oldmon = 0, newmon;
  double missval;
  field_t **vars1 = NULL;
  field_t field;
  int pn;
  HISTOGRAM_SET *hset = NULL;
  int season_start;

  cdoInitialize(argument);

  cdoOperatorAdd("seaspctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number");
  pn = atoi(operatorArgv()[0]);
      
  if ( pn < 1 || pn > 99 )
    cdoAbort("Illegal argument: percentile number %d is not in the range 1..99!", pn);

  season_start = get_season_start();

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));
  streamID3 = streamOpenRead(cdoStreamName(2));
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = streamInqVlist(streamID3);
  vlistID4 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, func_hrd);
  vlistCompare(vlistID1, vlistID3, func_hrd);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisID3 = vlistInqTaxis(vlistID3);
  /* TODO - check that time axes 2 and 3 are equal */
  
  taxisID4 = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID4, taxisID4);

  streamID4 = streamOpenWrite(cdoStreamName(3), cdoFiletype());
  if ( streamID4 < 0 ) cdiError(streamID4, "Open failed on %s", cdoStreamName(3));

  streamDefVlist(streamID4, vlistID4);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int *) malloc(nrecords * sizeof(int));
  recLevelID = (int *) malloc(nrecords * sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);

  field.ptr = (double *) malloc(gridsize*sizeof(double));

  vars1 = (field_t **) malloc(nvars * sizeof(field_t *));
  hset = hsetCreate(nvars);

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevels   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval  = vlistInqVarMissval(vlistID1, varID);

      vars1[varID] = (field_t *) malloc(nlevels * sizeof(field_t));
      hsetCreateVarLevels(hset, varID, nlevels, gridID);

      for ( levelID = 0; levelID < nlevels; levelID++ )
	{
	  vars1[varID][levelID].grid    = gridID;
	  vars1[varID][levelID].nmiss   = 0;
	  vars1[varID][levelID].missval = missval;
	  vars1[varID][levelID].ptr     = (double *) malloc(gridsize * sizeof(double));
	}
    }

  tsID    = 0;
  otsID   = 0;
  while ( TRUE )
    {
      nsets   = 0;
      newseas = FALSE;

      nrecs = streamInqTimestep(streamID2, otsID);
      if ( nrecs != streamInqTimestep(streamID3, otsID) )
        cdoAbort("Number of records in time step %d of %s and %s are different!", otsID+1, cdoStreamName(1), cdoStreamName(2));
      
      vdate2 = taxisInqVdate(taxisID2);
      vtime2 = taxisInqVtime(taxisID2);
      vdate3 = taxisInqVdate(taxisID3);
      vtime3 = taxisInqVtime(taxisID3);
      if ( vdate2 != vdate3 || vtime2 != vtime3 )
        cdoAbort("Verification dates for time step %d of %s and %s are different!", otsID+1, cdoStreamName(1), cdoStreamName(2));

      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID2, &varID, &levelID);
	  streamReadRecord(streamID2, vars1[varID][levelID].ptr, &nmiss);
          vars1[varID][levelID].nmiss = nmiss;
        }
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID3, &varID, &levelID);
	  streamReadRecord(streamID3, field.ptr, &nmiss);
          field.nmiss   = nmiss;
          field.grid    = vars1[varID][levelID].grid;
	  field.missval = vars1[varID][levelID].missval;
	  
	  hsetDefVarLevelBounds(hset, varID, levelID, &vars1[varID][levelID], &field);
        }

      while ( nrecs && (nrecs = streamInqTimestep(streamID1, tsID)) )
	{
	  vdate1 = taxisInqVdate(taxisID1);
	  vtime1 = taxisInqVtime(taxisID1);
	  cdiDecodeDate(vdate1, &year, &month, &day);
	  if ( month < 0 || month > 16 )
	    cdoAbort("Month %d out of range!", month);

	  newmon = month;

	  if ( season_start == START_DEC )
	    {
	      if ( newmon == 12 ) newmon = 0;

	      if ( month <= 12 )
		seas = (month % 12) / 3;
	      else
		seas = month - 13;
	    }
	  else
	    {
	      if ( month <= 12 )
		seas = (month - 1) / 3;
	      else
		seas = month - 13;
	    }

	  if ( seas < 0 || seas > 3 )
	    cdoAbort("Season %d out of range!", seas + 1);

	  if ( nsets == 0 )
	    {
	      seas0 = seas;
	      oldmon = newmon;
	    }

	  if ( newmon < oldmon ) newseas = TRUE;

	  if ( (seas != seas0) || newseas ) break;

	  oldmon = newmon;

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);
	      if ( tsID == 0 )
		{
		  recVarID[recID]   = varID;
		  recLevelID[recID] = levelID;
		}
	      streamReadRecord(streamID1, vars1[varID][levelID].ptr, &nmiss);
	      vars1[varID][levelID].nmiss = nmiss;
	      
	      hsetAddVarLevelValues(hset, varID, levelID, &vars1[varID][levelID]);
	    }

	  vdate4 = vdate1;
	  vtime4 = vtime1;
	  nsets++;
	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      if ( vdate2 != vdate4 )
        cdoAbort("Verification dates for time step %d of %s, %s and %s are different!", otsID+1, cdoStreamName(1), cdoStreamName(2), cdoStreamName(3));
      if ( vtime2 != vtime4 )
        cdoAbort("Verification times for time step %d of %s, %s and %s are different!", otsID+1, cdoStreamName(1), cdoStreamName(2), cdoStreamName(3));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	  nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  
	  for ( levelID = 0; levelID < nlevels; levelID++ )
            hsetGetVarLevelPercentiles(&vars1[varID][levelID], hset, varID, levelID, pn);
	}

      taxisDefVdate(taxisID4, vdate4);
      taxisDefVtime(taxisID4, vtime4);
      streamDefTimestep(streamID4, otsID++);

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  varID   = recVarID[recID];
	  levelID = recLevelID[recID];

	  if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
	    {
	      streamDefRecord(streamID4, varID, levelID);
	      streamWriteRecord(streamID4, vars1[varID][levelID].ptr, vars1[varID][levelID].nmiss);
	    }
	}

      if ( nrecs == 0 ) break;
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevels; levelID++ )
	free(vars1[varID][levelID].ptr);
      free(vars1[varID]);
    }

  free(vars1);
  hsetDestroy(hset);

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID4);
  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
