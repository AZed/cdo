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

      Ymonpctl   ymonpctl        Multi-year monthly percentiles
*/


#include <stdio.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "field.h"
#include "percentiles.h"

#define  NMONTH     17

void *Ymonpctl(void *argument)
{
  static char func[] = "Ymonpctl";
  int gridsize;
  int varID;
  int recID;
  int gridID;
  int vdate, vtime;
  int year, month, day;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  long nsets[NMONTH];
  int streamID1, streamID2, streamID3, streamID4;
  int vlistID1, vlistID2, vlistID3, vlistID4, taxisID1, taxisID2, taxisID3, taxisID4;
  int nmiss;
  int nvars, nlevels;
  int *recVarID, *recLevelID;
  int vdates1[NMONTH], vtimes1[NMONTH];
  int vdates2[NMONTH], vtimes2[NMONTH];
  double missval;
  FIELD **vars1[NMONTH];
  FIELD field;
  int pn;
  HISTOGRAM_SET *hsets[NMONTH];

  cdoInitialize(argument);
  cdoOperatorAdd("ymonpctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number");
  pn = atoi(operatorArgv()[0]);
      
  if ( pn < 1 || pn > 99 )
    cdoAbort("Illegal argument: percentile number %d is not in the range 1..99!", pn);

  for ( month = 0; month < NMONTH; month++ )
    {
      vars1[month] = NULL;
      hsets[month] = NULL;
      nsets[month] = 0;
    }

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

  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);
  field.ptr = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID2, tsID)) )
    {
      if ( nrecs != streamInqTimestep(streamID3, tsID) )
        cdoAbort("Number of records in time step %d of %s and %s are different!", tsID+1, cdoStreamName(1), cdoStreamName(2));
      
      vdate = taxisInqVdate(taxisID2);
      vtime = taxisInqVtime(taxisID2);
      
      if ( vdate != taxisInqVdate(taxisID3) || vtime != taxisInqVtime(taxisID3) )
        cdoAbort("Verification dates for time step %d of %s and %s are different!", tsID+1, cdoStreamName(1), cdoStreamName(2));
        
      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      decode_date(vdate, &year, &month, &day);
      if ( month < 0 || month >= NMONTH )
	cdoAbort("Month %d out of range!", month);

      vdates2[month] = vdate;
      vtimes2[month] = vtime;

      if ( vars1[month] == NULL )
	{
	  vars1[month] = (FIELD **) malloc(nvars*sizeof(FIELD *));
          hsets[month] = hsetCreate(nvars);

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevels  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      missval  = vlistInqVarMissval(vlistID1, varID);

	      vars1[month][varID] = (FIELD *)  malloc(nlevels*sizeof(FIELD));
              hsetCreateVarLevels(hsets[month], varID, nlevels, gridID);
	      
	      for ( levelID = 0; levelID < nlevels; levelID++ )
		{
		  vars1[month][varID][levelID].grid    = gridID;
		  vars1[month][varID][levelID].nmiss   = 0;
		  vars1[month][varID][levelID].missval = missval;
		  vars1[month][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		}
	    }
	}
      
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID2, &varID, &levelID);
	  streamReadRecord(streamID2, vars1[month][varID][levelID].ptr, &nmiss);
          vars1[month][varID][levelID].nmiss = nmiss;
        }
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID3, &varID, &levelID);
	  streamReadRecord(streamID3, field.ptr, &nmiss);
          field.nmiss   = nmiss;
          field.grid    = vars1[month][varID][levelID].grid;
	  field.missval = vars1[month][varID][levelID].missval;
	  
	  hsetDefVarLevelBounds(hsets[month], varID, levelID, &vars1[month][varID][levelID], &field);
        }
      
      tsID++;
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      decode_date(vdate, &year, &month, &day);
      if ( month < 0 || month >= NMONTH )
	cdoAbort("Month %d out of range!", month);

      vdates1[month] = vdate;
      vtimes1[month] = vtime;

      if ( vars1[month] == NULL )
        cdoAbort("No data for month %d in %s and %s", month, cdoStreamName(1), cdoStreamName(2));

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  recVarID[recID]   = varID;
	  recLevelID[recID] = levelID;

	  streamReadRecord(streamID1, vars1[month][varID][levelID].ptr, &nmiss);
	  vars1[month][varID][levelID].nmiss = nmiss;

	  hsetAddVarLevelValues(hsets[month], varID, levelID, &vars1[month][varID][levelID]);
	}

      nsets[month]++;
      tsID++;
    }

  otsID = 0;
  for ( month = 0; month < NMONTH; month++ )
    if ( nsets[month] )
      {
        if ( vdates1[month] != vdates2[month] )
          cdoAbort("Verification dates for month %d of %s, %s and %s are different!", month, cdoStreamName(1), cdoStreamName(2), cdoStreamName(3));
        if ( vtimes1[month] != vtimes2[month] )
          cdoAbort("Verification times for month %d of %s, %s and %s are different!", month, cdoStreamName(1), cdoStreamName(2), cdoStreamName(3));

	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      
	    for ( levelID = 0; levelID < nlevels; levelID++ )
	      hsetGetVarLevelPercentiles(&vars1[month][varID][levelID], hsets[month], varID, levelID, pn);
	  }

	taxisDefVdate(taxisID4, vdates1[month]);
	taxisDefVtime(taxisID4, vtimes1[month]);
	streamDefTimestep(streamID4, otsID++);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
	      {
		streamDefRecord(streamID4, varID, levelID);
		streamWriteRecord(streamID4, vars1[month][varID][levelID].ptr, vars1[month][varID][levelID].nmiss);
	      }
	  }
      }

  for ( month = 0; month < NMONTH; month++ )
    {
      if ( vars1[month] != NULL )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevels; levelID++ )
		free(vars1[month][varID][levelID].ptr);
	      free(vars1[month][varID]);
	    }
	  free(vars1[month]); 
	  hsetDestroy(hsets[month]);
	}
    }

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
