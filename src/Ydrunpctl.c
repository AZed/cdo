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

      Ydrunpctl    ydrunpctl         Multi-year daily running percentiles
*/

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "functs.h"
#include "field.h"
#include "percentiles.h"


#define NDAY 373


void *Ydrunpctl(void *argument)
{
  static char func[] = "Ydrunpctl";
  int gridsize;
  int varID;
  int recID;
  int gridID;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  int inp, its, ndates = 0;
  int streamID1, streamID2, streamID3, streamID4;
  int vlistID1, vlistID2, vlistID3, vlistID4;
  int nmiss;
  int nvars, nlevels;
  int *recVarID, *recLevelID;
  double missval;
  field_t ***vars1 = NULL, **vars2[NDAY];
  datetime_t *datetime;
  int taxisID1, taxisID2, taxisID3, taxisID4;
  int calendar, dpy;
  int vdate, vtime;
  int vdates1[NDAY], vtimes1[NDAY];
  int vdates2[NDAY], vtimes2[NDAY];
  int nsets[NDAY];
  int year, month, day, dayoy;
  field_t field;
  int pn;
  HISTOGRAM_SET *hsets[NDAY];
    
  cdoInitialize(argument);
  cdoOperatorAdd("ydrunpctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number, number of timesteps");
  operatorCheckArgc(2);
  pn     = atoi(operatorArgv()[0]);
  ndates = atoi(operatorArgv()[1]);

  if ( pn < 1 || pn > 99 )
    cdoAbort("Illegal argument: percentile number %d is not in the range 1..99!", pn);
  
  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    {
      vars2[dayoy] = NULL;
      hsets[dayoy] = NULL;
      nsets[dayoy] = 0;
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

  taxisID4 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID4, taxisID4);

  calendar = taxisInqCalendar(taxisID1);
  dpy      = calendar_dpy(calendar);

  streamID4 = streamOpenWrite(cdoStreamName(3), cdoFiletype());
  if ( streamID4 < 0 ) cdiError(streamID4, "Open failed on %s", cdoStreamName(3));

  streamDefVlist(streamID4, vlistID4);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);
  field.ptr = (double *) malloc(gridsize*sizeof(double));

  datetime = (datetime_t *) malloc((ndates+1)*sizeof(datetime_t));
  
  vars1 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));
  
  for ( its = 0; its < ndates; its++ )
    {
      vars1[its] = (field_t **) malloc(nvars*sizeof(field_t *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  nlevels  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  missval  = vlistInqVarMissval(vlistID1, varID);

	  vars1[its][varID] = (field_t *) malloc(nlevels*sizeof(field_t));

	  for ( levelID = 0; levelID < nlevels; levelID++ )
	    {
	      vars1[its][varID][levelID].grid    = gridID;
	      vars1[its][varID][levelID].nmiss   = 0;
	      vars1[its][varID][levelID].missval = missval;
	      vars1[its][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
	    }
	}
    }

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

      cdiDecodeDate(vdate, &year, &month, &day);

      if ( month >= 1 && month <= 12 )
	dayoy = (month-1)*31 + day;
      else
	dayoy = 0;

      if ( dayoy < 0 || dayoy >= NDAY )
	cdoAbort("Day %d out of range!", dayoy);

      vdates2[dayoy] = vdate;
      vtimes2[dayoy] = vtime;

      if ( vars2[dayoy] == NULL )
	{
	  vars2[dayoy] = (field_t **) malloc(nvars*sizeof(field_t *));
          hsets[dayoy] = hsetCreate(nvars);

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID   = vlistInqVarGrid(vlistID2, varID);
	      gridsize = gridInqSize(gridID);
	      nlevels  = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      missval  = vlistInqVarMissval(vlistID2, varID);

	      vars2[dayoy][varID] = (field_t *)  malloc(nlevels*sizeof(field_t));
              hsetCreateVarLevels(hsets[dayoy], varID, nlevels, gridID);
	      
	      for ( levelID = 0; levelID < nlevels; levelID++ )
		{
		  vars2[dayoy][varID][levelID].grid    = gridID;
		  vars2[dayoy][varID][levelID].nmiss   = 0;
		  vars2[dayoy][varID][levelID].missval = missval;
		  vars2[dayoy][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		}
	    }
	}
      
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID2, &varID, &levelID);
	  streamReadRecord(streamID2, vars2[dayoy][varID][levelID].ptr, &nmiss);
          vars2[dayoy][varID][levelID].nmiss = nmiss;
        }
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID3, &varID, &levelID);
	  streamReadRecord(streamID3, field.ptr, &nmiss);
          field.nmiss   = nmiss;
          field.grid    = vars2[dayoy][varID][levelID].grid;
	  field.missval = vars2[dayoy][varID][levelID].missval;
	  
	  hsetDefVarLevelBounds(hsets[dayoy], varID, levelID, &vars2[dayoy][varID][levelID], &field);
        }
      
      tsID++;
    }
  
  for ( tsID = 0; tsID < ndates; tsID++ )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 )
	cdoAbort("File has less then %d timesteps!", ndates);

      datetime[tsID].date = taxisInqVdate(taxisID1);
      datetime[tsID].time = taxisInqVtime(taxisID1);
	
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;
	    }
	  
	  streamReadRecord(streamID1, vars1[tsID][varID][levelID].ptr, &nmiss);
	  vars1[tsID][varID][levelID].nmiss = nmiss;
	}
    }
  
  while ( TRUE )
    {
      datetime_avg(dpy, ndates, datetime);
      
      vdate = datetime[ndates].date;
      vtime = datetime[ndates].time;
      
      cdiDecodeDate(vdate, &year, &month, &day);

      if ( month >= 1 && month <= 12 )
	dayoy = (month-1)*31 + day;
      else
	dayoy = 0;

      if ( dayoy < 0 || dayoy >= NDAY )
	cdoAbort("Day %d out of range!", dayoy);

      vdates1[dayoy] = vdate;
      vtimes1[dayoy] = vtime;
      
      if ( vars2[dayoy] == NULL )
        cdoAbort("No data for day %d in %s and %s", dayoy, cdoStreamName(1), cdoStreamName(2));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	  nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      
	  for ( levelID = 0; levelID < nlevels; levelID++ )
	    for ( inp = 0; inp < ndates; inp ++ )
	      hsetAddVarLevelValues(hsets[dayoy], varID, levelID, &vars1[inp][varID][levelID]);
	}
        
      datetime[ndates] = datetime[0];
      vars1[ndates] = vars1[0];

      for ( inp = 0; inp < ndates; inp++ )
	{
	  datetime[inp] = datetime[inp+1];
	  vars1[inp] = vars1[inp+1];
	}

      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      datetime[ndates-1].date = taxisInqVdate(taxisID1);
      datetime[ndates-1].time = taxisInqVtime(taxisID1);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  
	  streamReadRecord(streamID1, vars1[ndates-1][varID][levelID].ptr, &nmiss);
	  vars1[ndates-1][varID][levelID].nmiss = nmiss;
	}

      nsets[dayoy] += ndates;
      tsID++;
    }

  otsID = 0;
  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( nsets[dayoy] )
      {
        if ( vdates1[dayoy] != vdates2[dayoy] )
          cdoAbort("Verification dates for day %d of %s, %s and %s are different!", dayoy, cdoStreamName(1), cdoStreamName(2), cdoStreamName(3));
        if ( vtimes1[dayoy] != vtimes2[dayoy] )
          cdoAbort("Verification times for day %d of %s, %s and %s are different!", dayoy, cdoStreamName(1), cdoStreamName(2), cdoStreamName(3));

	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      
	    for ( levelID = 0; levelID < nlevels; levelID++ )
	      hsetGetVarLevelPercentiles(&vars2[dayoy][varID][levelID], hsets[dayoy], varID, levelID, pn);
	  }

	taxisDefVdate(taxisID4, vdates1[dayoy]);
	taxisDefVtime(taxisID4, vtimes1[dayoy]);
	streamDefTimestep(streamID4, otsID++);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
	      {
		streamDefRecord(streamID4, varID, levelID);
		streamWriteRecord(streamID4, vars2[dayoy][varID][levelID].ptr,
		  vars2[dayoy][varID][levelID].nmiss);
	      }
	  }
      }
  
  for ( its = 0; its < ndates; its++ )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevels; levelID++ )
	    free(vars1[its][varID][levelID].ptr);
	  free(vars1[its][varID]);
	}
      free(vars1[its]);
    }
  free(vars1);

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    {
      if ( vars2[dayoy] != NULL )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      for ( levelID = 0; levelID < nlevels; levelID++ )
		free(vars2[dayoy][varID][levelID].ptr);
	      free(vars2[dayoy][varID]);
	    }
	  free(vars2[dayoy]); 
	  hsetDestroy(hsets[dayoy]);
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
