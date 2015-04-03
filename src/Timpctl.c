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

      Timpctl    timpctl         Time percentiles
      Hourpctl   hourpctl        Hourly percentiles
      Daypctl    daypctl         Daily percentiles
      Monpctl    monpctl         Monthly percentiles
      Yearpctl   yearpctl        Yearly percentiles
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "percentiles.h"


static
void timpctl(int operatorID)
{
  int cmplen;
  char indate1[DATE_LEN+1], indate2[DATE_LEN+1];
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
  int streamID1, streamID2, streamID3, streamID4;
  int vlistID1, vlistID2, vlistID3, vlistID4, taxisID1, taxisID2, taxisID3, taxisID4;
  int nmiss;
  int nvars, nlevels;
  int *recVarID, *recLevelID;
  field_t **vars1 = NULL;
  field_t field;
  double pn;
  HISTOGRAM_SET *hset = NULL;
  
  operatorInputArg("percentile number");
  pn = atof(operatorArgv()[0]);
      
  if ( !(pn >= 0 && pn <= 100) )
    cdoAbort("Illegal argument: percentile number %g is not in the range 0..100!", pn);

  cmplen = DATE_LEN - cdoOperatorF2(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  streamID3 = streamOpenRead(cdoStreamName(2));
  
  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = streamInqVlist(streamID3);
  vlistID4 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);
  vlistCompare(vlistID1, vlistID3, CMP_ALL);
  
  if ( cdoOperatorF2(operatorID) == 16 ) vlistDefNtsteps(vlistID4, 1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisID3 = vlistInqTaxis(vlistID3);
  /* TODO - check that time axes 2 and 3 are equal */

  taxisID4 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID4, taxisID4);

  streamID4 = streamOpenWrite(cdoStreamName(3), cdoFiletype());

  streamDefVlist(streamID4, vlistID4);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int*) malloc(nrecords * sizeof(int));
  recLevelID = (int*) malloc(nrecords * sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);

  field_init(&field);
  field.ptr = (double*) malloc(gridsize * sizeof(double));

  vars1 = field_malloc(vlistID1, FIELD_PTR);
  hset = hsetCreate(nvars);
  
  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      nlevels  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

      hsetCreateVarLevels(hset, varID, nlevels, gridID);
    }

  tsID    = 0;
  otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      
      nrecs = streamInqTimestep(streamID2, otsID);
      if ( nrecs != streamInqTimestep(streamID3, otsID) )
        cdoAbort("Number of records at time step %d of %s and %s differ!", otsID+1, cdoStreamName(1)->args, cdoStreamName(2)->args);

      vdate2 = taxisInqVdate(taxisID2);
      vtime2 = taxisInqVtime(taxisID2);
      vdate3 = taxisInqVdate(taxisID3);
      vtime3 = taxisInqVtime(taxisID3);
      if ( vdate2 != vdate3 || vtime2 != vtime3 )
        cdoAbort("Verification dates at time step %d of %s and %s differ!", otsID+1, cdoStreamName(1)->args, cdoStreamName(2)->args);
      
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


	  if ( nsets == 0 ) SET_DATE(indate2, vdate1, vtime1);
	  SET_DATE(indate1, vdate1, vtime1);

	  if ( DATE_IS_NEQ(indate1, indate2, cmplen) ) break;

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
        cdoAbort("Verification dates at time step %d of %s, %s and %s differ!",
		 otsID+1, cdoStreamName(1)->args, cdoStreamName(2)->args, cdoStreamName(3)->args);
      if ( vtime2 != vtime4 )
        cdoAbort("Verification times at time step %d of %s, %s and %s differ!",
		 otsID+1, cdoStreamName(1)->args, cdoStreamName(2)->args, cdoStreamName(3)->args);
      
      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	  nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  
	  for ( levelID = 0; levelID < nlevels; levelID++ )
            hsetGetVarLevelPercentiles(&vars1[varID][levelID], hset, varID, levelID, pn);
	}

      taxisDefVdate(taxisID4, vdate4);
      taxisDefVtime(taxisID4, vtime4);
      streamDefTimestep(streamID4, otsID);

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  varID   = recVarID[recID];
	  levelID = recLevelID[recID];

	  if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	  streamDefRecord(streamID4, varID, levelID);
	  streamWriteRecord(streamID4, vars1[varID][levelID].ptr, vars1[varID][levelID].nmiss);
	}

      if ( nrecs == 0 ) break;
      otsID++;
    }

  field_free(vars1, vlistID1);
  hsetDestroy(hset);
  
  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID4);
  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);
}

void *Timpctl(void *argument)
{
  int operatorID;
  
  cdoInitialize(argument);

  cdoOperatorAdd("timpctl",  func_pctl, 31, NULL);
  cdoOperatorAdd("yearpctl", func_pctl, 10, NULL);
  cdoOperatorAdd("monpctl",  func_pctl,  8, NULL);
  cdoOperatorAdd("daypctl",  func_pctl,  6, NULL);
  cdoOperatorAdd("hourpctl", func_pctl,  4, NULL);

  operatorID = cdoOperatorID();
  
  timpctl(operatorID);

  cdoFinish();

  return (0);
}
