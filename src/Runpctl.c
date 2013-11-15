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

      Runpctl    runpctl         Running percentiles
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "nth_element.h"


void *Runpctl(void *argument)
{
  int gridsize;
  int varID;
  int recID;
  int gridID;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  int i, j, inp, its, ndates = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int nmiss;
  int nvars, nlevels;
  int *recVarID, *recLevelID;
  double missval, val;
  field_t ***vars1 = NULL;
  dtinfo_t *dtinfo;
  int taxisID1, taxisID2;
  int calendar;
  int pn;
  double *array;

  cdoInitialize(argument);

  cdoOperatorAdd("runpctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number, number of timesteps");
  operatorCheckArgc(2);
  pn     = atoi(operatorArgv()[0]);
  ndates = atoi(operatorArgv()[1]);

  if ( pn < 1 || pn > 99 )
    cdoAbort("Illegal argument: percentile number %d is not in the range 1..99!", pn);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  calendar = taxisInqCalendar(taxisID1);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  dtinfo = (dtinfo_t *) malloc((ndates+1)*sizeof(dtinfo_t));
  vars1 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));
  array = (double *) malloc(ndates*sizeof(double));
  
  for ( its = 0; its < ndates; its++ )
    {
      vars1[its] = field_malloc(vlistID1, FIELD_PTR);
    }

  for ( tsID = 0; tsID < ndates; tsID++ )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 )
        cdoAbort("File has less than %d timesteps!", ndates);

      taxisInqDTinfo(taxisID1, &dtinfo[tsID]);
        
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

  otsID = 0;
  while ( TRUE )
    {
      for ( varID = 0; varID < nvars; varID++ )
        {
          if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
          
          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          missval  = vlistInqVarMissval(vlistID1, varID);
          nlevels  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          
          for ( levelID = 0; levelID < nlevels; levelID++ )
            {
              nmiss = 0;  
              for ( i = 0; i < gridsize; i++ )
                {
                  for ( inp = 0, j = 0; inp < ndates; inp++ )
                    {
                      val = vars1[inp][varID][levelID].ptr[i];
                      if ( !DBL_IS_EQUAL(val, missval) )
                        array[j++] = val;
                    }
                  
                  if ( j > 0 )
                    {
                      vars1[0][varID][levelID].ptr[i] = nth_element(array, j, (int)ceil(j*(pn/100.0))-1);
                    }
                  else
                    {
                      vars1[0][varID][levelID].ptr[i] = missval;
                      nmiss++;
                    }
                }
              vars1[0][varID][levelID].nmiss = nmiss;  
            }
        }
     
      datetime_avg_dtinfo(calendar, ndates, dtinfo);

      if ( taxisHasBounds(taxisID2) )
	{
	  dtinfo[ndates].b[0] = dtinfo[0].b[0];
	  dtinfo[ndates].b[1] = dtinfo[ndates-1].b[1];
	}

      taxisDefDTinfo(taxisID2, dtinfo[ndates]);
      streamDefTimestep(streamID2, otsID);

      for ( recID = 0; recID < nrecords; recID++ )
        {
          varID    = recVarID[recID];
          levelID  = recLevelID[recID];

	  if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, vars1[0][varID][levelID].ptr, vars1[0][varID][levelID].nmiss);
        }

      otsID++;

      dtinfo[ndates] = dtinfo[0];
      vars1[ndates] = vars1[0];

      for ( inp = 0; inp < ndates; inp++ )
        {
          dtinfo[inp] = dtinfo[inp+1];
          vars1[inp] = vars1[inp+1];
        }

      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      taxisInqDTinfo(taxisID1, &dtinfo[ndates-1]);

      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID1, &varID, &levelID);
          
          streamReadRecord(streamID1, vars1[ndates-1][varID][levelID].ptr, &nmiss);
          vars1[ndates-1][varID][levelID].nmiss = nmiss;
        }

      tsID++;
    }

  for ( its = 0; its < ndates; its++ )
    {
      field_free(vars1[its], vlistID1);
    }

  free(vars1);
  free(array);
  
  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
