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

#include <stdio.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "functs.h"
#include "field.h"
#include "dmemory.h"
#include "nth_element.h"


void *Runpctl(void *argument)
{
  static char func[] = "Runpctl";
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
  FIELD ***vars1 = NULL;
  datetime_t *datetime;
  int taxisID1, taxisID2;
  int calendar, dpy;
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
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  calendar = taxisInqCalendar(taxisID1);
  dpy      = calendar_dpy(calendar);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  datetime = (datetime_t *) malloc((ndates+1)*sizeof(datetime_t));
  vars1 = (FIELD ***) malloc((ndates+1)*sizeof(FIELD **));
  array = (double *) malloc(ndates*sizeof(double));
  
  for ( its = 0; its < ndates; its++ )
    {
      vars1[its] = (FIELD **) malloc(nvars*sizeof(FIELD *));

      for ( varID = 0; varID < nvars; varID++ )
        {
          gridID   = vlistInqVarGrid(vlistID1, varID);
          gridsize = gridInqSize(gridID);
          nlevels  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          missval  = vlistInqVarMissval(vlistID1, varID);

          vars1[its][varID] = (FIELD *)  malloc(nlevels*sizeof(FIELD));

          for ( levelID = 0; levelID < nlevels; levelID++ )
            {
              vars1[its][varID][levelID].grid    = gridID;
              vars1[its][varID][levelID].nmiss   = 0;
              vars1[its][varID][levelID].missval = missval;
              vars1[its][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
            }
        }
    }

  for ( tsID = 0; tsID < ndates; tsID++ )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 )
        cdoAbort("File has less than %d timesteps!", ndates);

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

  otsID = 0;
  while ( TRUE )
    {
      for ( varID = 0; varID < nvars; varID++ )
        {
          if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
          
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
     
      datetime_avg(dpy, ndates, datetime);

      taxisDefVdate(taxisID2, datetime[ndates].date);
      taxisDefVtime(taxisID2, datetime[ndates].time);
      streamDefTimestep(streamID2, otsID++);

      for ( recID = 0; recID < nrecords; recID++ )
        {
          varID    = recVarID[recID];
          levelID  = recLevelID[recID];

          if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
            {
              streamDefRecord(streamID2, varID, levelID);
              streamWriteRecord(streamID2, vars1[0][varID][levelID].ptr, vars1[0][varID][levelID].nmiss);
            }
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

      tsID++;
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
  free(array);
  
  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
