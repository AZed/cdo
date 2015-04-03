/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2007 Brockmann Consult
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

      Seascount   seascount         Seasonal counts
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


void *Seascount(void *argument)
{
  static char func[] = "Seascount";
  int operatorID;
  int operfunc;
  int gridsize;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int nrecs, nrecords;
  int gridID, varID, levelID, recID;
  int tsID;
  int otsID;
  long nsets;
  int i;
  int year, month, day, seas, seas0 = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  int newseas, oldmon = 0, newmon;
  double missval;
  FIELD **vars1 = NULL;
  FIELD field;
  int season_start;

  cdoInitialize(argument);

  cdoOperatorAdd("seascount", 0, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  season_start = get_season_start();

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);

  field.ptr = (double *) malloc(gridsize*sizeof(double));

  vars1 = (FIELD **) malloc(nvars*sizeof(FIELD *));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval  = vlistInqVarMissval(vlistID1, varID);

      vars1[varID] = (FIELD *)  malloc(nlevel*sizeof(FIELD));

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
      newseas = FALSE;
      while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
        {
          vdate = taxisInqVdate(taxisID1);
          vtime = taxisInqVtime(taxisID1);
	  decode_date(vdate, &year, &month, &day);
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
            cdoAbort("Season %d out of range!", seas+1);

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
      streamDefTimestep(streamID2, otsID++);

      for ( recID = 0; recID < nrecords; recID++ )
        {
          varID   = recVarID[recID];
          levelID = recLevelID[recID];

          if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
            {
              streamDefRecord(streamID2, varID, levelID);
              streamWriteRecord(streamID2, vars1[varID][levelID].ptr,  vars1[varID][levelID].nmiss);
            }
        }

      if ( nrecs == 0 ) break;
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
