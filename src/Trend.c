/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2010 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Trend      trend           Trend
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Trend(void *argument)
{
  static char func[] = "Trend";
  int gridsize;
  int vdate = 0, vtime = 0;
  int nrecs, nrecords;
  int gridID, varID, levelID, recID;
  int tsID;
  int i, w;
  int streamID1, streamID2, streamID3;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  int nwork = 5;
  double temp1, temp2;
  double missval, missval1, missval2;
  field_t **work[5];
  field_t field1, field2;

  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  vlistDefNtsteps(vlistID2, 1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  for ( varID = 0; varID < nvars; varID++ )
    vlistDefVarDatatype(vlistID2, varID, DATATYPE_FLT64);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));

  streamDefVlist(streamID3, vlistID2);

  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);

  field1.ptr = (double *) malloc(gridsize*sizeof(double));
  field2.ptr = (double *) malloc(gridsize*sizeof(double));

  for ( w = 0; w < nwork; w++ )
    work[w] = (field_t **) malloc(nvars*sizeof(field_t *));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval  = vlistInqVarMissval(vlistID1, varID);

      for ( w = 0; w < nwork; w++ )
	work[w][varID] = (field_t *)  malloc(nlevel*sizeof(field_t));

      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  for ( w = 0; w < nwork; w++ )
	    {
	      work[w][varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));
	      for ( i = 0; i < gridsize; i++ ) work[w][varID][levelID].ptr[i] = 0;
	    }
	}
    }

  tsID    = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      tsID++; /* don't move this line !!! */

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 1 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;
	    }

	  streamReadRecord(streamID1, field1.ptr, &nmiss);

	  missval  = vlistInqVarMissval(vlistID1, varID);
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);

	  for ( i = 0; i < gridsize; i++ )
	    if ( !DBL_IS_EQUAL(field1.ptr[i], missval) )
	      {
		work[0][varID][levelID].ptr[i] += tsID;
		work[1][varID][levelID].ptr[i] += tsID * tsID;
		work[2][varID][levelID].ptr[i] += tsID * field1.ptr[i];
		work[3][varID][levelID].ptr[i] += field1.ptr[i];
		work[4][varID][levelID].ptr[i]++;
	      }      
	}
    }
	  

  taxisDefVdate(taxisID2, vdate);
  taxisDefVtime(taxisID2, vtime);
  streamDefTimestep(streamID2, 0);
  streamDefTimestep(streamID3, 0);

  for ( recID = 0; recID < nrecords; recID++ )
    {
      varID   = recVarID[recID];
      levelID = recLevelID[recID];

      missval  = vlistInqVarMissval(vlistID1, varID);
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);

      missval1  = missval;
      missval2  = missval;

      for ( i = 0; i < gridsize; i++ )
	{
	  temp1 = SUB(work[2][varID][levelID].ptr[i],
		      DIV(MUL(work[0][varID][levelID].ptr[i], work[3][varID][levelID].ptr[i]), work[4][varID][levelID].ptr[i]));
	  temp2 = SUB(work[1][varID][levelID].ptr[i],
		      DIV(MUL(work[0][varID][levelID].ptr[i], work[0][varID][levelID].ptr[i]), work[4][varID][levelID].ptr[i]));

	  field2.ptr[i] = DIV(temp1, temp2);
	  field1.ptr[i] = SUB(DIV(work[3][varID][levelID].ptr[i], work[4][varID][levelID].ptr[i]),
			      MUL(DIV(work[0][varID][levelID].ptr[i], work[4][varID][levelID].ptr[i]), field2.ptr[i]));
	}

      nmiss = 0;
      for ( i = 0; i < gridsize; i++ )
	if ( DBL_IS_EQUAL(field1.ptr[i], missval) ) nmiss++;

      streamDefRecord(streamID2, varID, levelID);
      streamWriteRecord(streamID2, field1.ptr, nmiss);

      nmiss = 0;
      for ( i = 0; i < gridsize; i++ )
	if ( DBL_IS_EQUAL(field2.ptr[i], missval) ) nmiss++;

      streamDefRecord(streamID3, varID, levelID);
      streamWriteRecord(streamID3, field2.ptr, nmiss);
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  for ( w = 0; w < nwork; w++ )
	    free(work[w][varID][levelID].ptr);
	}

      for ( w = 0; w < nwork; w++ )
	free(work[w][varID]);
    }

  for ( w = 0; w < nwork; w++ ) free(work[w]);

  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
