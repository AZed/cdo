/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2011 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Arith      add             Add two fields
      Arith      sub             Subtract two fields
      Arith      mul             Multiply two fields
      Arith      div             Divide two fields
      Arith      min             Minimum of two fields
      Arith      max             Maximum of two fields
      Arith      atan2           Arc tangent of two fields
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Arith(void *argument)
{
  int operatorID;
  int operfunc;
  enum {FILL_NONE, FILL_TS, FILL_REC, FILL_RECTS, FILL_FILE};
  int filltype = FILL_NONE;
  int streamIDx1, streamIDx2, streamID1, streamID2, streamID3;
  int gridsize;
  int nrecs, nrecs2, nvars = 0, nlev, recID;
  int tsID, tsID2;
  int varID, levelID;
  int offset;
  int ntsteps1, ntsteps2;
  int vlistIDx1, vlistIDx2, vlistID1, vlistID2, vlistID3;
  int taxisIDx1, taxisID1, taxisID2, taxisID3;
  field_t *fieldx1, *fieldx2, fieldrec, field1, field2;
  int **varnmiss = NULL;
  double **vardata = NULL;

  cdoInitialize(argument);

  cdoOperatorAdd("add",   func_add,   0, NULL);
  cdoOperatorAdd("sub",   func_sub,   0, NULL);
  cdoOperatorAdd("mul",   func_mul,   0, NULL);
  cdoOperatorAdd("div",   func_div,   0, NULL);
  cdoOperatorAdd("min",   func_min,   0, NULL);
  cdoOperatorAdd("max",   func_max,   0, NULL);
  cdoOperatorAdd("atan2", func_atan2, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));

  streamIDx1 = streamID1;
  streamIDx2 = streamID2;
  fieldx1 = &field1;
  fieldx2 = &field2;

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistIDx1 = vlistID1;
  vlistIDx2 = vlistID2;

  if ( cdoVerbose ) vlistPrint(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisIDx1 = taxisID1;

  ntsteps1 = vlistNtsteps(vlistID1);
  ntsteps2 = vlistNtsteps(vlistID2);
  if ( ntsteps1 == 0 ) ntsteps1 = 1;
  if ( ntsteps2 == 0 ) ntsteps2 = 1;

  if ( vlistNrecs(vlistID1) != 1 && vlistNrecs(vlistID2) == 1 )
    {
      if ( ntsteps1 != 1 && ntsteps2 == 1 )
	{
	  filltype = FILL_REC;
	  cdoPrint("Filling up stream2 >%s< by copying the first record.", cdoStreamName(1));
	}
      else
	{
	  filltype = FILL_RECTS;
	  cdoPrint("Filling up stream2 >%s< by copying the first record of each timestep.", cdoStreamName(1));
	}
    }
  else if ( vlistNrecs(vlistID1) == 1 && vlistNrecs(vlistID2) != 1 )
    {
      if ( ntsteps1 == 1 && ntsteps2 != 1 )
	{
	  filltype = FILL_REC;
	  cdoPrint("Filling up stream1 >%s< by copying the first record.", cdoStreamName(0));
	}
      else
	{
	  filltype = FILL_RECTS;
	  cdoPrint("Filling up stream1 >%s< by copying the first record of each timestep.", cdoStreamName(0));
	}
      streamIDx1 = streamID2;
      streamIDx2 = streamID1;
      vlistIDx1 = vlistID2;
      vlistIDx2 = vlistID1;
      taxisIDx1 = taxisID2;
      fieldx1 = &field2;
      fieldx2 = &field1;
    }

  if ( filltype == FILL_NONE )
    vlistCompare(vlistID1, vlistID2, CMP_ALL);

  gridsize = vlistGridsizeMax(vlistIDx1);

  field1.ptr = (double *) malloc(gridsize*sizeof(double));
  field2.ptr = (double *) malloc(gridsize*sizeof(double));
  fieldrec.ptr = NULL;
  fieldrec.nmiss = 0;
  if ( filltype == FILL_REC || filltype == FILL_RECTS )
    fieldrec.ptr = (double *) malloc(gridsize*sizeof(double));

  if ( cdoVerbose )
    cdoPrint("Number of timesteps: file1 %d, file2 %d", ntsteps1, ntsteps2);

  if ( filltype == FILL_NONE )
    {
      if ( ntsteps1 != 1 && ntsteps2 == 1 )
	{
	  filltype = FILL_TS;
	  cdoPrint("Filling up stream2 >%s< by copying the first timestep.", cdoStreamName(1));
	}
      else if ( ntsteps1 == 1 && ntsteps2 != 1 )
	{
	  filltype = FILL_TS;
	  cdoPrint("Filling up stream1 >%s< by copying the first timestep.", cdoStreamName(0));
	  streamIDx1 = streamID2;
          streamIDx2 = streamID1;
	  vlistIDx1 = vlistID2;
	  vlistIDx2 = vlistID1;
	  taxisIDx1 = taxisID2;
	  fieldx1 = &field2;
	  fieldx2 = &field1;
	}

      if ( filltype == FILL_TS )
	{
	  nvars  = vlistNvars(vlistIDx2);
	  vardata  = (double **) malloc(nvars*sizeof(double *));
	  varnmiss = (int **) malloc(nvars*sizeof(int *));
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
	      nlev     = zaxisInqSize(vlistInqVarZaxis(vlistIDx2, varID));
	      vardata[varID]  = (double *) malloc(nlev*gridsize*sizeof(double));
	      varnmiss[varID] = (int *) malloc(nlev*sizeof(int));
	    }
	}
    }

  vlistID3 = vlistDuplicate(vlistIDx1);
  if ( filltype == FILL_TS && vlistIDx1 != vlistID1 )
    {
      nvars  = vlistNvars(vlistID1);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarMissval(vlistID3, varID, vlistInqVarMissval(vlistID1, varID));
    }

  taxisID3 = taxisDuplicate(taxisIDx1);
  vlistDefTaxis(vlistID3, taxisID3);

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID3, vlistID3);

  tsID = 0;
  tsID2 = 0;
  while ( (nrecs = streamInqTimestep(streamIDx1, tsID)) )
    {
      if ( tsID == 0 || filltype == FILL_NONE || filltype == FILL_FILE || filltype == FILL_RECTS )
	{
	  nrecs2 = streamInqTimestep(streamIDx2, tsID2);
	  if ( nrecs2 == 0 )
	    {
	      if ( filltype == FILL_NONE && streamIDx2 == streamID2 )
		{
		  filltype = FILL_FILE;
		  cdoPrint("Filling up stream2 >%s< by copying all timesteps.", cdoStreamName(1));
		}

	      if ( filltype == FILL_FILE )
		{
		  tsID2 = 0;
		  streamClose(streamID2);
		  streamID2 = streamOpenRead(cdoStreamName(1));
		  streamIDx2 = streamID2;

		  vlistID2 = streamInqVlist(streamID2);
		  vlistIDx2 = vlistID2;

		  vlistCompare(vlistID1, vlistID2, CMP_DIM);

		  nrecs2 = streamInqTimestep(streamIDx2, tsID2);
		  if ( nrecs2 == 0 )
		    cdoAbort("Empty input stream %s!", cdoStreamName(1));
		}
	      else
		cdoAbort("Input streams have different number of timesteps!");
	    }
	}

      taxisCopyTimestep(taxisID3, taxisIDx1);

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamIDx1, &varID, &levelID);
	  streamReadRecord(streamIDx1, fieldx1->ptr, &fieldx1->nmiss);

	  if ( tsID == 0 || filltype == FILL_NONE || filltype == FILL_FILE || filltype == FILL_RECTS )
	    {
	      if ( recID == 0 || (filltype != FILL_REC && filltype != FILL_RECTS) )
		{
		  streamInqRecord(streamIDx2, &varID, &levelID);
		  streamReadRecord(streamIDx2, fieldx2->ptr, &fieldx2->nmiss);
		}

	      if ( filltype == FILL_TS )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
		  offset   = gridsize*levelID;
		  memcpy(vardata[varID]+offset, fieldx2->ptr, gridsize*sizeof(double));
		  varnmiss[varID][levelID] = fieldx2->nmiss;
		}
	      else if ( recID == 0 && (filltype == FILL_REC || filltype == FILL_RECTS) )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, 0));
		  memcpy(fieldrec.ptr, fieldx2->ptr, gridsize*sizeof(double));
		  fieldrec.nmiss = fieldx2->nmiss;
		}
	    }
	  else if ( filltype == FILL_TS )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
	      offset   = gridsize*levelID;
	      memcpy(fieldx2->ptr, vardata[varID]+offset, gridsize*sizeof(double));
	      fieldx2->nmiss = varnmiss[varID][levelID];
	    }

	  fieldx1->grid    = vlistInqVarGrid(vlistIDx1, varID);
	  fieldx1->missval = vlistInqVarMissval(vlistIDx1, varID);

	  if ( filltype == FILL_REC || filltype == FILL_RECTS )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, 0));
	      memcpy(fieldx2->ptr, fieldrec.ptr, gridsize*sizeof(double));
	      fieldx2->nmiss   = fieldrec.nmiss;
	      fieldx2->grid    = vlistInqVarGrid(vlistIDx2, 0);
	      fieldx2->missval = vlistInqVarMissval(vlistIDx2, 0);
	    }
	  else
	    {
	      fieldx2->grid    = vlistInqVarGrid(vlistIDx2, varID);
	      fieldx2->missval = vlistInqVarMissval(vlistIDx2, varID);
	    }

	  farfun(&field1, field2, operfunc);

	  streamDefRecord(streamID3, varID, levelID);
	  streamWriteRecord(streamID3, field1.ptr, field1.nmiss);
	}

      tsID++;
      tsID2++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  if ( vardata )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  free(vardata[varID]);
	  free(varnmiss[varID]);
	}

      free(vardata);
      free(varnmiss);
    }

  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);
  if ( filltype == FILL_REC || filltype == FILL_RECTS )
    if ( fieldrec.ptr ) free(fieldrec.ptr);

  cdoFinish();

  return (0);
}
