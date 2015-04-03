/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2009 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      writegrid Write grid
*/


#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void *Writegrid(void *argument)
{
  static char func[] = "Writegrid";
  int WRITEGRID, GRIDAREA;
  int streamID;
  int vlistID1;
  int gridID;
  int operatorID;
  int i;
  int gridtype, gridsize;
  int varID, levelID;
  int nrecs;
  int nmiss;
  int *imask = NULL;
  double missval;
  double *array = NULL;

  cdoInitialize(argument);

  WRITEGRID = cdoOperatorAdd("writegrid", 0, 0, NULL);
  GRIDAREA  = cdoOperatorAdd("gridarea",  0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID = streamOpenRead(cdoStreamName(0));
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID);

  if ( operatorID == WRITEGRID )
    {
      nrecs = streamInqTimestep(streamID, 0);

      streamInqRecord(streamID, &varID, &levelID);

      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridtype = gridInqType(gridID);
      gridsize = gridInqSize(gridID);

      if ( gridtype == GRID_GME ) gridID = gridToCell(gridID);

      if ( gridtype != GRID_CURVILINEAR && gridtype != GRID_CELL )
	gridID = gridToCurvilinear(gridID);

      if ( gridInqXbounds(gridID, NULL) == 0 || gridInqYbounds(gridID, NULL) == 0 )
	cdoAbort("Grid corner missing!");

      array = (double *) malloc(gridsize*sizeof(double));
      imask = (int *) malloc(gridsize*sizeof(int));
      streamReadRecord(streamID, array, &nmiss);

      missval = vlistInqVarMissval(vlistID1, varID);
      for ( i = 0; i < gridsize; i++ )
	{
	  if ( DBL_IS_EQUAL(array[i], missval) )
	    imask[i] = 0;
	  else
	    imask[i] = 1;
	}
      
      writeNCgrid(cdoStreamName(1), gridID, imask);
    }
  else if ( operatorID == GRIDAREA )
    {
      double *area, *xvals, *yvals;
      int zaxisID, vlistID2, streamID2;
      int nlon, nlat;

      gridID = vlistGrid(vlistID1, 0);
      gridtype = gridInqType(gridID);
      gridsize = gridInqSize(gridID);
      if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT )
	{
	  area  = (double *) malloc(gridsize*sizeof(double));
	  xvals = (double *) malloc(gridsize*sizeof(double));
	  yvals = (double *) malloc(gridsize*sizeof(double));

	  zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

	  vlistID2 = vlistCreate();
	  varID = vlistDefVar(vlistID2, gridID, zaxisID, TIME_CONSTANT);
	  /*	  vlistDefVarCode(vlistID2, varID, 11); */
	  vlistDefVarName(vlistID2, varID, "area");
	  vlistDefVarLongname(vlistID2, varID, "area");
	  vlistDefVarUnits(vlistID2, varID, "m**2");

	  gridID = gridToCurvilinear(gridID);

	  nlon = gridInqXsize(gridID);
	  nlat = gridInqYsize(gridID);
	  gridInqXvals(gridID, xvals);
	  gridInqYvals(gridID, yvals);

	  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
	  streamDefVlist(streamID2, vlistID2);

	  streamDefRecord(streamID2, varID, 0);
	  nmiss = 0;
	  memcpy(area, xvals, gridsize*sizeof(double));
	  streamWriteRecord(streamID2, area, nmiss);

	  streamClose(streamID2);

	  free(area);
	  free(xvals);
	  free(yvals);
	}
      else
	cdoAbort("%s grid unsupported!", gridNamePtr(gridtype));
    }

  streamClose(streamID);

  if (array) free(array);
  if (imask) free(imask);

  cdoFinish();

  return (0);
}
