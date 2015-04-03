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

      Gridcell   gridarea        Grid cell area in m^2
      Gridcell   gridweights     Grid cell weights
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Gridcell(void *argument)
{
  static char func[] = "Gridcell";
  int GRIDAREA, GRIDWGTS;
  int operatorID;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridID, zaxisID;
  int gridsize, gridtype;
  int status;
  int ngrids;
  int tsID, varID, levelID, taxisID;
  double *grid_area = NULL;
  double *grid_wgts = NULL;
  double *pdata;
  char *envstr;
  double  EarthRadius = 6371000; /* default radius of the earth in m */
  double PlanetRadius = EarthRadius;

  cdoInitialize(argument);

  envstr = getenv("PLANET_RADIUS");
  if ( envstr )
    {
      double fval;
      fval = atof(envstr);
      if ( fval > 0 )
	{
	  PlanetRadius = fval;
	  if ( cdoVerbose )
	    cdoPrint("Set PlanetRadius to %g", PlanetRadius);
	}
    }

  GRIDAREA = cdoOperatorAdd("gridarea",     0,  0, NULL);
  GRIDWGTS = cdoOperatorAdd("gridweights",  0,  0, NULL);

  operatorID = cdoOperatorID();

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  ngrids = vlistNgrids(vlistID1);

  if ( ngrids > 1 )
    cdoWarning("Found more than 1 grid, using the first one!");

  gridID  = vlistGrid(vlistID1, 0);

  zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  vlistID2 = vlistCreate();
  varID    = vlistDefVar(vlistID2, gridID, zaxisID, TIME_CONSTANT);

  if ( operatorID == GRIDAREA )
    {
      vlistDefVarName(vlistID2, varID, "cell_area");
      vlistDefVarStdname(vlistID2, varID, "area");
      vlistDefVarLongname(vlistID2, varID, "area of grid cell");
      vlistDefVarUnits(vlistID2, varID, "m2");
    }
  else if ( operatorID == GRIDWGTS )
    vlistDefVarName(vlistID2, varID, "cell_weights");

  taxisID = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);


  gridsize = gridInqSize(gridID);
  grid_area = (double *) malloc(gridsize*sizeof(double));
  grid_wgts = (double *) malloc(gridsize*sizeof(double));


  tsID = 0;
  streamDefTimestep(streamID2, tsID);

  varID   = 0;
  levelID = 0;
  streamDefRecord(streamID2, varID, levelID);


  if ( operatorID == GRIDAREA )
    {
      gridtype = gridInqType(gridID);
      if ( gridtype != GRID_LONLAT      &&
	   gridtype != GRID_GAUSSIAN    &&
	   gridtype != GRID_LCC     &&
	   gridtype != GRID_GME         &&
	   gridtype != GRID_CURVILINEAR &&
	   gridtype != GRID_CELL )
	{
	  if ( gridInqType(gridID) == GRID_GAUSSIAN_REDUCED )
	    cdoAbort("Gridarea for %s data failed! Use CDO option -R to convert reduced to regular grid!",
		     gridNamePtr(gridtype));
	  else
	    cdoAbort("Gridarea for %s data failed!", gridNamePtr(gridtype));
	}
      else
	{
	  if ( gridHasArea(gridID) )
	    {
	      if ( cdoVerbose ) cdoPrint("Using existing grid cell area!");
	      gridInqArea(gridID, grid_area);
	    }
	  else
	    {
	      int i;

	      status = gridGenArea(gridID, grid_area);
	      if ( status == 1 )
		cdoAbort("Grid corner missing!");
	      else if ( status == 2 )
		cdoAbort("Can't compute grid cell areas for this grid!");

	      for ( i = 0; i < gridsize; ++i )
		grid_area[i] *= PlanetRadius*PlanetRadius;
	    }
	}

      pdata = grid_area;
    }
  else
    {
      status = gridWeights(gridID, grid_wgts);
      if ( status != 0 )
	  cdoWarning("Using constant grid cell area weights!");

      pdata = grid_wgts;
    }

  streamWriteRecord(streamID2, pdata, 0);


  streamClose(streamID2);
  streamClose(streamID1);

  if ( grid_area ) free(grid_area);
  if ( grid_wgts ) free(grid_wgts);

  cdoFinish();

  return (0);
}
