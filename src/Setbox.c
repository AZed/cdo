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

      Setbox     setclonlatbox   Set lon/lat box to constant
      Setbox     setcindexbox    Set index box to constant
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void genlonlatbox(double xlon1, double xlon2, double xlat1, double xlat2,
		  int nlon1, int nlat1, double *xvals1, double *yvals1,
		  int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22);

static
void genlonlatgrid(int gridID1, int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22, double *constant)
{
  int nlon1, nlat1;
  double *xvals1, *yvals1;
  double xlon1, xlon2, xlat1, xlat2;

  operatorCheckArgc(5);

  *constant = atof(operatorArgv()[0]);
  xlon1 = atof(operatorArgv()[1]);
  xlon2 = atof(operatorArgv()[2]);
  xlat1 = atof(operatorArgv()[3]);
  xlat2 = atof(operatorArgv()[4]);

  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  xvals1 = (double *) malloc(nlon1*sizeof(double));
  yvals1 = (double *) malloc(nlat1*sizeof(double));

  gridInqXvals(gridID1, xvals1);
  gridInqYvals(gridID1, yvals1);

  /* Convert lat/lon units if required */
  {
    char units[CDI_MAX_NAME];
    gridInqXunits(gridID1, units);
    gridToDegree(units, "grid center lon", nlon1, xvals1);
    gridInqYunits(gridID1, units);
    gridToDegree(units, "grid center lat", nlat1, yvals1);
  }

  genlonlatbox(xlon1, xlon2, xlat1, xlat2,
	       nlon1, nlat1, xvals1, yvals1,
	       lat1, lat2, lon11, lon12, lon21, lon22);

  free(xvals1);
  free(yvals1);
}


static
void genindexbox(int gridID1, int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22, double *constant)
{
  int nlon1, nlat1;
  int temp;

  operatorCheckArgc(5);

  *constant = atof(operatorArgv()[0]);
  *lon11 = atoi(operatorArgv()[1]);
  *lon12 = atoi(operatorArgv()[2]);
  *lat1  = atoi(operatorArgv()[3]);
  *lat2  = atoi(operatorArgv()[4]);

  if ( *lat1 > *lat2 )
    {
      temp = *lat1;
      *lat1 = *lat2;
      *lat2 = temp;
    }

  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  if ( *lat1 < 1 )
    {
      cdoWarning("first latitude index out of range. Set to 1.");
      *lat1 = 1;
    }
  if ( *lat2 > nlat1 )
    {
      cdoWarning("last latitude index out of range. Set to %d.", nlat1);
      *lat2 = nlat1;
    }
  if ( *lon11 < 1 )
    {
      cdoWarning("first longitude index out of range. Set to 1.");
      *lon11 = 1;
    }
  if ( *lon12 > nlon1+1 )
    {
      cdoWarning("last longitude index out of range. Set to %d.", nlon1);
      *lon12 = nlon1;
    }

  (*lon11)--;
  (*lon12)--;
  (*lat1)--;
  (*lat2)--;

  if ( *lon11 > *lon12 )
    {
      *lon21 = *lon11;
      *lon22 = nlon1 - 1;
      *lon11 = 0;
    }
  else
    {
      if ( *lon12 > nlon1-1 )
	{
	  *lon21 = *lon11;
	  *lon22 = nlon1 - 1;
	  *lon11 = 0;
	  *lon12 = 0;
	}
      else
	{
	  *lon21 = 0;
	  *lon22 = -1;
	}
    }
}


static
void setcbox(double constant, double *array, int gridID,
	     int lat1, int lat2, int lon11, int lon12, int lon21, int lon22)
{
  int nlon, nlat;
  int ilat, ilon;

  nlon = gridInqXsize(gridID);
  nlat = gridInqYsize(gridID);

  for ( ilat = 0; ilat < nlat; ilat++ )
    for ( ilon = 0; ilon < nlon; ilon++ )
      if ( (lat1 <= ilat && ilat <= lat2 && 
	    ((lon11 <= ilon && ilon <= lon12) || (lon21 <= ilon && ilon <= lon22))) ) {
	array[nlon*ilat + ilon] = constant;
      }
}


void *Setbox(void *argument)
{
  int SETCLONLATBOX, SETCINDEXBOX;
  int operatorID;
  int streamID1, streamID2;
  int nrecs, nvars;
  int tsID, recID, varID, levelID;
  int gridsize;
  int vlistID1, vlistID2;
  int gridID = -1;
  int index, ngrids, gridtype;
  int nmiss;
  int *vars;
  int i;
  int ndiffgrids;
  int lat1, lat2, lon11, lon12, lon21, lon22;
  double missval;
  double constant;
  double *array;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  SETCLONLATBOX = cdoOperatorAdd("setclonlatbox", 0, 0, "constant, western and eastern longitude and southern and northern latitude");
  SETCINDEXBOX  = cdoOperatorAdd("setcindexbox",  0, 0, "constant, index of first and last longitude and index of first and last latitude");


  operatorID = cdoOperatorID();

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  ngrids = vlistNgrids(vlistID1);
  ndiffgrids = 0;
  for ( index = 1; index < ngrids; index++ )
    if ( vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index))
      ndiffgrids++;

  for ( index = 0; index < ngrids; index++ )
    {
      gridID   = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID);
      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN ) break;
      if ( operatorID == SETCINDEXBOX && gridtype == GRID_CURVILINEAR ) break;
      if ( operatorID == SETCINDEXBOX && gridtype == GRID_GENERIC &&
	   gridInqXsize(gridID) > 0 && gridInqYsize(gridID) > 0 ) break;
    }

  if ( gridInqType(gridID) == GRID_GAUSSIAN_REDUCED )
    cdoAbort("Gaussian reduced grid found. Use option -R to convert it to a regular grid!");

  if ( index == ngrids ) cdoAbort("No regular grid found!");
  if ( ndiffgrids > 0 )  cdoAbort("Too many different grids!");

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorID == SETCLONLATBOX )
    genlonlatgrid(gridID, &lat1, &lat2, &lon11, &lon12, &lon21, &lon22, &constant);
  else
    genindexbox(gridID, &lat1, &lat2, &lon11, &lon12, &lon21, &lon22, &constant);

  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nvars = vlistNvars(vlistID1);
  vars  = (int *) malloc(nvars*sizeof(int));
  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( gridID == vlistInqVarGrid(vlistID1, varID) )
	vars[varID] = TRUE;
      else
	vars[varID] = FALSE;
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = gridInqSize(gridID);
  array = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( vars[varID] )
	    {
	      streamReadRecord(streamID1, array, &nmiss);

	      missval = vlistInqVarMissval(vlistID1, varID);
	      setcbox(constant, array, gridID, lat1, lat2, lon11, lon12, lon21, lon22);

	      nmiss = 0;
	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(array[i], missval) ) nmiss++;

	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, array, nmiss);
	    }
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( vars  ) free(vars);
  if ( array ) free(array);

  cdoFinish();

  return (0);
}
