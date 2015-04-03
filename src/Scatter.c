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

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

static
void genGrids(int gridID1, int *gridIDs, int nxvals, int nyvals, int nsx, int nsy,
	      int **gridindex, int nsplit, int gridsize2)
{
  static char *func = "genGrids";
  int gridID2;
  int gridtype;
  int gridsize, nx, ny;
  int index, i, j, ix, iy, ival, offset;
  double *xvals = NULL, *yvals = NULL;

  gridtype = gridInqType(gridID1);
  if ( !(gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GENERIC) )
    cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

  gridsize = gridInqSize(gridID1);
  nx = gridInqXsize(gridID1);
  ny = gridInqYsize(gridID1);

  xvals = (double *) malloc(nx*sizeof(double));
  yvals = (double *) malloc(ny*sizeof(double));

  gridInqXvals(gridID1, xvals);
  gridInqYvals(gridID1, yvals);

  if ( gridsize2 != nxvals*nyvals )
    cdoAbort("Internal problem, gridsize2 differ!");

  index = 0;
  for ( iy = 0; iy < nsy; ++iy )
    for ( ix = 0; ix < nsx; ++ix )
      {
	offset = iy*nyvals*nx + ix*nxvals;
	ival = 0;
	// printf("iy %d, ix %d offset %d\n", iy, ix,  offset);
	for ( j = 0; j < nyvals; ++j )
	  {
	    for ( i = 0; i < nxvals; ++i )
	      {
		//	printf(">> %d %d %d\n", j, i, offset + j*nx + i);
		gridindex[index][ival++] = offset + j*nx + i;
	      }
	  }

	gridID2 = gridCreate(gridtype, gridsize2);
	gridDefXsize(gridID2, nxvals);
	gridDefYsize(gridID2, nyvals);
	gridDefXvals(gridID2, xvals+ix*nxvals);
	gridDefYvals(gridID2, yvals+iy*nyvals);

	gridIDs[index] = gridID2;

	index++;
	if ( index > nsplit )
	  cdoAbort("Internal problem, index exceeded bounds!");
      }

  free(xvals);
  free(yvals);
}

static
void window_cell(double *array1, int gridID1, double *array2, int gridsize2, int *cellidx)
{
  int i;

  for ( i = 0; i < gridsize2; ++i )
    array2[i] = array1[cellidx[i]];
}

typedef struct
{
  int gridID;
  int *gridIDs;
  int **gridindex;
} sgrid_t;


void *Scatter(void *argument)
{
  static char func[] = "Scatter";
  int nchars;
  int streamID1;
  int *vlistIDs = NULL, *streamIDs = NULL;
  int gridID1 = -1, varID;
  int nrecs, ngrids;
  int tsID, recID, levelID;
  int vlistID1;
  char filesuffix[32];
  char filename[4096];
  int index;
  int nsplit, nsx, nsy, ix, iy;
  int xinc, yinc;
  int gridsize, gridsize2;
  int gridtype = -1;
  int nmiss;
  int nx, ny, i;
  double *array1 = NULL, *array2 = NULL;
  sgrid_t *grids;

  cdoInitialize(argument);

  operatorInputArg("xinc, yinc");
  operatorCheckArgc(2);
  xinc = atoi(operatorArgv()[0]);
  yinc = atoi(operatorArgv()[1]);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  ngrids = vlistNgrids(vlistID1);

  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);
      if ( gridtype == GRID_LONLAT   || gridtype == GRID_GAUSSIAN ||
	   (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0) )
	   break;
    }

  if ( index == ngrids )
    cdoAbort("No Lon/Lat, Gaussian or Generic grid found (%s data unsupported)!",
	     gridNamePtr(gridtype));

  gridID1 = vlistGrid(vlistID1, 0);
  gridsize = gridInqSize(gridID1);
  nx = gridInqXsize(gridID1);
  ny = gridInqYsize(gridID1);
  for ( i = 1; i < ngrids; i++ )
    {
      gridID1 = vlistGrid(vlistID1, i);
      if ( gridsize != gridInqSize(gridID1) )
	cdoAbort("Gridsize must not change!");
    }

  if ( nx%xinc != 0 ) cdoAbort("xsize is not multiple of xinc!");
  nsx = nx/xinc;
  if ( ny%yinc != 0 ) cdoAbort("ysize is not multiple of yinc!");
  nsy = ny/yinc;

  nsplit = nsx*nsy;
  gridsize2 = gridsize/nsplit;

  array1 = (double *) malloc(gridsize*sizeof(double));
  array2 = (double *) malloc(gridsize2*sizeof(double));
  
  vlistIDs  = (int *) malloc(nsplit*sizeof(int));
  streamIDs = (int *) malloc(nsplit*sizeof(int));

  grids = (sgrid_t *) malloc(ngrids*sizeof(sgrid_t));
  for ( i = 0; i < ngrids; i++ )
    {  
      gridID1 = vlistGrid(vlistID1, i);
      grids[i].gridID = vlistGrid(vlistID1, i);
      grids[i].gridIDs = (int *) malloc(nsplit*sizeof(int));
      grids[i].gridindex = (int **) malloc(nsplit*sizeof(int*));
      for ( index = 0; index < nsplit; index++ )
	grids[i].gridindex[index] = (int *) malloc(gridsize2*sizeof(int));
    }

  for ( index = 0; index < nsplit; index++ )
    vlistIDs[index] = vlistDuplicate(vlistID1);

  for ( i = 0; i < ngrids; i++ )
    {
      gridID1 = vlistGrid(vlistID1, i);
      genGrids(gridID1, grids[i].gridIDs, xinc, yinc, nsx, nsy, grids[i].gridindex, nsplit, gridsize2);

      for ( index = 0; index < nsplit; index++ )
	vlistChangeGridIndex(vlistIDs[index], i, grids[i].gridIDs[index]);
    }


  strcpy(filename, cdoStreamName(1));
  nchars = strlen(filename);

  filesuffix[0] = 0;
  if ( cdoDisableFilesuffix == FALSE )
    {
      strcat(filesuffix, streamFilesuffix(cdoDefaultFileType));
      if ( cdoDefaultFileType == FILETYPE_GRB )
	if ( vlistIsSzipped(vlistID1) || cdoZtype == COMPRESS_SZIP )
	  strcat(filesuffix, ".sz");
    }


  for ( index = 0; index < nsplit; index++ )
    {
      iy = index/nsx;
      ix = index - iy*nsx;
      //printf("index %d, ix %d, iy %d\n", index, ix, iy);

      sprintf(filename+nchars, "%05d", index);
      if ( filesuffix[0] )
	sprintf(filename+nchars+5, "%s", filesuffix);
      streamIDs[index] = streamOpenWrite(filename, cdoFiletype());
      if ( streamIDs[index] < 0 ) cdiError(streamIDs[index], "Open failed on %s", filename);

      streamDefVlist(streamIDs[index], vlistIDs[index]);
    }

  printf("Bausstelle: i=0 + nmiss!\n");
  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      for ( index = 0; index < nsplit; index++ )
	streamDefTimestep(streamIDs[index], tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  for ( index = 0; index < nsplit; index++ )
	    {
	      i = 0;
	      window_cell(array1, gridID1, array2, gridsize2, grids[i].gridindex[index]);
	      streamDefRecord(streamIDs[index], varID, levelID);
	      streamWriteRecord(streamIDs[index], array2, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID1);

  for ( index = 0; index < nsplit; index++ )
    {
      streamClose(streamIDs[index]);
      vlistDestroy(vlistIDs[index]);
    }

  if ( array1 ) free(array1);
  if ( array2 ) free(array2);

  if ( vlistIDs  ) free(vlistIDs);
  if ( streamIDs ) free(streamIDs);

  for ( i = 0; i < ngrids; i++ )
    {
      for ( index = 0; index < nsplit; index++ )
	gridDestroy(grids[i].gridIDs[index]);
      free(grids[i].gridIDs);

      for ( index = 0; index < nsplit; index++ )
	free(grids[i].gridindex[index]);
      free(grids[i].gridindex);
    }
  free(grids);

  cdoFinish();

  return (0);
}
