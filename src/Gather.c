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

#if defined (_OPENMP)
#  include <omp.h>
#endif

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


typedef struct
{
  int streamID;
  int vlistID;
  int gridID;
  double *array;
} ens_file_t;


typedef struct
{
  double x, y;
  int id;
}
xyinfo_t;

static
int cmpx(const void *s1, const void *s2)
{
  int cmp = 0;
  xyinfo_t *xy1 = (xyinfo_t *) s1;
  xyinfo_t *xy2 = (xyinfo_t *) s2;

  if      ( xy1->x < xy2->x ) cmp = -1;
  else if ( xy1->x > xy2->x ) cmp =  1;

  return (cmp);
}

static
int cmpxy(const void *s1, const void *s2)
{
  int cmp = 0;
  xyinfo_t *xy1 = (xyinfo_t *) s1;
  xyinfo_t *xy2 = (xyinfo_t *) s2;

  if      ( xy1->y < xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x < xy2->x) ) cmp = -1;
  else if ( xy1->y > xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x > xy2->x) ) cmp =  1;

  return (cmp);
}

static
int genGrid(int nfiles, ens_file_t *ef, int **gridindex, int igrid)
{
  static char *func = "genGrid";
  int fileID;
  int gridID;
  int gridID2 = -1;
  int gridtype = -1;
  int xsize = 0, ysize = 0;
  int xsize2, ysize2;
  int nx, ny, ix, iy, i, j, ij, offset;
  double **xvals, **yvals;
  double *xvals2, *yvals2;
  xyinfo_t *xyinfo;

  xyinfo = (xyinfo_t *) malloc(nfiles*sizeof(xyinfo_t));
  xvals = (double **) malloc(nfiles*sizeof(double));
  yvals = (double **) malloc(nfiles*sizeof(double));

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      gridID = vlistGrid(ef[fileID].vlistID, igrid);
      gridtype = gridInqType(gridID);
      if ( !(gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN ||
	     (gridtype == GRID_GENERIC && gridInqXsize(gridID) > 0 && gridInqYsize(gridID) > 0)) )
	cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

      if ( xsize == 0 ) xsize = gridInqXsize(gridID);
      if ( ysize == 0 ) ysize = gridInqYsize(gridID);
      if ( xsize != gridInqXsize(gridID) ) cdoAbort("xsize differ!");
      if ( ysize != gridInqYsize(gridID) ) cdoAbort("ysize differ!");

      xvals[fileID] = (double *) malloc(xsize*sizeof(double));
      yvals[fileID] = (double *) malloc(ysize*sizeof(double));
      gridInqXvals(gridID, xvals[fileID]);
      gridInqYvals(gridID, yvals[fileID]);

      // printf("fileID %d, gridID %d\n", fileID, gridID);

      xyinfo[fileID].x  = xvals[fileID][0];
      xyinfo[fileID].y  = yvals[fileID][0];
      xyinfo[fileID].id = fileID;
    }
  /*
  for ( fileID = 0; fileID < nfiles; fileID++ )
    printf("1 %d %g %g \n",  xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);
  */
  qsort(xyinfo, nfiles, sizeof(xyinfo_t), cmpx);  	      
  /*
  for ( fileID = 0; fileID < nfiles; fileID++ )
    printf("2 %d %g %g \n",  xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);
  */
  qsort(xyinfo, nfiles, sizeof(xyinfo_t), cmpxy);  	      
  /*
  for ( fileID = 0; fileID < nfiles; fileID++ )
    printf("3 %d %g %g \n",  xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);
  */
  nx = 1;
  for ( fileID = 1; fileID < nfiles; fileID++ )
    {
      if ( DBL_IS_EQUAL(xyinfo[0].y, xyinfo[fileID].y) ) nx++;
      else break;
    }
  ny = nfiles/nx;
  // printf("nx %d  ny %d\n", nx, ny);
  xsize2 = nx*xsize;
  ysize2 = ny*ysize;

  xvals2 = (double *) malloc(xsize2*sizeof(double));
  yvals2 = (double *) malloc(ysize2*sizeof(double));

  for ( i = 0; i < nx; ++i )
    memcpy(xvals2+i*xsize, xvals[xyinfo[i].id], xsize*sizeof(double));

  for ( j = 0; j < ny; ++j )
    memcpy(yvals2+j*ysize, yvals[xyinfo[j*nx].id], ysize*sizeof(double));

  if ( igrid == 0 )
    {
      for ( fileID = 0; fileID < nfiles; fileID++ )
	{
	  iy = fileID/nx;
	  ix = fileID - iy*nx;

          offset = iy*xsize2*ysize + ix*xsize;
	  /*
	  printf("fileID %d %d, iy %d, ix %d, offset %d\n",
		 fileID, xyinfo[fileID].id, iy, ix, offset);
	  */
	  ij = 0;
	  for ( j = 0; j < ysize; ++j )
	    for ( i = 0; i < xsize; ++i )
	      {
		gridindex[xyinfo[fileID].id][ij++] = offset+j*xsize2+i;
	      }
	}
    }

  gridID2 = gridCreate(gridtype, xsize2*ysize2);
  gridDefXsize(gridID2, xsize2);
  gridDefYsize(gridID2, ysize2);
  gridDefXvals(gridID2, xvals2);
  gridDefYvals(gridID2, yvals2);

  free(xvals2);
  free(yvals2);

  {
    char string[1024];
    string[0] = 0;
    gridID = vlistGrid(ef[0].vlistID, igrid);
    gridInqXname(gridID, string);
    gridDefXname(gridID2, string);
    gridInqYname(gridID, string);
    gridDefYname(gridID2, string);
    gridInqXlongname(gridID, string);
    gridDefXlongname(gridID2, string);
    gridInqYlongname(gridID, string);
    gridDefYlongname(gridID2, string);
    gridInqXunits(gridID, string);
    gridDefXunits(gridID2, string);
    gridInqYunits(gridID, string);
    gridDefYunits(gridID2, string);
  }

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      free(xvals[fileID]);
      free(yvals[fileID]);
    }
  free(xvals);
  free(yvals);
  free(xyinfo);

  return (gridID2);
}


void *Gather(void *argument)
{
  static char func[] = "Gather";
  int i;
  int nvars;
  int cmpfunc;
  int varID, recID;
  int gridsize = 0;
  int gridsize2;
  int gridID2;
  int *gridIDs = NULL;
  int nrecs, nrecs0;
  int ngrids;
  int levelID;
  int tsID;
  int streamID = 0, streamID2;
  int vlistID, vlistID1, vlistID2;
  int nmiss;
  int taxisID1, taxisID2;
  int **gridindex;
  double missval;
  double *array2 = NULL;
  int fileID, nfiles;
  const char *ofilename;
  ens_file_t *ef = NULL;

  cdoInitialize(argument);
    
  nfiles = cdoStreamCnt() - 1;

  ofilename = cdoStreamName(nfiles);

  if ( !cdoSilentMode )
    if ( fileExist(ofilename) )
      if ( !userFileOverwrite(ofilename) )
	cdoAbort("Outputfile %s already exist!", ofilename);

  ef = (ens_file_t *) malloc(nfiles*sizeof(ens_file_t));

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID = streamOpenRead(cdoStreamName(fileID));
      if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(fileID));

      vlistID = streamInqVlist(streamID);

      ef[fileID].streamID = streamID;
      ef[fileID].vlistID  = vlistID;
    }

  /* check that the contents is always the same */
  nvars = vlistNvars(ef[0].vlistID);
  if ( nvars == 1 ) 
    cmpfunc = func_sftn;
  else
    cmpfunc = func_sftn;

  for ( fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(ef[0].vlistID, ef[fileID].vlistID, cmpfunc);

  vlistID1 = ef[0].vlistID;
  gridsize = vlistGridsizeMax(vlistID1);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    ef[fileID].array = (double *) malloc(gridsize*sizeof(double));

  ngrids = vlistNgrids(ef[0].vlistID);
  gridIDs = (int *) malloc(ngrids*sizeof(int));
  gridindex = (int **) malloc(nfiles*sizeof(int*));
  for ( fileID = 0; fileID < nfiles; fileID++ )
    gridindex[fileID] = (int *) malloc(gridsize*sizeof(int));

  for ( i = 0; i < ngrids; ++i )
    gridIDs[i] = genGrid(nfiles, ef, gridindex, i);


  vlistID2 = vlistDuplicate(vlistID1);
  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  gridID2 = gridIDs[0];
  gridsize2 = gridInqSize(gridID2);
  // printf("gridsize2 %d\n", gridsize2);
  for ( i = 1; i < ngrids; ++i )
    {
      if ( gridsize2 != gridInqSize(gridIDs[i]) )
	cdoAbort("gridsize differ!");
    }
  for ( i = 0; i < ngrids; ++i )
    {
      vlistChangeGridIndex(vlistID2, i, gridIDs[i]);
    }

  streamID2 = streamOpenWrite(ofilename, cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", ofilename);

  streamDefVlist(streamID2, vlistID2);
	  
  array2 = (double *) malloc(gridsize2*sizeof(double));

  tsID = 0;
  do
    {
      nrecs0 = streamInqTimestep(ef[0].streamID, tsID);
      for ( fileID = 1; fileID < nfiles; fileID++ )
	{
	  streamID = ef[fileID].streamID;
	  nrecs = streamInqTimestep(streamID, tsID);
	  if ( nrecs != nrecs0 )
	    cdoAbort("Number of records changed from %d to %d", nrecs0, nrecs);
	}

      taxisCopyTimestep(taxisID2, taxisID1);

      if ( nrecs0 > 0 ) streamDefTimestep(streamID2, tsID);
      
      for ( recID = 0; recID < nrecs0; recID++ )
	{
	  missval = vlistInqVarMissval(vlistID1, 0);
	  for ( i = 0; i < gridsize2; i++ ) array2[i] = missval;

#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(fileID, streamID, nmiss, i) \
                                     lastprivate(varID, levelID)
#endif
	  for ( fileID = 0; fileID < nfiles; fileID++ )
	    {
	      streamID = ef[fileID].streamID;
	      streamInqRecord(streamID, &varID, &levelID);
	      streamReadRecord(streamID, ef[fileID].array, &nmiss);

	      for ( i = 0; i < gridsize; ++i )
		array2[gridindex[fileID][i]] = ef[fileID].array[i];
	    }

	  nmiss = 0;
	  for ( i = 0; i < gridsize2; i++ )
	    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss);
	}

      tsID++;
    }
  while ( nrecs0 > 0 );

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID = ef[fileID].streamID;
      streamClose(streamID);
    }

  streamClose(streamID2);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    if ( ef[fileID].array ) free(ef[fileID].array);

  if ( ef ) free(ef);
  if ( array2 ) free(array2);

  free(gridIDs);

  cdoFinish();

  return (0);
}
