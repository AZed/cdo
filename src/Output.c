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

      Output     output          ASCII output
      Output     outputf         Formatted output
      Output     outputint       Integer output
      Output     outputsrv       SERVICE output
      Output     outputext       EXTRA output
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"


void *Output(void *argument)
{
  static char func[] = "Output";
  int OUTPUT, OUTPUTINT, OUTPUTSRV, OUTPUTEXT, OUTPUTF, OUTPUTTS, OUTPUTFLD, OUTPUTARR, OUTPUTXYZ;
  int operatorID;
  int i;
  int indf;
  int varID, recID;
  int gridsize = 0;
  int gridID, zaxisID, code, vdate, vtime;
  int ngrids;
  int nrecs;
  int levelID;
  int tsID, taxisID;
  int streamID = 0;
  int vlistID;
  int nmiss, nout;
  int nlon, nlat;
  int hour, minute, second;
  int year, month, day;
  int nelem = 0;
  int index;
  int ndiffgrids;
  const char *format = NULL;
  double level;
  double *grid_center_lon = NULL, *grid_center_lat = NULL;
  double *array = NULL;
  double xdate;
  double missval;

  cdoInitialize(argument);

  OUTPUT    = cdoOperatorAdd("output",    0, 0, NULL);
  OUTPUTINT = cdoOperatorAdd("outputint", 0, 0, NULL);
  OUTPUTSRV = cdoOperatorAdd("outputsrv", 0, 0, NULL);
  OUTPUTEXT = cdoOperatorAdd("outputext", 0, 0, NULL);
  OUTPUTF   = cdoOperatorAdd("outputf",   0, 0, NULL);
  OUTPUTTS  = cdoOperatorAdd("outputts",  0, 0, NULL);
  OUTPUTFLD = cdoOperatorAdd("outputfld", 0, 0, NULL);
  OUTPUTARR = cdoOperatorAdd("outputarr", 0, 0, NULL);
  OUTPUTXYZ = cdoOperatorAdd("outputxyz", 0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == OUTPUTF )
    {
      operatorInputArg("format and number of elements");
      operatorCheckArgc(2);
      format = operatorArgv()[0];
      nelem  = atoi(operatorArgv()[1]);
    }

  for ( indf = 0; indf < cdoStreamCnt(); indf++ )
    {
      streamID = streamOpenRead(cdoStreamName(indf));
      if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(indf));

      vlistID = streamInqVlist(streamID);

      ngrids = vlistNgrids(vlistID);
      ndiffgrids = 0;
      for ( index = 1; index < ngrids; index++ )
	if ( vlistGrid(vlistID, 0) != vlistGrid(vlistID, index) )
	  ndiffgrids++;

      if ( ndiffgrids > 0 ) cdoAbort("Too many different grids!");

      gridID = vlistGrid(vlistID, 0);
      gridsize = gridInqSize(gridID);

      array = (double *) malloc(gridsize*sizeof(double));

      if ( operatorID == OUTPUTFLD || operatorID == OUTPUTXYZ )
	{
	  char units[128];

	  if ( gridInqType(gridID) == GRID_GME ) gridID = gridToCell(gridID);

	  if ( gridInqType(gridID) != GRID_CELL && gridInqType(gridID) != GRID_CURVILINEAR )
	    gridID = gridToCurvilinear(gridID);

	  grid_center_lon = (double *) malloc(gridsize*sizeof(double));
	  grid_center_lat = (double *) malloc(gridsize*sizeof(double));
	  gridInqXvals(gridID, grid_center_lon);
	  gridInqYvals(gridID, grid_center_lat);

	  /* Convert lat/lon units if required */
	  gridInqYunits(gridID, units);

	  gridToDegree(units, "grid center lon", gridsize, grid_center_lon);
	  gridToDegree(units, "grid center lat", gridsize, grid_center_lat);
	}

      tsID = 0;
      taxisID = vlistInqTaxis(vlistID);
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{
	  vdate = taxisInqVdate(taxisID);
	  vtime = taxisInqVtime(taxisID);

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID, &varID, &levelID);

	      code     = vlistInqVarCode(vlistID, varID);
	      gridID   = vlistInqVarGrid(vlistID, varID);
	      zaxisID  = vlistInqVarZaxis(vlistID, varID);
	      missval  = vlistInqVarMissval(vlistID, varID);
	      gridsize = gridInqSize(gridID);
	      nlon     = gridInqXsize(gridID);
	      nlat     = gridInqYsize(gridID);
	      level    = zaxisInqLevel(zaxisID, levelID);
	      
	      if ( nlon*nlat != gridsize ) { nlon = gridsize; nlat = 1; }

	      streamReadRecord(streamID, array, &nmiss);

	      if ( operatorID == OUTPUTSRV )
		fprintf(stdout, "%4d %8g %8d %4d %8d %8d %d %d\n", code, level, vdate, vtime, nlon, nlat, 0, 0);

	      if ( operatorID == OUTPUTEXT )
		fprintf(stdout, "%8d %4d %8g %8d\n", vdate, code, level, gridsize);
		
	      if ( operatorID == OUTPUTINT )
		{
		  nout = 0;
		  for ( i = 0; i < gridsize; i++ )
		    {
		      if ( nout == 8 )
			{
			  nout = 0;
			  fprintf(stdout, "\n");
			}
		      fprintf(stdout, " %8d", (int) array[i]);
		      nout++;
		    }
		  fprintf(stdout, "\n");
		}
	      else if ( operatorID == OUTPUTF )
		{
		  nout = 0;
		  for ( i = 0; i < gridsize; i++ )
		    {
		      if ( nout == nelem )
			{
			  nout = 0;
			  fprintf(stdout, "\n");
			}
		      fprintf(stdout, format, array[i]);
		      nout++;
		    }
		  fprintf(stdout, "\n");
		}
	      else if ( operatorID == OUTPUTTS )
		{
		  if ( gridsize > 1 )
		    cdoAbort("operator works only with one gridpoint!");

		  decode_date(vdate, &year, &month, &day);
		  decode_time(vtime, &hour, &minute, &second);

		  fprintf(stdout, DATE_FORMAT" "TIME_FORMAT" %12.12g\n",
			  year, month, day, hour, minute, second, array[0]);
		}
	      else if ( operatorID == OUTPUTFLD )
		{
		  decode_time(vtime, &hour, &minute, &second);
		  xdate  = vdate - (vdate/100)*100 + (hour*3600 + minute*60 + second)/86400.;
		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(array[i], missval) )
		      fprintf(stdout, "%g\t%g\t%g\t%g\n", xdate, 
			      grid_center_lat[i], grid_center_lon[i], array[i]);
		}
	      else if ( operatorID == OUTPUTXYZ )
		{
		  if ( tsID == 0 && recID == 0 )
		    {
		      char *fname = "frontplane.xyz";
		      FILE *fp;
		      double fmin = 0;
		      double dx, x0, y0, z0, x, y, z;
		      for ( i = 0; i < gridsize; i++ )
			if ( !DBL_IS_EQUAL(array[i], missval) )
			  {
			    if ( array[i] < fmin ) fmin = array[i];
			    fprintf(stdout, "%g\t%g\t%g\t%g\n",
				    grid_center_lon[i], grid_center_lat[i], array[i], array[i]);
			  }
		      fp = fopen(fname, "w");
		      if ( fp == NULL ) cdoAbort("Open failed on %s", fname);
		      // first front plane
		      dx = (grid_center_lon[1] - grid_center_lon[0]);
		      x0 = grid_center_lon[0]-dx/2;
		      y0 = grid_center_lat[0]-dx/2;
		      z0 = fmin;
		      fprintf(fp, ">\n");
		      for ( i = 0; i < nlon; ++i )
			{
			  x = x0;  y = y0; z = z0;
			  fprintf(fp, "%g %g %g\n", x, y, z);
			  x = x0;  y = y0; z = array[i];
			  fprintf(fp, "%g %g %g\n", x, y, z);
			  x = x0+dx;  y = y0;
			  fprintf(fp, "%g %g %g\n", x, y, z);
			  x0 = x; y0 = y0; z0 = z;
			}
		      x = x0;  y = y0; z = fmin;
		      fprintf(fp, "%g %g %g\n", x, y, z);
		      x = grid_center_lon[0]-dx/2;
		      fprintf(fp, "%g %g %g\n", x, y, z);

		      // second front plane
		      x0 = grid_center_lon[0]-dx/2;
		      y0 = grid_center_lat[0]-dx/2;
		      z0 = fmin;
		      fprintf(fp, ">\n");
		      for ( i = 0; i < nlat; ++i )
			{
			  x = x0;  y = y0; z = z0;
			  fprintf(fp, "%g %g %g\n", x, y, z);
			  x = x0;  y = y0; z = array[i*nlon];
			  fprintf(fp, "%g %g %g\n", x, y, z);
			  x = x0;  y = y0+dx;
			  fprintf(fp, "%g %g %g\n", x, y, z);
			  x0 = x0; y0 = y; z0 = z;
			}
		      x = x0;  y = y0; z = fmin;
		      fprintf(fp, "%g %g %g\n", x, y, z);
		      y = grid_center_lat[0]-dx/2;
		      fprintf(fp, "%g %g %g\n", x, y, z);

		      fclose(fp);
		    }
		}
	      else if ( operatorID == OUTPUTARR )
		{
		  for ( i = 0; i < gridsize; i++ )
		    {
		      fprintf(stdout, "  arr[%d] = %12.6g;\n", i, array[i]);
		      nout++;
		    }
		}
	      else
		{
		  nout = 0;
		  for ( i = 0; i < gridsize; i++ )
		    {
		      if ( nout == 6 )
			{
			  nout = 0;
			  fprintf(stdout, "\n");
			}
		      fprintf(stdout, " %12.6g", array[i]);
		      nout++;
		    }
		  fprintf(stdout, "\n");
		}
	    }
	  tsID++;
	}
      streamClose(streamID);

      if ( array ) free(array);
      if ( grid_center_lon ) free(grid_center_lon);
      if ( grid_center_lat ) free(grid_center_lat);
    }

  cdoFinish();

  return (0);
}
