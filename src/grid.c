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

#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdarg.h> /* va_list */

#if defined (HAVE_LIBPROJ)
#  include "projects.h"
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "error.h"
#include "grid.h"


void gridToDegree(const char *units, const char *string, int gridsize, double *array)
{
  long i;

  if ( memcmp(units, "radian", 6) == 0 )
    {
      for ( i = 0; i < gridsize; i++ ) array[i] *= rad2deg;
    }
  else if ( memcmp(units, "degree", 6) == 0 )
    {
      /* No conversion necessary */
    }
  else
    {
      cdoWarning("Unknown units supplied for %s: %s", string, "proceeding assuming degrees!");
    }
}


int gridToZonal(int gridID1)
{
  int gridID2;
  int gridtype, gridsize;
  double  xval = 0;
  double *yvals;

  gridtype = gridInqType(gridID1);
  gridsize = gridInqYsize(gridID1);
  gridID2  = gridCreate(gridtype, gridsize);
	  
  if ( gridtype == GRID_LONLAT   ||
       gridtype == GRID_GAUSSIAN ||
       gridtype == GRID_GENERIC )
    {
      gridDefXsize(gridID2, 1);
      gridDefYsize(gridID2, gridsize);

      gridDefXvals(gridID2, &xval);

      if ( gridInqYvals(gridID1, NULL) )
	{
	  yvals = (double *) malloc(gridsize*sizeof(double));

	  gridInqYvals(gridID1, yvals);
	  gridDefYvals(gridID2, yvals);

	  free(yvals);
	}
    }
  else
    {
      Error("Gridtype %s unsupported!", gridNamePtr(gridtype));
    }

  return (gridID2);
}


int gridToMeridional(int gridID1)
{
  int gridID2;
  int gridtype, gridsize;
  double *xvals;
  double  yval = 0;

  gridtype = gridInqType(gridID1);
  gridsize = gridInqXsize(gridID1);
  gridID2  = gridCreate(gridtype, gridsize);
	  
  if ( gridtype == GRID_LONLAT   ||
       gridtype == GRID_GAUSSIAN ||
       gridtype == GRID_GENERIC )
    {
      gridDefXsize(gridID2, gridsize);
      gridDefYsize(gridID2, 1);

      if ( gridInqXvals(gridID1, NULL) )
	{
	  xvals = (double *) malloc(gridsize*sizeof(double));

	  gridInqXvals(gridID1, xvals);
	  gridDefXvals(gridID2, xvals);

	  free(xvals);
	}

      gridDefYvals(gridID2, &yval);
    }
  else
    {
      Error("Gridtype %s unsupported!", gridNamePtr(gridtype));
    }

  return (gridID2);
}


void gridGenXbounds(int nx, double *xvals, double *xbounds)
{
  int i;

  for ( i = 0; i < nx-1; i++ )
    {
      xbounds[2*i+1]   = 0.5*(xvals[i] + xvals[i+1]);
      xbounds[2*(i+1)] = 0.5*(xvals[i] + xvals[i+1]);
    }

  xbounds[0]      = 2*xvals[0] - xbounds[1];
  xbounds[2*nx-1] = 2*xvals[nx-1] - xbounds[2*(nx-1)];
}


void gridGenYbounds(int ny, double *yvals, double *ybounds)
{
  int i;

  for ( i = 0; i < ny-1; i++ )
    {
      ybounds[2*i+1]   = 0.5*(yvals[i] + yvals[i+1]);
      ybounds[2*(i+1)] = 0.5*(yvals[i] + yvals[i+1]);
    }

  ybounds[0]      = 2*yvals[0] - ybounds[1];
  ybounds[2*ny-1] = 2*yvals[ny-1] - ybounds[2*(ny-1)];

  if ( yvals[0] > yvals[ny-1] )
    {
      if ( ybounds[0]      >  88 ) ybounds[0]      =  90;
      if ( ybounds[2*ny-1] < -88 ) ybounds[2*ny-1] = -90;
    }
  else
    {
      if ( ybounds[0]      < -88 ) ybounds[0]      = -90;
      if ( ybounds[2*ny-1] >  88 ) ybounds[2*ny-1] =  90;
    }
}


void gridGenYboundsM(int ny, double *yvals, double *ybounds)
{
  int i;

  for ( i = 0; i < ny-1; i++ )
    {
      ybounds[2*i+1]   = 0.5*(yvals[i] + yvals[i+1]);
      ybounds[2*(i+1)] = 0.5*(yvals[i] + yvals[i+1]);
    }

  ybounds[0]      = 2*yvals[0] - ybounds[1];
  ybounds[2*ny-1] = 2*yvals[ny-1] - ybounds[2*(ny-1)];
}


void gridGenRotBounds(int gridID, int nx, int ny,
		      double *xbounds, double *ybounds, double *xbounds2D, double *ybounds2D)
{
  long i, j, index;
  double minlon, maxlon;
  double minlat, maxlat;
  double xpole, ypole, angle;

  xpole = gridInqXpole(gridID);
  ypole = gridInqYpole(gridID);
  angle = gridInqAngle(gridID);

  for ( j = 0; j < ny; j++ )
    {
      if ( ybounds[0] > ybounds[1] )
	{
	  maxlat = ybounds[2*j];
	  minlat = ybounds[2*j+1];
	}
      else
	{
	  maxlat = ybounds[2*j+1];
	  minlat = ybounds[2*j];
	}

      for ( i = 0; i < nx; i++ )
	{
	  minlon = xbounds[2*i];
	  maxlon = xbounds[2*i+1];

	  index = j*4*nx + 4*i;
	  xbounds2D[index+0] = lamrot_to_lam(minlat, minlon, ypole, xpole, angle);
	  xbounds2D[index+1] = lamrot_to_lam(minlat, maxlon, ypole, xpole, angle);
	  xbounds2D[index+2] = lamrot_to_lam(maxlat, maxlon, ypole, xpole, angle);
	  xbounds2D[index+3] = lamrot_to_lam(maxlat, minlon, ypole, xpole, angle);

	  ybounds2D[index+0] = phirot_to_phi(minlat, minlon, ypole, angle);
	  ybounds2D[index+1] = phirot_to_phi(minlat, maxlon, ypole, angle);
	  ybounds2D[index+2] = phirot_to_phi(maxlat, maxlon, ypole, angle);
	  ybounds2D[index+3] = phirot_to_phi(maxlat, minlon, ypole, angle);
	}
    }
}

static
void gridGenXbounds2D(int nx, int ny, const double * restrict xbounds, double * restrict xbounds2D)
{
  long i, j, index;
  double minlon, maxlon;

#if defined (_OPENMP)
#pragma omp parallel for default(none)        \
  shared(nx, ny, xbounds, xbounds2D)	      \
  private(i, j, minlon, maxlon, index)
#endif
  for ( i = 0; i < nx; ++i )
    {
      minlon = xbounds[2*i  ];
      maxlon = xbounds[2*i+1];

      for ( j = 0; j < ny; ++j )
	{
	  index = j*4*nx + 4*i;
	  xbounds2D[index  ] = minlon;
	  xbounds2D[index+1] = maxlon;
	  xbounds2D[index+2] = maxlon;
	  xbounds2D[index+3] = minlon;
	}
    }
}

static
void gridGenYbounds2D(int nx, int ny, const double * restrict ybounds, double * restrict ybounds2D)
{
  long i, j, index;
  double minlat, maxlat;

#if defined (_OPENMP)
#pragma omp parallel for default(none)        \
  shared(nx, ny, ybounds, ybounds2D)	      \
  private(i, j, minlat, maxlat, index)
#endif
  for ( j = 0; j < ny; ++j )
    {
      if ( ybounds[0] > ybounds[1] )
	{
	  maxlat = ybounds[2*j  ];
	  minlat = ybounds[2*j+1];
	}
      else
	{
	  maxlat = ybounds[2*j+1];
	  minlat = ybounds[2*j  ];
	}

      for ( i = 0; i < nx; ++i )
	{
	  index = j*4*nx + 4*i;
	  ybounds2D[index  ] = minlat;
	  ybounds2D[index+1] = minlat;
	  ybounds2D[index+2] = maxlat;
	  ybounds2D[index+3] = maxlat;
	}
    }
}

static
char *gen_param(const char *fmt, ...)
{
  va_list args;
  char str[256];
  char *rstr;
  int len;

  va_start(args, fmt);

  len = vsprintf(str, fmt, args);

  va_end(args);

  len++;
  rstr = (char *) malloc(len*sizeof(char));
  memcpy(rstr, str, len*sizeof(char));

  return (rstr);
}

static
void lcc_to_geo(int gridID, int gridsize, double *xvals, double *yvals)
{
  double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
  double zlat, zlon;
  double xi, xj;
  int projflag, scanflag;
  int status;
  long i;
  proj_info_t proj;

  gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
	     &projflag, &scanflag);
  /*
    while ( originLon < 0 ) originLon += 360;
    while ( lonParY   < 0 ) lonParY   += 360;
  */
  if ( IS_NOT_EQUAL(xincm, yincm) )
    Warning("X and Y increment must be equal on Lambert Conformal grid (Xinc = %g, Yinc = %g)\n", 
	    xincm, yincm);
  /*
  if ( IS_NOT_EQUAL(lat1, lat2) )
    Warning("Lat1 and Lat2 must be equal on Lambert Conformal grid (Lat1 = %g, Lat2 = %g)\n", 
	    lat1, lat2);
  */
  map_set(PROJ_LC, originLat, originLon, xincm, lonParY, lat1, lat2, &proj);

  for ( i = 0; i < gridsize; i++ )
    {
      xi = xvals[i];
      xj = yvals[i];
      // status = W3FB12(xi, xj, originLat, originLon, xincm, lonParY, lat1, &zlat, &zlon);
      ijll_lc(xi, xj, proj, &zlat, &zlon);
      xvals[i] = zlon;
      yvals[i] = zlat;
    }
}

static
void sinusoidal_to_geo(int gridsize, double *xvals, double *yvals)
{
#if defined (HAVE_LIBPROJ)
  PJ   *libProj;
  char *params[20];
  int nbpar=0;
  projUV data, res;
  long i;

  nbpar = 0;
  params[nbpar++] = (char*) "proj=sinu";
  params[nbpar++] = (char*) "ellps=WGS84";

  if ( cdoVerbose )
    for ( i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%ld] = %s", i+1, params[i]);

  libProj = pj_init(nbpar, params);
  if ( !libProj )
    cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  /* libProj->over = 1; */		/* allow longitude > 180° */

  for ( i = 0; i < gridsize; i++ )
    {
      data.u = xvals[i];
      data.v = yvals[i];
      res = pj_inv(data, libProj);
      xvals[i] = res.u*rad2deg;
      yvals[i] = res.v*rad2deg;
    }
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}

static
void laea_to_geo(int gridID, int gridsize, double *xvals, double *yvals)
{
#if defined (HAVE_LIBPROJ)
  PJ   *libProj;
  char *params[20];
  int nbpar=0;
  projUV data, res;
  double a, lon_0, lat_0;
  long i;

  gridInqLaea(gridID, &a , &lon_0, &lat_0);

  nbpar = 0;
  params[nbpar++] = gen_param("proj=laea");
  params[nbpar++] = gen_param("a=%g", a);
  params[nbpar++] = gen_param("lon_0=%g", lon_0);
  params[nbpar++] = gen_param("lat_0=%g", lat_0);

  if ( cdoVerbose )
    for ( i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

  libProj = pj_init(nbpar, &params[0]);
  if ( !libProj )
    cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( i = 0; i < nbpar; ++i ) free(params[i]);

  /* libProj->over = 1; */		/* allow longitude > 180° */

  for ( i = 0; i < gridsize; i++ )
    {
      data.u = xvals[i];
      data.v = yvals[i];
      res = pj_inv(data, libProj);
      xvals[i] = res.u*rad2deg;
      yvals[i] = res.v*rad2deg;
    }
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}

static
void lcc2_to_geo(int gridID, int gridsize, double *xvals, double *yvals)
{
#if defined (HAVE_LIBPROJ)
  PJ   *libProj;
  char *params[20];
  int nbpar=0;
  projUV data, res;
  double a, lon_0, lat_0, lat_1, lat_2;
  long i;

  gridInqLcc2(gridID, &a , &lon_0, &lat_0, &lat_1, &lat_2);

  nbpar = 0;
  params[nbpar++] = gen_param("proj=lcc");
  if ( a > 0 ) params[nbpar++] = gen_param("a=%g", a);
  params[nbpar++] = gen_param("lon_0=%g", lon_0);
  params[nbpar++] = gen_param("lat_0=%g", lat_0);
  params[nbpar++] = gen_param("lat_1=%g", lat_1);
  params[nbpar++] = gen_param("lat_2=%g", lat_2);

  if ( cdoVerbose )
    for ( i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);
  
  libProj = pj_init(nbpar, &params[0]);
  if ( !libProj )
    cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( i = 0; i < nbpar; ++i ) free(params[i]);

  /* libProj->over = 1; */		/* allow longitude > 180° */
  
  for ( i = 0; i < gridsize; i++ )
    {
      data.u = xvals[i];
      data.v = yvals[i];
      res = pj_inv(data, libProj);
      xvals[i] = res.u*rad2deg;
      yvals[i] = res.v*rad2deg;
    }
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}

int    qu2reg3(double *pfield, int *kpoint, int klat, int klon,
	       double msval, int *kret, int omisng, int operio, int oveggy);

void field2regular(int gridID1, int gridID2, double missval, double *array, int nmiss)
{
  int nlon, nlat;
  int gridtype;
  int lmiss, lperio, lveggy;
  int iret;
  int *rowlonptr;

  gridtype = gridInqType(gridID1);

  if ( gridtype != GRID_GAUSSIAN_REDUCED ) Error("Not a reduced gaussian grid!");

  nlat = gridInqYsize(gridID1);
  nlon = 2*nlat;

  rowlonptr = (int *) malloc(nlat*sizeof(int));

  if ( gridInqSize(gridID2) != nlon*nlat ) Error("Gridsize differ!");

  gridInqRowlon(gridID1, rowlonptr);

  lmiss = nmiss > 0;
  lperio = 1;
  lveggy = 0;

  (void) qu2reg3(array, rowlonptr, nlat, nlon, missval, &iret, lmiss, lperio, lveggy);

  free(rowlonptr);
}


int gridToRegular(int gridID1)
{
  int gridID2;
  int gridtype, gridsize;
  int nx, ny;
  long i;
  double *xvals = NULL, *yvals = NULL;

  gridtype = gridInqType(gridID1);

  if ( gridtype != GRID_GAUSSIAN_REDUCED ) Error("Not a reduced gaussian grid!");

  ny = gridInqYsize(gridID1);
  nx = 2*ny;
  gridsize = nx*ny;

  gridID2  = gridCreate(GRID_GAUSSIAN, gridsize);
	  
  gridDefXsize(gridID2, nx);
  gridDefYsize(gridID2, ny);
  
  xvals = (double *) malloc(nx*sizeof(double));
  yvals = (double *) malloc(ny*sizeof(double));

  for ( i = 0; i < nx; ++i ) xvals[i] = i * 360./nx;
  gridInqYvals(gridID1, yvals);

  gridDefXvals(gridID2, xvals);
  gridDefYvals(gridID2, yvals);

  free(xvals);
  free(yvals);

  return (gridID2);
}

static
void gridCopyMask(int gridID1, int gridID2, long gridsize)
{
  if ( gridInqMask(gridID1, NULL) )
    {
      int *mask;
      mask = (int *) malloc(gridsize*sizeof(int));
      gridInqMask(gridID1, mask);
      gridDefMask(gridID2, mask);
      free(mask);
    }
}

static
int check_range(long n, double *vals, double valid_min, double valid_max)
{
  int status = 0;
  long i;

  for ( i = 0; i < n; ++i )
    {
      if ( vals[i] < valid_min || vals[i] > valid_max )
	{
	  status = 1;
	  break;
	}
    }

  return (status);
}

int gridToCurvilinear(int gridID1)
{
  int gridID2;
  int gridtype, gridsize;
  long index;

  gridtype = gridInqType(gridID1);
  gridsize = gridInqSize(gridID1);
  gridID2  = gridCreate(GRID_CURVILINEAR, gridsize);
  gridDefPrec(gridID2, DATATYPE_FLT32);
	  
  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_LCC:
    case GRID_LCC2:
    case GRID_LAEA:
    case GRID_SINUSOIDAL:
      {
	long i, j;
	int nx, ny;
	double *xvals = NULL, *yvals = NULL;
	double *xvals2D, *yvals2D;
	double *xbounds = NULL, *ybounds = NULL;
	double *xbounds2D, *ybounds2D;
	char xunits[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
	double xscale = 1, yscale = 1;

	nx = gridInqXsize(gridID1);
	ny = gridInqYsize(gridID1);

	gridInqXunits(gridID1, xunits);
	gridInqYunits(gridID1, yunits);

	if ( memcmp(xunits, "km", 2) == 0 ) xscale = 1000;
	if ( memcmp(yunits, "km", 2) == 0 ) yscale = 1000;

	gridDefXsize(gridID2, nx);
	gridDefYsize(gridID2, ny);

	xvals2D = (double *) malloc(gridsize*sizeof(double));
	yvals2D = (double *) malloc(gridsize*sizeof(double));


	if ( gridtype == GRID_LCC )
	  {
	    for ( j = 0; j < ny; j++ )
	      for ( i = 0; i < nx; i++ )
		{
		  xvals2D[j*nx+i] = i+1;
		  yvals2D[j*nx+i] = j+1;
		}

	    lcc_to_geo(gridID1, gridsize, xvals2D, yvals2D);
	  }
	else
	  {
	    if ( ! (gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL)) )
	      Error("Grid has no values");

	    xvals = (double *) malloc(nx*sizeof(double));
	    yvals = (double *) malloc(ny*sizeof(double));

	    gridInqXvals(gridID1, xvals);
	    gridInqYvals(gridID1, yvals);

	    if ( gridIsRotated(gridID1) )
	      {
		double xpole, ypole, angle;
		
		xpole = gridInqXpole(gridID1);
		ypole = gridInqYpole(gridID1);
		angle = gridInqAngle(gridID1);
		
		for ( j = 0; j < ny; j++ )
		  for ( i = 0; i < nx; i++ )
		    {
		      xvals2D[j*nx+i] = lamrot_to_lam(yvals[j], xvals[i], ypole, xpole, angle);
		      yvals2D[j*nx+i] = phirot_to_phi(yvals[j], xvals[i], ypole, angle);
		    }	    
	      }
	    else
	      {
		for ( j = 0; j < ny; j++ )
		  for ( i = 0; i < nx; i++ )
		    {
		      xvals2D[j*nx+i] = xscale*xvals[i];
		      yvals2D[j*nx+i] = yscale*yvals[j];
		    }

		if ( gridtype == GRID_SINUSOIDAL )
		  {
		    sinusoidal_to_geo(gridsize, xvals2D, yvals2D);
		    /* correct_sinxvals(nx, ny, xvals2D); */
		  }
		else if ( gridtype == GRID_LAEA )
		  laea_to_geo(gridID1, gridsize, xvals2D, yvals2D);
		else if ( gridtype == GRID_LCC2 )
		  lcc2_to_geo(gridID1, gridsize, xvals2D, yvals2D);
	      }
	  }

	gridDefXvals(gridID2, xvals2D);
	gridDefYvals(gridID2, yvals2D);

	if ( xvals2D ) free(xvals2D);
	if ( yvals2D ) free(yvals2D);

	if ( gridtype == GRID_LCC )
	  {		
	    xbounds2D = (double *) malloc(4*gridsize*sizeof(double));
	    ybounds2D = (double *) malloc(4*gridsize*sizeof(double));

	    for ( j = 0; j < ny; j++ )
	      for ( i = 0; i < nx; i++ )
		{
		  index = j*4*nx + 4*i;

		  xbounds2D[index+0] = i+1.5;
		  ybounds2D[index+0] = j+1.5;

		  xbounds2D[index+1] = i+0.5;
		  ybounds2D[index+1] = j+1.5;

		  xbounds2D[index+2] = i+0.5;
		  ybounds2D[index+2] = j+0.5;

		  xbounds2D[index+3] = i+1.5;
		  ybounds2D[index+3] = j+0.5;
		}

	    lcc_to_geo(gridID1, 4*gridsize, xbounds2D, ybounds2D);

	    gridDefXbounds(gridID2, xbounds2D);
	    gridDefYbounds(gridID2, ybounds2D);

	    free(xbounds2D);
	    free(ybounds2D);
	  }
	else
	  {
	    if ( gridInqXbounds(gridID1, NULL) )
	      {
		xbounds = (double *) malloc(2*nx*sizeof(double));
		gridInqXbounds(gridID1, xbounds);
		if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
		  if ( check_range(2*nx, xbounds, -720, 720) )
		    {
		      cdoWarning("longitude bounds out of range, skipped!");
		      free(xbounds);
		      xbounds = NULL;
		    }
	      }
	    else if ( nx > 1 )
	      {
		xbounds = (double *) malloc(2*nx*sizeof(double));
		gridGenXbounds(nx, xvals, xbounds);
	      }

	    if ( gridInqYbounds(gridID1, NULL) )
	      {
		ybounds = (double *) malloc(2*ny*sizeof(double));
		gridInqYbounds(gridID1, ybounds);
		if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
		  if ( check_range(2*ny, ybounds, -180, 180) )
		    {
		      cdoWarning("latitude bounds out of range, skipped!");
		      free(ybounds);
		      ybounds = NULL;
		    }
	      }
	    else if ( ny > 1 )
	      {
		ybounds = (double *) malloc(2*ny*sizeof(double));
		if ( gridtype == GRID_SINUSOIDAL || 
		     gridtype == GRID_LAEA       || 
		     gridtype == GRID_LCC2 )
		  gridGenYboundsM(ny, yvals, ybounds);
		else
		  gridGenYbounds(ny, yvals, ybounds);
	      }

	    if ( xvals ) free(xvals);
	    if ( yvals ) free(yvals);

	    if ( xbounds && ybounds )
	      {
		xbounds2D = (double *) malloc(4*gridsize*sizeof(double));
		ybounds2D = (double *) malloc(4*gridsize*sizeof(double));

		if ( gridIsRotated(gridID1) )
		  {
		    gridGenRotBounds(gridID1, nx, ny, xbounds, ybounds, xbounds2D, ybounds2D);
		  }
		else
		  {
		    if ( gridtype == GRID_SINUSOIDAL ||
			 gridtype == GRID_LAEA       || 
			 gridtype == GRID_LCC2 )
		      {
			for ( j = 0; j < ny; j++ )
			  for ( i = 0; i < nx; i++ )
			    {
			      index = j*4*nx + 4*i;

			      xbounds2D[index+0] = xscale*xbounds[2*i];
			      ybounds2D[index+0] = yscale*ybounds[2*j];

			      xbounds2D[index+1] = xscale*xbounds[2*i];
			      ybounds2D[index+1] = yscale*ybounds[2*j+1];

			      xbounds2D[index+2] = xscale*xbounds[2*i+1];
			      ybounds2D[index+2] = yscale*ybounds[2*j+1];

			      xbounds2D[index+3] = xscale*xbounds[2*i+1];
			      ybounds2D[index+3] = yscale*ybounds[2*j];
			    }
			
			if ( gridtype == GRID_SINUSOIDAL )
			  {
			    sinusoidal_to_geo(4*gridsize, xbounds2D, ybounds2D);
			    /*
			    xvals2D = (double *) malloc(gridsize*sizeof(double));
			    for ( j = 0; j < 4; ++j )
			      {
				for ( i = 0; i < gridsize; ++i ) xvals2D[i] = xbounds2D[i*4+j];
				correct_sinxvals(nx, ny, xvals2D);
				for ( i = 0; i < gridsize; ++i ) xbounds2D[i*4+j] = xvals2D[i];
			      }
			    free(xvals2D);
			    */
			  }
			else if ( gridtype == GRID_LAEA )
			  laea_to_geo(gridID1, 4*gridsize, xbounds2D, ybounds2D);
			else if ( gridtype == GRID_LCC2 )
			  lcc2_to_geo(gridID1, 4*gridsize, xbounds2D, ybounds2D);
		      }
		    else
		      {
			gridGenXbounds2D(nx, ny, xbounds, xbounds2D);
			gridGenYbounds2D(nx, ny, ybounds, ybounds2D);
		      }
		  }
		
		gridDefXbounds(gridID2, xbounds2D);
		gridDefYbounds(gridID2, ybounds2D);
		
		if ( xbounds )  free(xbounds);
		if ( ybounds )  free(ybounds);
		if ( xbounds2D) free(xbounds2D);
		if ( ybounds2D) free(ybounds2D);
	      }
	  }

	gridCopyMask(gridID1, gridID2, gridsize);

	break;
      }
    default:
      {
	Error("Grid type >%s< unsupported!", gridNamePtr(gridtype));
	break;
      }
    }

  return (gridID2);
}


int gridToUnstructured(int gridID1)
{
  int gridID2;
  int gridtype, gridsize;

  gridtype = gridInqType(gridID1);
  gridsize = gridInqSize(gridID1);
  gridID2  = gridCreate(GRID_UNSTRUCTURED, gridsize);
  gridDefPrec(gridID2, DATATYPE_FLT32);
	  
  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
      {
	long i, j;
	int nx, ny;
	double *xvals, *yvals;
	double *xvals2D, *yvals2D;
	double *xbounds = NULL, *ybounds = NULL;
	double *xbounds2D, *ybounds2D;

	gridDefXname(gridID2, "lon");
	gridDefYname(gridID2, "lat");
	gridDefXlongname(gridID2, "longitude");
	gridDefYlongname(gridID2, "latitude");
	gridDefXunits(gridID2, "degrees_east");
	gridDefYunits(gridID2, "degrees_north");

	gridDefNvertex(gridID2, 4);

	nx = gridInqXsize(gridID1);
	ny = gridInqYsize(gridID1);
	 
	gridDefXsize(gridID2, gridsize);
	gridDefYsize(gridID2, gridsize);

	xvals = (double *) malloc(nx*sizeof(double));
	yvals = (double *) malloc(ny*sizeof(double));

	xvals2D = (double *) malloc(gridsize*sizeof(double));
	yvals2D = (double *) malloc(gridsize*sizeof(double));

	gridInqXvals(gridID1, xvals);
	gridInqYvals(gridID1, yvals);

	if ( gridIsRotated(gridID1) )
	  {
	    double xpole, ypole, angle;
	    
	    xpole = gridInqXpole(gridID1);
	    ypole = gridInqYpole(gridID1);
	    angle = gridInqAngle(gridID1);
		
	    for ( j = 0; j < ny; j++ )
	      for ( i = 0; i < nx; i++ )
		{
		  xvals2D[j*nx+i] = lamrot_to_lam(yvals[j], xvals[i], ypole, xpole, angle);
		  yvals2D[j*nx+i] = phirot_to_phi(yvals[j], xvals[i], ypole, angle);
		}
	  }
	else
	  {
	    for ( j = 0; j < ny; j++ )
	      for ( i = 0; i < nx; i++ )
		{
		  xvals2D[j*nx+i] = xvals[i];
		  yvals2D[j*nx+i] = yvals[j];
		}
	  }

	gridDefXvals(gridID2, xvals2D);
	gridDefYvals(gridID2, yvals2D);

	free(xvals2D);
	free(yvals2D);

	if ( gridInqXbounds(gridID1, NULL) )
	  {
	    xbounds = (double *) malloc(2*nx*sizeof(double));
	    gridInqXbounds(gridID1, xbounds);
	  }
	else if ( nx > 1 )
	  {
	    xbounds = (double *) malloc(2*nx*sizeof(double));
	    gridGenXbounds(nx, xvals, xbounds);
	  }

	if ( gridInqYbounds(gridID1, NULL) )
	  {
	    ybounds = (double *) malloc(2*ny*sizeof(double));
	    gridInqYbounds(gridID1, ybounds);
	  }
	else if ( ny > 1 )
	  {
	    ybounds = (double *) malloc(2*ny*sizeof(double));
	    gridGenYbounds(ny, yvals, ybounds);
	  }

	free(xvals);
	free(yvals);

	if ( xbounds && ybounds )
	  {
	    xbounds2D = (double *) malloc(4*gridsize*sizeof(double));
	    ybounds2D = (double *) malloc(4*gridsize*sizeof(double));

	    if ( gridIsRotated(gridID1) )
	      {
		gridGenRotBounds(gridID1, nx, ny, xbounds, ybounds, xbounds2D, ybounds2D);
	      }
	    else
	      {
		gridGenXbounds2D(nx, ny, xbounds, xbounds2D);
		gridGenYbounds2D(nx, ny, ybounds, ybounds2D);
	      }

	    gridDefXbounds(gridID2, xbounds2D);
	    gridDefYbounds(gridID2, ybounds2D);

	    free(xbounds);
	    free(ybounds);
	    free(xbounds2D);
	    free(ybounds2D);
	  }

	gridCopyMask(gridID1, gridID2, gridsize);

	break;
      }
    case GRID_CURVILINEAR:
      {
	gridID2 = gridDuplicate(gridID1);
	gridChangeType(gridID2, GRID_UNSTRUCTURED);
	gridDefXsize(gridID2, gridsize);
	gridDefYsize(gridID2, gridsize);

	break;
      }
    case GRID_GME:
      {
	int nd, ni, ni2, ni3;
	long i, j;
	int nv = 6;
	int *imask;
	double *xvals, *yvals;
	double *xbounds, *ybounds;

	nd  = gridInqGMEnd(gridID1);
	ni  = gridInqGMEni(gridID1);
	ni2 = gridInqGMEni2(gridID1);
	ni3 = gridInqGMEni3(gridID1);

	imask   = (int *) malloc(gridsize*sizeof(int));
	xvals   = (double *) malloc(gridsize*sizeof(double));
	yvals   = (double *) malloc(gridsize*sizeof(double));
	xbounds = (double *) malloc(nv*gridsize*sizeof(double));
	ybounds = (double *) malloc(nv*gridsize*sizeof(double));

	gme_grid(gridsize, xvals, yvals, xbounds, ybounds, imask, ni, nd, ni2, ni3);
	
	for ( i = 0; i < gridsize; i++ )
	  {
	    xvals[i] *= RAD2DEG;
	    yvals[i] *= RAD2DEG;

	    for ( j = 0; j < nv; j++ )
	      {
		xbounds[i*nv + j] *= RAD2DEG;
		ybounds[i*nv + j] *= RAD2DEG;
	      }
	    /* printf("%d %g %g\n", i, xvals[i], yvals[i]); */
	  }
	
	gridDefXsize(gridID2, gridsize);
	gridDefYsize(gridID2, gridsize);

	gridDefXvals(gridID2, xvals);
	gridDefYvals(gridID2, yvals);

	gridDefMaskGME(gridID2, imask);

	gridDefNvertex(gridID2, nv);

	gridDefXbounds(gridID2, xbounds);
	gridDefYbounds(gridID2, ybounds);

	gridDefXunits(gridID2, "degrees_east");
	gridDefYunits(gridID2, "degrees_north");

	free (imask);
	free (xvals);
	free (yvals);
	free (xbounds);
	free (ybounds);
	
	gridCopyMask(gridID1, gridID2, gridsize);

	break;
      }
    default:
      {
	Error("Grid type >%s< unsupported!", gridNamePtr(gridtype));
	break;
      }
    }

  return (gridID2);
}


int referenceToGrid(int gridID1)
{
  int gridID2 = -1;
  int gridtype, gridsize;
  int offset = 7;
  char gridfile[8912];

  gridsize = gridInqSize(gridID1);

  if ( gridInqReference(gridID1, NULL) )
    {
      int streamID;
      int number, position;

      number = gridInqNumber(gridID1);
      position = gridInqPosition(gridID1);
      gridInqReference(gridID1, gridfile);

      if ( gridfile[offset] != '/' ) offset--;
      streamID = streamOpenRead(&gridfile[offset]);
      if ( streamID >= 0 )
	{
	  int vlistID, gridID = -1;
	  int ngrids;
	  vlistID = streamInqVlist(streamID);
	  ngrids = vlistNgrids(vlistID);
	  if ( position > 0 && position <= ngrids )
	    {
	      gridID = vlistGrid(vlistID, position-1);
	      if ( gridInqSize(gridID) == gridsize )
		gridID2 = gridDuplicate(gridID);
	      else
		cdoWarning("Grid size %d on position %d do not match! Reference=%s", gridsize, position, gridfile);
	    }
	  else
	    cdoWarning("Grid position %d not available! Reference=%s", position, gridfile);

	  streamClose(streamID);
	}
      else
	cdoWarning("Reference to grid not found! Path=%s", gridfile);
    }
  else
    {
      cdoWarning("No reference to grid found!");
    }

  return (gridID2);
}


static
double areas(struct cart *dv1, struct cart *dv2, struct cart *dv3)
{
  double a1, a2, a3;
  double ca1, ca2, ca3;
  double s12, s23, s31;

  struct cart u12, u23, u31;

  double areas;

  /* compute cross products Uij = Vi X Vj */

  u12.x[0] = dv1->x[1]*dv2->x[2] - dv1->x[2]*dv2->x[1];
  u12.x[1] = dv1->x[2]*dv2->x[0] - dv1->x[0]*dv2->x[2];
  u12.x[2] = dv1->x[0]*dv2->x[1] - dv1->x[1]*dv2->x[0];
  
  u23.x[0] = dv2->x[1]*dv3->x[2] - dv2->x[2]*dv3->x[1];
  u23.x[1] = dv2->x[2]*dv3->x[0] - dv2->x[0]*dv3->x[2];
  u23.x[2] = dv2->x[0]*dv3->x[1] - dv2->x[1]*dv3->x[0];
  
  u31.x[0] = dv3->x[1]*dv1->x[2] - dv3->x[2]*dv1->x[1];
  u31.x[1] = dv3->x[2]*dv1->x[0] - dv3->x[0]*dv1->x[2];
  u31.x[2] = dv3->x[0]*dv1->x[1] - dv3->x[1]*dv1->x[0];
  
  /* normalize Uij to unit vectors */
  
  s12 = u12.x[0]*u12.x[0]+u12.x[1]*u12.x[1]+u12.x[2]*u12.x[2];
  s23 = u23.x[0]*u23.x[0]+u23.x[1]*u23.x[1]+u23.x[2]*u23.x[2];
  s31 = u31.x[0]*u31.x[0]+u31.x[1]*u31.x[1]+u31.x[2]*u31.x[2];

  /* test for a degenerate triangle associated with collinear vertices */
  
  if ( !(fabs(s12) > 0.0) || !(fabs(s23) > 0.0) || !(fabs(s31) > 0.0) ) {
    areas = 0.0;
    return areas;
  }

  s12 = sqrt(s12);
  s23 = sqrt(s23);
  s31 = sqrt(s31);
  
  u12.x[0] = u12.x[0]/s12; u12.x[1] = u12.x[1]/s12; u12.x[2] = u12.x[2]/s12;
  u23.x[0] = u23.x[0]/s23; u23.x[1] = u23.x[1]/s23; u23.x[2] = u23.x[2]/s23;
  u31.x[0] = u31.x[0]/s31; u31.x[1] = u31.x[1]/s31; u31.x[2] = u31.x[2]/s31;
  
  /*
   *  Compute interior angles Ai as the dihedral angles between planes:
   *  CA1 = cos(A1) = -<U12,U31>
   *  CA2 = cos(A2) = -<U23,U12>
   *  CA3 = cos(A3) = -<U31,U23>
   */

  ca1 = -( u12.x[0]*u31.x[0]+u12.x[1]*u31.x[1]+u12.x[2]*u31.x[2] );
  ca2 = -( u23.x[0]*u12.x[0]+u23.x[1]*u12.x[1]+u23.x[2]*u12.x[2] );
  ca3 = -( u31.x[0]*u23.x[0]+u31.x[1]*u23.x[1]+u31.x[2]*u23.x[2] );

#if ! defined (FMAX)
#define  FMAX(a,b)  ((a) > (b) ? (a) : (b))
#endif
#if ! defined (FMIN)
#define  FMIN(a,b)  ((a) < (b) ? (a) : (b))
#endif

  ca1 = FMAX(ca1, -1.0);
  ca1 = FMIN(ca1, +1.0);
  ca2 = FMAX(ca2, -1.0);
  ca2 = FMIN(ca2, +1.0);
  ca3 = FMAX(ca3, -1.0);
  ca3 = FMIN(ca3, +1.0);
  
  a1 = acos(ca1);
  a2 = acos(ca2);
  a3 = acos(ca3);
  
  /* compute AREAS = A1 + A2 + A3 - PI */
  
  areas = a1 + a2 + a3 - M_PI;

  if ( areas < 0.0 ) {
    areas = 0.0;
  }

  return areas;
}

static
double cell_area(long i, long nv, double *grid_center_lon, double *grid_center_lat,
		 double *grid_corner_lon, double *grid_corner_lat, int *status)
{
  long k;
  double xa;
  double area;
  double lim = 179*deg2rad;
  struct geo p1, p2, p3;
  struct cart c1, c2, c3;

  area = 0;
      
  p3.lon = grid_center_lon[i]*deg2rad; 
  p3.lat = grid_center_lat[i]*deg2rad;
  c3 = gc2cc(&p3);
      
  for ( k = 1; k < nv; ++k )
    {
      p1.lon = grid_corner_lon[i*nv+k-1]*deg2rad; 
      p1.lat = grid_corner_lat[i*nv+k-1]*deg2rad;
      c1 = gc2cc(&p1);
      p2.lon = grid_corner_lon[i*nv+k]*deg2rad; 
      p2.lat = grid_corner_lat[i*nv+k]*deg2rad;
      c2 = gc2cc(&p2);

      xa = areas(&c1, &c2, &c3);
      /*
      if ( xa > 0.001 )
	if ( (fabs(p1.lon - p2.lon) > lim) ||
	     (fabs(p2.lon - p3.lon) > lim) ||
	     (fabs(p3.lon - p1.lon) > lim) )
	  {
	    printf("%ld %ld xa %g\n", i, k, xa);
	    printf("  dlon %g %g %g %g\n", lim, fabs(p1.lon - p2.lon), fabs(p2.lon - p3.lon), fabs(p3.lon - p1.lon));
	    *status = 2;
	  }
      */
      area += xa;
    }

  p1.lon = grid_corner_lon[i*nv+0]*deg2rad; 
  p1.lat = grid_corner_lat[i*nv+0]*deg2rad;
  c1 = gc2cc(&p1);
  p2.lon = grid_corner_lon[i*nv+nv-1]*deg2rad; 
  p2.lat = grid_corner_lat[i*nv+nv-1]*deg2rad;
  c2 = gc2cc(&p2);

  xa = areas(&c1, &c2, &c3);
  /*
  if ( xa > 0.001 )
    if ( (fabs(p1.lon - p2.lon) > lim) ||
	 (fabs(p2.lon - p3.lon) > lim) ||
	 (fabs(p3.lon - p1.lon) > lim) )
      {
	printf("%ld %ld xa %g\n", i, k, xa);
	printf("  dlon %g %g %g %g\n", lim, fabs(p1.lon - p2.lon), fabs(p2.lon - p3.lon), fabs(p3.lon - p1.lon));
	*status = 2;
      }
  */
  area += xa;

  return (area);
}


int gridGenArea(int gridID, double *area)
{
  int status = 0;
  int gridtype;
  int lgrid_gen_bounds = FALSE;
  int lgriddestroy = FALSE;
  long i;
  long nv, gridsize;
  double total_area;
  double *grid_center_lon = NULL;
  double *grid_center_lat = NULL;
  double *grid_corner_lon = NULL;
  double *grid_corner_lat = NULL;
  int *grid_mask = NULL;

  gridsize = gridInqSize(gridID);
  gridtype = gridInqType(gridID);

  if ( gridtype != GRID_LONLAT      &&
       gridtype != GRID_GAUSSIAN    &&
       gridtype != GRID_LCC         &&
       gridtype != GRID_LCC2        &&
       gridtype != GRID_LAEA        &&
       gridtype != GRID_SINUSOIDAL  &&
       gridtype != GRID_GME         &&
       gridtype != GRID_REFERENCE   &&
       gridtype != GRID_CURVILINEAR &&
       gridtype != GRID_UNSTRUCTURED )
    {
      cdoAbort("Internal error! Unsupported gridtype: %s", gridNamePtr(gridtype)); 
    }

  if ( gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR )
    {
      if ( gridtype == GRID_GME )
	{
	  lgriddestroy = TRUE;
	  gridID = gridToUnstructured(gridID);
	  grid_mask = (int *) malloc(gridsize*sizeof(int));
	  gridInqMaskGME(gridID, grid_mask);
	}
      else if ( gridtype == GRID_REFERENCE )
	{
	  lgriddestroy = TRUE;
	  gridID = referenceToGrid(gridID);
	  if ( gridID == -1 ) return (1);
	}
      else
	{
	  lgriddestroy = TRUE;
	  gridID = gridToCurvilinear(gridID);
	  lgrid_gen_bounds = TRUE;
	}
    }

  gridtype = gridInqType(gridID);

  if ( gridtype == GRID_UNSTRUCTURED )
    nv = gridInqNvertex(gridID);
  else
    nv = 4;
  
  grid_center_lon = (double *) malloc(gridsize*sizeof(double));
  grid_center_lat = (double *) malloc(gridsize*sizeof(double));

  gridInqXvals(gridID, grid_center_lon);
  gridInqYvals(gridID, grid_center_lat);

  grid_corner_lon = (double *) malloc(nv*gridsize*sizeof(double));
  grid_corner_lat = (double *) malloc(nv*gridsize*sizeof(double));

  if ( gridInqYbounds(gridID, NULL) && gridInqXbounds(gridID, NULL) )
    {
      gridInqXbounds(gridID, grid_corner_lon);
      gridInqYbounds(gridID, grid_corner_lat);
    }
  else
    {
      if ( lgrid_gen_bounds )
	{
	  int nlon = gridInqXsize(gridID);
	  int nlat = gridInqYsize(gridID);
	  double dlon = 0;
	  if ( nlon == 1 )
	    {
	      dlon = 1;
	    }
	  genXbounds(nlon, nlat, grid_center_lon, grid_corner_lon, dlon);
	  genYbounds(nlon, nlat, grid_center_lat, grid_corner_lat);
	}
      else
	{
	  status = 1;
	  return (status);
	}
    }
  
  /* Convert lat/lon units if required */
  {
    char units[CDI_MAX_NAME];

    gridInqXunits(gridID, units);

    if ( memcmp(units, "degree", 6) == 0 )
      {
	/* No conversion necessary */
      }
    else if ( memcmp(units, "radian", 6) == 0 )
      {
	for ( i = 0; i < gridsize; ++i )
	  {
	    grid_center_lon[i] *= rad2deg;
	    grid_center_lat[i] *= rad2deg;
	  }
	for ( i = 0; i < gridsize*nv; ++i )
	  {
	    grid_corner_lon[i] *= rad2deg;
	    grid_corner_lat[i] *= rad2deg;
	  }
      }
    else
      {
	cdoWarning("Unknown units supplied for grid1 center lat/lon: proceeding assuming radians");
      }
  }

  if ( lgriddestroy ) gridDestroy(gridID);

  total_area = 0;
#if defined (_OPENMP)
#pragma omp parallel for default(none)        \
  shared(gridsize, area, nv, grid_center_lon, grid_center_lat, grid_corner_lon, grid_corner_lat, status) \
  private(i)
#endif
  for ( i = 0; i < gridsize; ++i )
    {
      area[i] = cell_area(i, nv, grid_center_lon, grid_center_lat, grid_corner_lon, grid_corner_lat, &status);
      //     total_area += area[i];
    }

  //  if ( cdoVerbose ) cdoPrint("Total area = %g", total_area);

  free(grid_center_lon);
  free(grid_center_lat);
  free(grid_corner_lon);
  free(grid_corner_lat);
  if ( grid_mask ) free(grid_mask);

  return (status);
}


int gridGenWeights(int gridID, double *grid_area, double *grid_wgts)
{
  int i, nvals, gridsize, gridtype;
  int status = 0;
  int *grid_mask = NULL;
  double total_area;

  gridtype = gridInqType(gridID);
  gridsize = gridInqSize(gridID);
  
  if ( gridtype == GRID_GME )
    {
      gridID = gridToUnstructured(gridID);	  
      grid_mask = (int *) malloc(gridsize*sizeof(int));
      gridInqMaskGME(gridID, grid_mask);
    }

  total_area = 0;
  nvals = 0;
  for ( i = 0; i < gridsize; i++ )
    {
      if ( grid_mask )
	if ( grid_mask[i] == 0 ) continue;
      total_area += grid_area[i];
      nvals++;
    }

  if ( cdoVerbose ) cdoPrint("Total area = %g", total_area);

  for ( i = 0; i < gridsize; i++ )
    {
      if ( grid_mask )
	if ( grid_mask[i] == 0 )
	  {
	    grid_wgts[i] = 0;
	    continue;
	  }
      
      grid_wgts[i] = grid_area[i] / total_area;
    }
  
  if ( grid_mask ) free(grid_mask);

  return (status);
}


static
int gridWeightsOld(int gridID, double *weights)
{
  int status = FALSE;
  long i, j;
  int len;

  len = gridInqSize(gridID);

  if ( gridHasArea(gridID) )
    {
      gridInqArea(gridID, weights);
    }
  else
    {
      int gridtype = gridInqType(gridID);

      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
	{
	  int     nlat, nlon;
	  int     datapoint;
	  double *lats = NULL, *lons = NULL;
	  double sumw;
	  double phi1, phi2, theta1, theta2, sindphi;

	  nlon = gridInqXsize(gridID);
	  nlat = gridInqYsize(gridID);

	  lons = 1 + (double *) malloc((nlon+2)*sizeof(double));
	  lats = 1 + (double *) malloc((nlat+2)*sizeof(double));

	  gridInqXvals(gridID, lons);
	  gridInqYvals(gridID, lats);

	  /* Interpolate to find latitudes outside boundaries. */
	  lats[-1]   = 2*lats[0] - lats[1];
	  lats[nlat] = 2*lats[nlat-1] - lats[nlat-2];
	  lons[-1]   = 2*lons[0] - lons[1];
	  lons[nlon] = 2*lons[nlon-1] - lons[nlon-2];
  
	  /*  Calculate weights.  */
	  /*  phi 1 and 2 and theta 1 and 2 represent respectively the boundary */
	  /*  latitudes and longitudes of a particular grid square.             */
	  datapoint = 0;
	  sumw = 0;
	  for ( j = 0; j < nlat; j++ )
	    {
	      phi1 = (lats[j-1]+lats[j])/2*deg2rad;
	      phi2 = (lats[j+1]+lats[j])/2*deg2rad;
	      if ( phi1 < (-1*M_PI/2) ) phi1 = -1*M_PI/2;
	      if ( phi1 > (   M_PI/2) ) phi1 =    M_PI/2;
	      if ( phi2 > (   M_PI/2) ) phi2 =    M_PI/2;
	      if ( phi2 < (-1*M_PI/2) ) phi2 = -1*M_PI/2;
	      sindphi = sin(phi2)-sin(phi1);
	      for( i = 0; i < nlon; i++ )
		{
		  if ( lons[i] >= lons[0]+360 || fabs(lats[j]) > 90 )
		    weights[datapoint] = 0;
		  else
		    {
		      theta1 = (lons[i-1]+lons[i])/2*deg2rad;
		      theta2 = (lons[i+1]+lons[i])/2*deg2rad;
		      weights[datapoint] = fabs((theta2-theta1)*sindphi);
		      sumw += weights[datapoint];
		    }
		  datapoint++;
		}
	    }

	  /* Normalise weights.  */
	  if( IS_NOT_EQUAL(sumw, 0) )
	    for( i = 0; i < datapoint; i++ ) weights[i] /= sumw;

	  if ( lons-1 ) free(lons-1);
	  if ( lats-1 ) free(lats-1);
	}
      else
	{
	  status = TRUE;

	  for ( i = 0; i < len; i++ ) weights[i] = 1./len;
	}
    }

  return (status);
}


int gridWeights(int gridID, double *grid_wgts)
{
  int i, gridsize, gridtype;
  int a_status, w_status;
  double *grid_area;

  gridtype = gridInqType(gridID);
  gridsize = gridInqSize(gridID);
  
  grid_area = (double *) malloc(gridsize*sizeof(double));

  a_status = 0;

  if ( gridHasArea(gridID) )
    {
      if ( cdoVerbose ) cdoPrint("Using existing grid cell area!");
      gridInqArea(gridID, grid_area);
    }
  else
    {
      if ( gridtype == GRID_LONLAT      ||
	   gridtype == GRID_GAUSSIAN    ||
	   gridtype == GRID_LCC         ||
	   gridtype == GRID_LCC2        ||
	   gridtype == GRID_LAEA        ||
	   gridtype == GRID_SINUSOIDAL  ||
	   gridtype == GRID_GME         ||
	   gridtype == GRID_REFERENCE   ||
	   gridtype == GRID_CURVILINEAR ||
	   gridtype == GRID_UNSTRUCTURED )
	{
	  a_status = gridGenArea(gridID, grid_area);
	}
      else
	{
	  a_status = 1;
	}
    }

  if ( a_status == 0 )
    {
      w_status = gridGenWeights(gridID, grid_area, grid_wgts);
    }
  else
    {
      for ( i = 0; i < gridsize; ++i )
	grid_wgts[i] = 1./gridsize;

      w_status = 1;
    }
  /*
  for ( i = 0; i < gridsize; ++i ) 
    printf("weights: %d %d %d %g %g\n", a_status, w_status, i, grid_area[i], grid_wgts[i]);
  */
  free(grid_area);

  return (w_status);
}
