/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Selbox     sellonlatbox    Select lon/lat box
      Selbox     selindexbox     Select index box
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"


static
void correct_xvals(long nlon, long inc, double *xvals)
{
  long i; 

  if ( IS_EQUAL(xvals[0], xvals[(nlon-1)*inc]) ) xvals[(nlon-1)*inc] += 360;

  if ( xvals[0] > xvals[(nlon-1)*inc] )
    for ( i = 0; i < nlon; i++ )
      if ( xvals[i*inc] >= 180 ) xvals[i*inc] -= 360;

  for ( i = 0; i < nlon; i++ )
    {
      if ( xvals[i*inc] < -180 ) xvals[i*inc] += 360;
      if ( xvals[i*inc] >  360 ) xvals[i*inc] -= 360;
    }

  if ( xvals[0] > xvals[(nlon-1)*inc] )
    for ( i = 1; i < nlon; i++ )
      if ( xvals[i*inc] < xvals[(i-1)*inc] ) xvals[i*inc] += 360;
}

static
int gengrid(int gridID1, int lat1, int lat2, int lon11, int lon12, int lon21, int lon22)
{
  int gridtype, gridID2;
  int nlon1, nlat1;
  int nlon2, nlat2;
  int nlon21, nlon22;
  int i;
  int prec;
  int ilat, ilon;
  int lxvals, lyvals;
  char xname[CDI_MAX_NAME], xlongname[CDI_MAX_NAME], xunits[CDI_MAX_NAME];
  char yname[CDI_MAX_NAME], ylongname[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
  double *xvals1 = NULL, *yvals1 = NULL;
  double *xvals2 = NULL, *yvals2 = NULL;
  double *xbounds1 = NULL, *ybounds1 = NULL;
  double *xbounds2 = NULL, *ybounds2 = NULL;
  double *pxvals2 = NULL, *pyvals2 = NULL;
  double *pxbounds2 = NULL, *pybounds2 = NULL;

  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  nlon21 = lon12 - lon11 + 1;
  nlon22 = lon22 - lon21 + 1;
  nlon2 = nlon21 + nlon22;
  nlat2 = lat2 - lat1 + 1;

  gridtype = gridInqType(gridID1);
  prec     = gridInqPrec(gridID1);

  gridID2 = gridCreate(gridtype, nlon2*nlat2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  gridDefNP(gridID2, gridInqNP(gridID1));

  gridDefPrec(gridID2, prec);

  gridInqXname(gridID1, xname);
  gridInqXlongname(gridID1, xlongname);
  gridInqXunits(gridID1, xunits);
  gridInqYname(gridID1, yname);
  gridInqYlongname(gridID1, ylongname);
  gridInqYunits(gridID1, yunits);

  gridDefXname(gridID2, xname);
  gridDefXlongname(gridID2, xlongname);
  gridDefXunits(gridID2, xunits);
  gridDefYname(gridID2, yname);
  gridDefYlongname(gridID2, ylongname);
  gridDefYunits(gridID2, yunits);

  if ( gridIsRotated(gridID1) )
    {
      gridDefXpole(gridID2, gridInqXpole(gridID1));
      gridDefYpole(gridID2, gridInqYpole(gridID1));
    }

  lxvals = gridInqXvals(gridID1, NULL);
  lyvals = gridInqYvals(gridID1, NULL);

  if ( gridtype == GRID_CURVILINEAR )
    {
      if ( lxvals && lyvals )
	{
	  xvals1 = malloc(nlon1*nlat1*sizeof(double));
	  yvals1 = malloc(nlon1*nlat1*sizeof(double));
	  xvals2 = malloc(nlon2*nlat2*sizeof(double));
	  yvals2 = malloc(nlon2*nlat2*sizeof(double));
	}
    }
  else
    {
      if ( lxvals ) xvals1 = malloc(nlon1*sizeof(double));
      if ( lyvals ) yvals1 = malloc(nlat1*sizeof(double));
      if ( lxvals ) xvals2 = malloc(nlon2*sizeof(double));
      if ( lyvals ) yvals2 = malloc(nlat2*sizeof(double));
    }

  pxvals2 = xvals2;
  pyvals2 = yvals2;

  if ( xvals1 ) gridInqXvals(gridID1, xvals1);
  if ( yvals1 ) gridInqYvals(gridID1, yvals1);

  if ( gridtype == GRID_CURVILINEAR )
    {
      if ( lxvals && lyvals )
	for ( ilat = lat1; ilat <= lat2; ilat++ )
	  {
	    for ( ilon = lon21; ilon <= lon22; ilon++ )
	      {
		*pxvals2++ = xvals1[ilat*nlon1 + ilon];
		*pyvals2++ = yvals1[ilat*nlon1 + ilon];
	      }
	    for ( ilon = lon11; ilon <= lon12; ilon++ )
	      {
		*pxvals2++ = xvals1[ilat*nlon1 + ilon];
		*pyvals2++ = yvals1[ilat*nlon1 + ilon];
	      }
	  }
    }
  else
    {
      if ( lxvals )
	{
	  for ( i = lon21; i <= lon22; i++ ) *pxvals2++ = xvals1[i];
	  for ( i = lon11; i <= lon12; i++ ) *pxvals2++ = xvals1[i];
	  if ( memcmp(xunits, "degree", 6) == 0 ) correct_xvals(nlon2, 1, xvals2);
	}
      
      if ( lyvals ) for ( i = lat1;  i <= lat2;  i++ ) *pyvals2++ = yvals1[i];
    }
  /*
    for ( i = 0; i < nlat2; i++ ) printf("lat : %d %g\n", i+1, yvals2[i]);
    for ( i = 0; i < nlon2; i++ ) printf("lon : %d %g\n", i+1, xvals2[i]);
  */
  if ( xvals2 ) gridDefXvals(gridID2, xvals2);
  if ( yvals2 ) gridDefYvals(gridID2, yvals2);

  if ( xvals1 ) free(xvals1);
  if ( yvals1 ) free(yvals1);
  if ( xvals2 ) free(xvals2);
  if ( yvals2 ) free(yvals2);

  if ( gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL) )
    {
      if ( gridtype == GRID_CURVILINEAR )
	{
	  xbounds1 = malloc(4*nlon1*nlat1*sizeof(double));
	  ybounds1 = malloc(4*nlon1*nlat1*sizeof(double));
	  xbounds2 = malloc(4*nlon2*nlat2*sizeof(double));
	  ybounds2 = malloc(4*nlon2*nlat2*sizeof(double));
	}
      else
	{
	  xbounds1 = malloc(2*nlon1*sizeof(double));
	  ybounds1 = malloc(2*nlat1*sizeof(double));
	  xbounds2 = malloc(2*nlon2*sizeof(double));
	  ybounds2 = malloc(2*nlat2*sizeof(double));
	}

      pxbounds2 = xbounds2;
      pybounds2 = ybounds2;

      gridInqXbounds(gridID1, xbounds1);
      gridInqYbounds(gridID1, ybounds1);

      if ( gridtype == GRID_CURVILINEAR )
	{
	  gridDefNvertex(gridID2, 4);
	  for ( ilat = lat1; ilat <= lat2; ilat++ )
	    {
	      for ( ilon = 4*lon21; ilon < 4*(lon22+1); ilon++ )
		{
		  *pxbounds2++ = xbounds1[4*ilat*nlon1 + ilon];
		  *pybounds2++ = ybounds1[4*ilat*nlon1 + ilon];
		}
	      for ( ilon = 4*lon11; ilon < 4*(lon12+1); ilon++ )
		{
		  *pxbounds2++ = xbounds1[4*ilat*nlon1 + ilon];
		  *pybounds2++ = ybounds1[4*ilat*nlon1 + ilon];
		}
	    }
	}
      else
	{
	  gridDefNvertex(gridID2, 2);
	  for ( i = 2*lon21; i < 2*(lon22+1); i++ ) *pxbounds2++ = xbounds1[i];
	  for ( i = 2*lon11; i < 2*(lon12+1); i++ ) *pxbounds2++ = xbounds1[i];
	  for ( i = 2*lat1;  i < 2*(lat2+1);  i++ ) *pybounds2++ = ybounds1[i];

	  if ( memcmp(xunits, "degree", 6) == 0 )
	    {
	      correct_xvals(nlon2, 2, xbounds2);
	      correct_xvals(nlon2, 2, xbounds2+1);
	    }
	}

      gridDefXbounds(gridID2, xbounds2);
      gridDefYbounds(gridID2, ybounds2);

      free(xbounds1);
      free(ybounds1);
      free(xbounds2);
      free(ybounds2);
    }

  return (gridID2);
}

static
int gengridcell(int gridID1, int gridsize2, int *cellidx)
{
  int gridtype, gridID2;
  int gridsize1;
  int i, k, nv;
  int prec;
  char xname[CDI_MAX_NAME], xlongname[CDI_MAX_NAME], xunits[CDI_MAX_NAME];
  char yname[CDI_MAX_NAME], ylongname[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
  double *xvals1 = NULL, *yvals1 = NULL;
  double *xvals2 = NULL, *yvals2 = NULL;
  double *xbounds1 = NULL, *ybounds1 = NULL;
  double *xbounds2 = NULL, *ybounds2 = NULL;

  gridsize1 = gridInqSize(gridID1);

  /* printf("gridsize1 = %d, gridsize2 = %d\n", gridsize1, gridsize2); */

  gridtype = gridInqType(gridID1);
  prec     = gridInqPrec(gridID1);

  gridID2 = gridCreate(gridtype, gridsize2);

  gridDefPrec(gridID2, prec);

  gridInqXname(gridID1, xname);
  gridInqXlongname(gridID1, xlongname);
  gridInqXunits(gridID1, xunits);
  gridInqYname(gridID1, yname);
  gridInqYlongname(gridID1, ylongname);
  gridInqYunits(gridID1, yunits);

  gridDefXname(gridID2, xname);
  gridDefXlongname(gridID2, xlongname);
  gridDefXunits(gridID2, xunits);
  gridDefYname(gridID2, yname);
  gridDefYlongname(gridID2, ylongname);
  gridDefYunits(gridID2, yunits);

  if ( gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL) )
    {
      xvals1 = malloc(gridsize1*sizeof(double));
      yvals1 = malloc(gridsize1*sizeof(double));
      xvals2 = malloc(gridsize2*sizeof(double));
      yvals2 = malloc(gridsize2*sizeof(double));

      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);

      for ( i = 0; i < gridsize2; ++i )
	{
	  xvals2[i] = xvals1[cellidx[i]];
	  yvals2[i] = yvals1[cellidx[i]];
	}

      if ( cdoVerbose )
	for ( i = 0; i < gridsize2; i++ ) printf("lat/lon : %d %g %g\n", i+1, yvals2[i], xvals2[i]);

      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);

      free(xvals1);
      free(yvals1);
      free(xvals2);
      free(yvals2);
    }

  if ( gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL) )
    {
      nv = gridInqNvertex(gridID1);

      xbounds1 = malloc(nv*gridsize1*sizeof(double));
      ybounds1 = malloc(nv*gridsize1*sizeof(double));
      xbounds2 = malloc(nv*gridsize2*sizeof(double));
      ybounds2 = malloc(nv*gridsize2*sizeof(double));

      gridInqXbounds(gridID1, xbounds1);
      gridInqYbounds(gridID1, ybounds1);

      gridDefNvertex(gridID2, nv);

      for ( i = 0; i < gridsize2; ++i )
	{
	  for ( k = 0; k < nv; ++k )
	    {
	      xbounds2[i*nv+k] = xbounds1[cellidx[i]*nv+k];
	      ybounds2[i*nv+k] = ybounds1[cellidx[i]*nv+k];
	    }
	}

      gridDefXbounds(gridID2, xbounds2);
      gridDefYbounds(gridID2, ybounds2);

      free(xbounds1);
      free(ybounds1);
      free(xbounds2);
      free(ybounds2);
    }

  return (gridID2);
}


void genlonlatbox_reg(double xlon1, double xlon2, double xlat1, double xlat2,
		      int nlon1, int nlat1, double *xvals1, double *yvals1,
		      int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{
  if ( IS_NOT_EQUAL(xlon1, xlon2) )
    {
      xlon2 -= 360 * floor ((xlon2 - xlon1) / 360);
      if ( IS_EQUAL(xlon1, xlon2) ) xlon2 += 360;
    }
  else
    {
      xlon2 += 0.00001;
    }

  xlon2 -= 360 * floor ((xlon1 - xvals1[0]) / 360);
  xlon1 -= 360 * floor ((xlon1 - xvals1[0]) / 360);

  // while ( nlon1 == 1 || (xvals1[nlon1-1] - xvals1[0]) >= 360 ) nlon1--;

  for ( *lon21 = 0; *lon21 < nlon1 && xvals1[*lon21] < xlon1; (*lon21)++ );
  for ( *lon22 = *lon21; *lon22 < nlon1 && xvals1[*lon22] < xlon2; (*lon22)++ );

  if ( *lon22 >= nlon1 || xvals1[*lon22] > xlon2 ) (*lon22)--;

  xlon1 -= 360;
  xlon2 -= 360;

  for ( *lon11 = 0; xvals1[*lon11] < xlon1; (*lon11)++ );
  for ( *lon12 = *lon11; *lon12 < nlon1 && xvals1[*lon12] < xlon2; (*lon12)++ );
  
  // (*lon12)--;
  if ( *lon12 >= nlon1 || xvals1[*lon12] > xlon2 ) (*lon12)--;
  if ( *lon12 >= 0 )
    if ( IS_EQUAL(xvals1[*lon12], xvals1[*lon21]) ) (*lon12)--;

  if ( *lon12 - *lon11 + 1 + *lon22 - *lon21 + 1 <= 0 )
    cdoAbort("Longitudinal dimension is too small!");
  
  if ( yvals1[0] > yvals1[nlat1 - 1] )
    {
      if ( xlat1 > xlat2 )
	{
	  for ( *lat1 = 0; *lat1 < nlat1 && yvals1[*lat1] > xlat1; (*lat1)++ );
	  for ( *lat2 = nlat1 - 1; *lat2 && yvals1[*lat2] < xlat2; (*lat2)-- );
	}
      else
	{
	  for ( *lat1 = 0; *lat1 < nlat1 && yvals1[*lat1] > xlat2; (*lat1)++ );
	  for ( *lat2 = nlat1 - 1; *lat2 && yvals1[*lat2] < xlat1; (*lat2)-- );
	}
    }
  else
    {
      if ( xlat1 < xlat2 )
	{
	  for ( *lat1 = 0; *lat1 < nlat1 && yvals1[*lat1] < xlat1; (*lat1)++ );
	  for ( *lat2 = nlat1 - 1; *lat2 && yvals1[*lat2] > xlat2; (*lat2)-- );
	}
      else
	{
	  for ( *lat1 = 0; *lat1 < nlat1 && yvals1[*lat1] < xlat2; (*lat1)++ );
	  for ( *lat2 = nlat1 - 1; *lat2 && yvals1[*lat2] > xlat1; (*lat2)-- );
	}
    }
  
  if ( *lat2 - *lat1 + 1 <= 0 )
    cdoAbort("Latitudinal dimension is too small!");
}


void genlonlatbox(int argc_offset, int gridID1, int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{
  int ilon, ilat;
  int nlon1, nlat1;
  int gridtype;
  double *xvals1, *yvals1;
  double xlon1, xlon2, xlat1, xlat2;
  int grid_is_circular;
  char xunits[CDI_MAX_NAME];
  char yunits[CDI_MAX_NAME];
  double xfact = 1, yfact = 1;

  operatorCheckArgc(argc_offset+4);

  xlon1 = atof(operatorArgv()[argc_offset+0]);
  xlon2 = atof(operatorArgv()[argc_offset+1]);
  xlat1 = atof(operatorArgv()[argc_offset+2]);
  xlat2 = atof(operatorArgv()[argc_offset+3]);

  gridtype = gridInqType(gridID1);

  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  grid_is_circular = gridIsCircular(gridID1);

  if ( gridtype == GRID_CURVILINEAR )
    {
      xvals1 = malloc(nlon1*nlat1*sizeof(double));
      yvals1 = malloc(nlon1*nlat1*sizeof(double));
    }
  else
    {
      xvals1 = malloc(nlon1*sizeof(double));
      yvals1 = malloc(nlat1*sizeof(double));
    }

  gridInqXvals(gridID1, xvals1);
  gridInqYvals(gridID1, yvals1);

  gridInqXunits(gridID1, xunits);
  gridInqYunits(gridID1, yunits);

  if ( memcmp(xunits, "radian", 6) == 0 ) xfact = RAD2DEG;
  if ( memcmp(yunits, "radian", 6) == 0 ) yfact = RAD2DEG;

  if ( gridtype == GRID_CURVILINEAR )
    {
      double xval, yval, xfirst, xlast, ylast;
      int lp2 = FALSE;

      if ( xlon1 > xlon2 ) 
	cdoAbort("The second longitude have to be greater than the first one!");

      if ( xlat1 > xlat2 )
	{
	  double xtemp = xlat1;
	  xlat1 = xlat2;
	  xlat2 = xtemp;
	}
	  
      *lat1 = nlat1-1;
      *lat2 = 0;
      *lon11 = 0;
      *lon12 = -1;
      *lon21 = nlon1-1;
      *lon22 = 0;

      for ( ilat = 0; ilat < nlat1; ilat++ )
	{
	  xlast = xfact * xvals1[ilat*nlon1 + nlon1-1];
	  ylast = yfact * yvals1[ilat*nlon1 + nlon1-1];
	  if ( ylast >= xlat1 && ylast <= xlat2 )
	    if ( grid_is_circular && xlon1 <= xlast && xlon2 > xlast && (xlon2-xlon1) < 360 )
	      {
		*lon11 = nlon1-1;
		*lon12 = 0;
		lp2 = TRUE;
	      }
	}

      for ( ilat = 0; ilat < nlat1; ilat++ )
	{
	  for ( ilon = 0; ilon < nlon1; ilon++ )
	    {
	      xval = xvals1[ilat*nlon1 + ilon];
	      yval = yvals1[ilat*nlon1 + ilon];

	      xval *= xfact;
	      yval *= yfact;

	      if ( yval >= xlat1 && yval <= xlat2 )
		{
		  if ( lp2 )
		    {
		      xfirst = xfact * xvals1[ilat*nlon1];
		      if ( xfirst < xlon1 ) xfirst = xlon1;

		      xlast = xfact * xvals1[ilat*nlon1 + nlon1-1];
		      if ( xlast > xlon2 ) xlast = xlon2;

		      if ( xval >= xlon1 && xval <= xlast )
			{
			  if ( ilon < *lon21 ) *lon21 = ilon;
			  if ( ilon > *lon22 ) *lon22 = ilon;
			  if ( ilat < *lat1 ) *lat1 = ilat;
			  if ( ilat > *lat2 ) *lat2 = ilat;
			}
		      else if ( xval >= xfirst && xval <= xlon2 )
			{
			  if ( ilon < *lon11 ) *lon11 = ilon;
			  if ( ilon > *lon12 ) *lon12 = ilon;
			  if ( ilat < *lat1 ) *lat1 = ilat;
			  if ( ilat > *lat2 ) *lat2 = ilat;
			}
		    }
		  else
		    {
		      if ( ((xval     >= xlon1 && xval     <= xlon2) ||
			    (xval-360 >= xlon1 && xval-360 <= xlon2) ||
			    (xval+360 >= xlon1 && xval+360 <= xlon2)) )
			{
			  if ( ilon < *lon21 ) *lon21 = ilon;
			  if ( ilon > *lon22 ) *lon22 = ilon;
			  if ( ilat < *lat1 ) *lat1 = ilat;
			  if ( ilat > *lat2 ) *lat2 = ilat;
			}
		    }
		}
	    }
	}

      // printf("lon11, lon12, lon21, lon22, lat1, lat2 %d %d %d %d %d %d\n", *lon11, *lon12, *lon21, *lon22, *lat1, *lat2);
      if ( *lon12 == 0 && *lon11 > 0 ) *lon11 = -1;

      if ( *lat2 - *lat1 + 1 <= 0 )
	cdoAbort("Latitudinal dimension is too small!");
    }
  else
    {
      for ( ilat = 0; ilat < nlat1; ilat++ ) yvals1[ilat] *= yfact;
      for ( ilon = 0; ilon < nlon1; ilon++ ) xvals1[ilon] *= xfact;

      genlonlatbox_reg(xlon1, xlon2, xlat1, xlat2,
		       nlon1, nlat1, xvals1, yvals1,
		       lat1, lat2, lon11, lon12, lon21, lon22);
    }

  free(xvals1);
  free(yvals1);
}

static
int genlonlatgrid(int gridID1, int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{
  int gridID2;

  genlonlatbox(0, gridID1, lat1, lat2, lon11, lon12, lon21, lon22);

  gridID2 = gengrid(gridID1, *lat1, *lat2, *lon11, *lon12, *lon21, *lon22);

  return (gridID2);
}

static
int gencellgrid(int gridID1, int *gridsize2, int **cellidx)
{
  int gridtype, gridID2;
  double *xvals1, *yvals1;
  double xlon1, xlon2, xlat1, xlat2, x, xval, yval;
  int i, gridsize1;
  int nvals = 0;
  int maxcell = 0;
  int cellinc = 4096;
  char xunits[CDI_MAX_NAME];
  char yunits[CDI_MAX_NAME];
  double xfact, yfact;
  int argc_offset = 0;

  operatorCheckArgc(argc_offset+4);

  xlon1 = atof(operatorArgv()[argc_offset+0]);
  xlon2 = atof(operatorArgv()[argc_offset+1]);
  xlat1 = atof(operatorArgv()[argc_offset+2]);
  xlat2 = atof(operatorArgv()[argc_offset+3]);

  if ( xlon1 >= xlon2 ) { x = xlon1; xlon1 = xlon2; xlon2 = x; }
  if ( xlat1 >= xlat2 ) { x = xlat1; xlat1 = xlat2; xlat2 = x; }

  gridtype = gridInqType(gridID1);

  gridsize1 = gridInqSize(gridID1);

  if ( gridtype != GRID_UNSTRUCTURED ) cdoAbort("Internal problem, wrong grid type!");

  xvals1 = malloc(gridsize1*sizeof(double));
  yvals1 = malloc(gridsize1*sizeof(double));

  gridInqXvals(gridID1, xvals1);
  gridInqYvals(gridID1, yvals1);

  gridInqXunits(gridID1, xunits);
  gridInqYunits(gridID1, yunits);

  if ( memcmp(xunits, "radian", 6) == 0 )
    xfact = RAD2DEG;
  else
    xfact = 1;

  if ( memcmp(yunits, "radian", 6) == 0 )
    yfact = RAD2DEG;
  else
    yfact = 1;

  /* find gridsize2 */
  *cellidx = NULL;
  for ( i = 0; i < gridsize1; ++i )
    {
      xval = xvals1[i]*xfact;
      yval = yvals1[i]*yfact;
      if ( yval >= xlat1 && yval <= xlat2 )
	if ( (xval >= xlon1 && xval <= xlon2) ||
	     (xval+360 >= xlon1 && xval+360 <= xlon2) ||
	     (xval-360 >= xlon1 && xval-360 <= xlon2)  )
	  {
	    nvals++;
	    if ( nvals > maxcell )
	      {
		maxcell += cellinc;
		*cellidx = realloc(*cellidx, maxcell*sizeof(int));
	      }
	    (*cellidx)[nvals-1] = i;
	  }
    }

  if ( nvals == 0 ) cdoAbort("No grid points found!");

  *gridsize2 = nvals;

  free(xvals1);
  free(yvals1);

  gridID2 = gengridcell(gridID1, *gridsize2, *cellidx);

  return (gridID2);
}


void genindexbox(int argc_offset, int gridID1, int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{
  int nlon1, nlat1;
  int temp;

  operatorCheckArgc(argc_offset+4);

  *lon11 = atoi(operatorArgv()[argc_offset+0]);
  *lon12 = atoi(operatorArgv()[argc_offset+1]);
  *lat1  = atoi(operatorArgv()[argc_offset+2]);
  *lat2  = atoi(operatorArgv()[argc_offset+3]);

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
      cdoWarning("First latitude index out of range, set to 1!");
      *lat1 = 1;
    }
  if ( *lat2 > nlat1 )
    {
      cdoWarning("First latitude index out of range, set to %d!", nlat1);
      *lat1 = nlat1;
    }
  if ( *lat2 < 1 )
    {
      cdoWarning("First latitude index out of range, set to 1!");
      *lat2 = 1;
    }
  if ( *lat2 > nlat1 )
    {
      cdoWarning("Last latitude index out of range, set to %d!", nlat1);
      *lat2 = nlat1;
    }
  if ( *lon11 < 1 )
    {
      cdoWarning("First longitude index out of range, set to 1!");
      *lon11 = 1;
    }
  if ( *lon12 > nlon1+1 )
    {
      cdoWarning("Last longitude index out of range, set to %d!", nlon1);
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
int genindexgrid(int gridID1, int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{
  int gridID2;

  genindexbox(0, gridID1, lat1, lat2, lon11, lon12, lon21, lon22);

  gridID2 = gengrid(gridID1, *lat1, *lat2, *lon11, *lon12, *lon21, *lon22);

  return (gridID2);
}

static
void window(int nwpv, double *array1, int gridID1, double *array2,
	    long lat1, long lat2, long lon11, long lon12, long lon21, long lon22)
{
  long nlon1;
  long ilat, ilon;

  nlon1 = gridInqXsize(gridID1);

  if ( nwpv == 2 )
    {
      for ( ilat = lat1; ilat <= lat2; ilat++ )
	{
	  for ( ilon = lon21; ilon <= lon22; ilon++ )
	    {
	      *array2++ = array1[ilat*nlon1*2 + ilon*2];
	      *array2++ = array1[ilat*nlon1*2 + ilon*2+1];
	    }
	  for ( ilon = lon11; ilon <= lon12; ilon++ )
	    {
	      *array2++ = array1[ilat*nlon1*2 + ilon*2];
	      *array2++ = array1[ilat*nlon1*2 + ilon*2+1];
	    }
	}
    }
  else
    {
      for ( ilat = lat1; ilat <= lat2; ilat++ )
	{
	  for ( ilon = lon21; ilon <= lon22; ilon++ )
	    *array2++ = array1[ilat*nlon1 + ilon];
	  for ( ilon = lon11; ilon <= lon12; ilon++ )
	    *array2++ = array1[ilat*nlon1 + ilon];
	}
    }
}

static
void window_cell(int nwpv, double *array1, int gridID1, double *array2, long gridsize2, int *cellidx)
{
  long i;

  if ( nwpv == 2 )
    {
      for ( i = 0; i < gridsize2; ++i )
	{
	  array2[i*2]   = array1[cellidx[i]*2];
	  array2[i*2+1] = array1[cellidx[i]*2+1];
	}
    }
  else
    {
      for ( i = 0; i < gridsize2; ++i )
	array2[i] = array1[cellidx[i]];
    }
}


void *Selbox(void *argument)
{
  int SELLONLATBOX, SELINDEXBOX;
  int operatorID;
  int streamID1, streamID2;
  int nrecs, nvars;
  int tsID, recID, varID, levelID;
  int gridsize, gridsize2;
  int vlistID1, vlistID2;
  int gridID1 = -1, gridID2;
  int index, ngrids, gridtype = -1;
  int nmiss;
  int *vars = NULL;
  int i;
  int nwpv; // number of words per value; real:1  complex:2
  double missval;
  double *array1 = NULL, *array2 = NULL;
  int taxisID1, taxisID2;
  typedef struct {
    int gridID1, gridID2;
    int *cellidx, nvals;
    int lat1, lat2, lon11, lon12, lon21, lon22; 
  } sbox_t;
  sbox_t *sbox = NULL;

  cdoInitialize(argument);

  SELLONLATBOX = cdoOperatorAdd("sellonlatbox", 0, 0, "western and eastern longitude and southern and northern latitude");
  SELINDEXBOX  = cdoOperatorAdd("selindexbox",  0, 0, "index of first and last longitude and index of first and last latitude");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nvars = vlistNvars(vlistID1);
  vars  = malloc(nvars*sizeof(int));
  for ( varID = 0; varID < nvars; varID++ ) vars[varID] = FALSE;

  ngrids = vlistNgrids(vlistID1);
  sbox = malloc(ngrids*sizeof(sbox_t));

  for ( index = 0; index < ngrids; index++ )
    {
      gridID1  = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);

      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR ||
	   (operatorID == SELINDEXBOX && (gridtype == GRID_GENERIC || gridtype == GRID_SINUSOIDAL) && 
	    gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0) ||
	   (operatorID == SELLONLATBOX && gridtype == GRID_UNSTRUCTURED) )
	{
	  if ( operatorID == SELLONLATBOX )
	    {
	      if ( gridtype == GRID_UNSTRUCTURED )
		gridID2 = gencellgrid(gridID1, &sbox[index].nvals, &sbox[index].cellidx);
	      else
		gridID2 = genlonlatgrid(gridID1, &sbox[index].lat1, &sbox[index].lat2, &sbox[index].lon11, 
					&sbox[index].lon12, &sbox[index].lon21, &sbox[index].lon22);
	    }
	  else
	    gridID2 = genindexgrid(gridID1, &sbox[index].lat1, &sbox[index].lat2, &sbox[index].lon11, 
				   &sbox[index].lon12, &sbox[index].lon21, &sbox[index].lon22);
	  
	  sbox[index].gridID1 = gridID1;
	  sbox[index].gridID2 = gridID2;

	  vlistChangeGridIndex(vlistID2, index, gridID2);

	  for ( varID = 0; varID < nvars; varID++ )
	    if ( gridID1 == vlistInqVarGrid(vlistID1, varID) )
	      vars[varID] = TRUE;
	}
      else if ( gridtype == GRID_GENERIC && gridInqXsize(gridID1) <= 1 && gridInqYsize(gridID1) <=1 )
	{
	}
      else
	{
	  cdoPrint("Unsupported grid type: %s", gridNamePtr(gridtype));
	  if ( gridtype == GRID_GAUSSIAN_REDUCED )
	    cdoPrint("Use option -R to convert Gaussian reduced grid to a regular grid!");
	  cdoAbort("Unsupported grid type!");
	}
    }

  if ( cdoVerbose )
    {
      if ( gridtype != GRID_UNSTRUCTURED )
	{
	  cdoPrint("box1 - idx1,idx2,idy1,idy2: %d,%d,%d,%d",
		   sbox[0].lon21+1, sbox[0].lon22+1, sbox[0].lat1+1, sbox[0].lat2+1);
	  cdoPrint("box2 - idx1,idx2,idy1,idy2: %d,%d,%d,%d",
		   sbox[0].lon11+1, sbox[0].lon12+1, sbox[0].lat1+1, sbox[0].lat2+1);
	}
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  array1 = malloc(gridsize*sizeof(double));

  gridsize2 = vlistGridsizeMax(vlistID2);
  if ( vlistNumber(vlistID2) != CDI_REAL ) gridsize2 *= 2;
  array2 = malloc(gridsize2*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  streamDefRecord(streamID2, varID, levelID);

	  if ( vars[varID] )
	    {
	      nwpv    = vlistInqNWPV(vlistID1, varID);
	      
	      gridID1 = vlistInqVarGrid(vlistID1, varID);

	      for ( index = 0; index < ngrids; index++ )
		if ( gridID1 == sbox[index].gridID1 ) break;

	      if ( index == ngrids ) cdoAbort("Internal problem, grid not found!");

	      gridsize2 = gridInqSize(sbox[index].gridID2);

	      if ( operatorID == SELLONLATBOX && gridtype == GRID_UNSTRUCTURED )
		window_cell(nwpv, array1, gridID1, array2, gridsize2, sbox[index].cellidx);
 	      else
		window(nwpv, array1, gridID1, array2, sbox[index].lat1, sbox[index].lat2, sbox[index].lon11, 
		       sbox[index].lon12, sbox[index].lon21, sbox[index].lon22);

	      if ( nmiss )
		{
		  nmiss = 0;
		  missval = vlistInqVarMissval(vlistID2, varID);
		  for ( i = 0; i < gridsize2; i++ )
		    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
		}

	      streamWriteRecord(streamID2, array2, nmiss);
	    }
	  else
	    {
	      streamWriteRecord(streamID2, array1, nmiss);
	    }
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  if ( vars   ) free(vars);
  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  if ( sbox )
    {
      if ( operatorID == SELLONLATBOX  && gridtype == GRID_UNSTRUCTURED )
	{
	  for ( index = 0; index < ngrids; index++ )
	    if ( sbox[index].cellidx ) free(sbox[index].cellidx);
	}
      free(sbox);
    }

  cdoFinish();

  return (0);
}
