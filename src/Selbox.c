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

      Selbox     sellonlatbox    Select lon/lat box
      Selbox     selindexbox     Select index box
*/

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


static
int gengrid(int gridID1, int lat1, int lat2, int lon11, int lon12, int lon21, int lon22)
{
  static char func[] = "gengrid";  
  int gridtype, gridID2;
  int nlon1, nlat1;
  int nlon2, nlat2;
  int nlon21, nlon22;
  int i;
  int prec;
  int ilat, ilon;
  char xname[128], xlongname[128], xunits[128];
  char yname[128], ylongname[128], yunits[128];
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

  if ( gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL) )
    {
      if ( gridtype == GRID_CURVILINEAR )
	{
	  xvals1 = (double *) malloc(nlon1*nlat1*sizeof(double));
	  yvals1 = (double *) malloc(nlon1*nlat1*sizeof(double));
	  xvals2 = (double *) malloc(nlon2*nlat2*sizeof(double));
	  yvals2 = (double *) malloc(nlon2*nlat2*sizeof(double));
	}
      else
	{
	  xvals1 = (double *) malloc(nlon1*sizeof(double));
	  yvals1 = (double *) malloc(nlat1*sizeof(double));
	  xvals2 = (double *) malloc(nlon2*sizeof(double));
	  yvals2 = (double *) malloc(nlat2*sizeof(double));
	}

      pxvals2 = xvals2;
      pyvals2 = yvals2;

      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);

      if ( gridtype == GRID_CURVILINEAR )
	{
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
	  for ( i = lon21; i <= lon22; i++ ) *pxvals2++ = xvals1[i];
	  for ( i = lon11; i <= lon12; i++ ) *pxvals2++ = xvals1[i];
	  for ( i = lat1;  i <= lat2;  i++ ) *pyvals2++ = yvals1[i];

	  if ( memcmp(xunits, "degree", 6) == 0 )
	    {
	      if ( IS_EQUAL(xvals2[0], xvals2[nlon2-1]) ) xvals2[nlon2-1] += 360;

	      if ( xvals2[0] > xvals2[nlon2-1] )
		for ( i = 0; i < nlon2; i++ )
		  if ( xvals2[i] >= 180 ) xvals2[i] -= 360;

	      for ( i = 0; i < nlon2; i++ )
		{
		  if ( xvals2[i] < -180 ) xvals2[i] += 360;
		  if ( xvals2[i] >  360 ) xvals2[i] -= 360;
		}

	      for ( i = 1; i < nlon2; i++ )
		{
		  if ( xvals2[i] < xvals2[i-1] ) xvals2[i] += 360;
		}
	    }
	}
      /*
      for ( i = 0; i < nlat2; i++ ) printf("lat : %d %g\n", i+1, yvals2[i]);
      for ( i = 0; i < nlon2; i++ ) printf("lon : %d %g\n", i+1, xvals2[i]);
      */
      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);

      free(xvals1);
      free(yvals1);
      free(xvals2);
      free(yvals2);
    }

  if ( gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL) )
    {
      if ( gridtype == GRID_CURVILINEAR )
	{
	  xbounds1 = (double *) malloc(4*nlon1*nlat1*sizeof(double));
	  ybounds1 = (double *) malloc(4*nlon1*nlat1*sizeof(double));
	  xbounds2 = (double *) malloc(4*nlon2*nlat2*sizeof(double));
	  ybounds2 = (double *) malloc(4*nlon2*nlat2*sizeof(double));
	}
      else
	{
	  xbounds1 = (double *) malloc(2*nlon1*sizeof(double));
	  ybounds1 = (double *) malloc(2*nlat1*sizeof(double));
	  xbounds2 = (double *) malloc(2*nlon2*sizeof(double));
	  ybounds2 = (double *) malloc(2*nlat2*sizeof(double));
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
  static char func[] = "gengridcell";  
  int gridtype, gridID2;
  int gridsize1;
  int i, k, nv;
  int prec;
  char xname[128], xlongname[128], xunits[128];
  char yname[128], ylongname[128], yunits[128];
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
      xvals1 = (double *) malloc(gridsize1*sizeof(double));
      yvals1 = (double *) malloc(gridsize1*sizeof(double));
      xvals2 = (double *) malloc(gridsize2*sizeof(double));
      yvals2 = (double *) malloc(gridsize2*sizeof(double));

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

      xbounds1 = (double *) malloc(nv*gridsize1*sizeof(double));
      ybounds1 = (double *) malloc(nv*gridsize1*sizeof(double));
      xbounds2 = (double *) malloc(nv*gridsize2*sizeof(double));
      ybounds2 = (double *) malloc(nv*gridsize2*sizeof(double));

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

static
int genlonlatgrid(int gridID1, int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{
  static char func[] = "genlonlatgrid";  
  int ilon, ilat;
  int nlon1, nlat1;
  int gridtype, gridID2;
  double *xvals1, *yvals1;
  double xlon1, xlon2, xlat1, xlat2;
  int grid_is_circular;
  char xunits[128];
  char yunits[128];
  double xfact, yfact;

  operatorCheckArgc(4);

  xlon1 = atof(operatorArgv()[0]);
  xlon2 = atof(operatorArgv()[1]);
  xlat1 = atof(operatorArgv()[2]);
  xlat2 = atof(operatorArgv()[3]);

  gridtype = gridInqType(gridID1);

  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  grid_is_circular = gridIsCircular(gridID1);

  if ( gridtype == GRID_CURVILINEAR )
    {
      xvals1 = (double *) malloc(nlon1*nlat1*sizeof(double));
      yvals1 = (double *) malloc(nlon1*nlat1*sizeof(double));
    }
  else
    {
      xvals1 = (double *) malloc(nlon1*sizeof(double));
      yvals1 = (double *) malloc(nlat1*sizeof(double));
    }

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

  if ( gridtype == GRID_CURVILINEAR )
    {
      double xval, yval, xlast, ylast;
      int ixmin, ixmax, iymin, iymax;
      int lp2 = FALSE;

      if ( xlon1 > xlon2 ) 
	cdoAbort("The second longitude have to be greater than the first one!");
      if ( xlat1 > xlat2 )
	cdoAbort("The second latitude have to be greater than the first one!");

      *lat1 = nlat1-1;
      *lat2 = 0;
      *lon11 = 0;
      *lon12 = -1;
      *lon21 = nlon1-1;
      *lon22 = 0;

      ixmin = nlon1-1;
      ixmax = 0;
      iymin = nlat1-1;
      iymax = 0;

      for ( ilat = 0; ilat < nlat1; ilat++ )
	{
	  xlast = xfact * xvals1[ilat*nlon1 + nlon1-1];
	  ylast = yfact * yvals1[ilat*nlon1 + nlon1-1];
	  if ( ylast >= xlat1 && ylast <= xlat2 )
	    if ( grid_is_circular && xlon1 <= xlast && xlon2 > xlast && (xlon2-xlon1) < 360)
	      {
		*lon11 = nlon1-1;
		*lon12 = 0;
		lp2 = TRUE;
	      }
	}

      for ( ilat = 0; ilat < nlat1; ilat++ )
	{
	  xlast = xfact * xvals1[ilat*nlon1 + nlon1-1];

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
		      if ( xval >= xlon1 && xval <= xlast )
			{
			  if ( ilon < *lon21 ) *lon21 = ilon;
			  if ( ilon > *lon22 ) *lon22 = ilon;
			  if ( ilat < *lat1 ) *lat1 = ilat;
			  if ( ilat > *lat2 ) *lat2 = ilat;
			}
		      else if ( xval > xlast && xval <= xlon2 )
			{
			  if ( ilon < *lon11 ) *lon11 = ilon;
			  if ( ilon > *lon12 ) *lon12 = ilon;
			  if ( ilat < *lat1 ) *lat1 = ilat;
			  if ( ilat > *lat2 ) *lat2 = ilat;
			}
		    }
		  else
		    {
		      if ( ((xval >= xlon1 && xval <= xlon2) ||
			    (xval-360 >= xlon1 && xval-360 <= xlon2) ||
			    (xval+360 >= xlon1 && xval+360 <= xlon2)) )
			{
			  if ( ilon < ixmin ) ixmin = ilon;
			  if ( ilon > ixmax ) ixmax = ilon;
			  if ( ilat < iymin ) iymin = ilat;
			  if ( ilat > iymax ) iymax = ilat;
			}
		    }
		}
	    }
	}

      if ( ! lp2 )
	{
	  *lat1 = iymin;
	  *lat2 = iymax;
	  *lon21 = ixmin;
	  *lon22 = ixmax;
	}
      
      if ( cdoVerbose )
	fprintf(stderr, " ixmin=%d, ixmax=%d, iymin=%d, iymax=%d\n", *lon21, *lon22, *lat1, *lat2);
      
    }
  else
    {
      for ( ilat = 0; ilat < nlat1; ilat++ ) yvals1[ilat] *= yfact;
      for ( ilon = 0; ilon < nlon1; ilon++ ) xvals1[ilon] *= xfact;

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

      for ( *lon21 = 0; *lon21 < nlon1 && xvals1[*lon21] < xlon1; (*lon21)++ );
      for ( *lon22 = *lon21; *lon22 < nlon1 && xvals1[*lon22] < xlon2; (*lon22)++ );

      if ( *lon22 >= nlon1 ) (*lon22)--;

      xlon1 -= 360;
      xlon2 -= 360;

      for ( *lon11 = 0; xvals1[*lon11] < xlon1; (*lon11)++ );
      for ( *lon12 = *lon11; *lon12 < nlon1 && xvals1[*lon12] < xlon2; (*lon12)++ );

      (*lon12)--;

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
    }

  if ( *lat2 - *lat1 + 1 <= 0 )
    cdoAbort("Latitudinal dimension is too small!");

  /*
  if ( (*lon22 - *lon21 + 1 <= 0)  && (*lon12 - *lon11 + 1 <= 0)  )
    cdoAbort("Longitudinal dimension is too small!");
  */
  free(xvals1);
  free(yvals1);

  gridID2 = gengrid(gridID1, *lat1, *lat2, *lon11, *lon12, *lon21, *lon22);

  return (gridID2);
}


static
int gencellgrid(int gridID1, int *gridsize2, int **cellidx)
{
  static char func[] = "gencellgrid";  
  int gridtype, gridID2;
  double *xvals1, *yvals1;
  double xlon1, xlon2, xlat1, xlat2, x, xval, yval;
  int i, gridsize1;
  int nvals = 0;
  int maxcell = 0;
  int cellinc = 4096;
  char xunits[128];
  char yunits[128];
  double xfact, yfact;

  operatorCheckArgc(4);

  xlon1 = atof(operatorArgv()[0]);
  xlon2 = atof(operatorArgv()[1]);
  xlat1 = atof(operatorArgv()[2]);
  xlat2 = atof(operatorArgv()[3]);

  if ( xlon1 >= xlon2 ) { x = xlon1; xlon1 = xlon2; xlon2 = x; }
  if ( xlat1 >= xlat2 ) { x = xlat1; xlat1 = xlat2; xlat2 = x; }

  gridtype = gridInqType(gridID1);

  gridsize1 = gridInqSize(gridID1);

  if ( gridtype != GRID_CELL ) cdoAbort("Internal problem, wrong grid type!");

  xvals1 = (double *) malloc(gridsize1*sizeof(double));
  yvals1 = (double *) malloc(gridsize1*sizeof(double));

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
		*cellidx = (int *) realloc(*cellidx, maxcell*sizeof(int));
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


static
int genindexgrid(int gridID1, int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{
  int gridID2;
  int nlon1, nlat1;
  int temp;

  operatorCheckArgc(4);

  *lon11 = atoi(operatorArgv()[0]);
  *lon12 = atoi(operatorArgv()[1]);
  *lat1  = atoi(operatorArgv()[2]);
  *lat2  = atoi(operatorArgv()[3]);

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

  gridID2 = gengrid(gridID1, *lat1, *lat2, *lon11, *lon12, *lon21, *lon22);

  return (gridID2);
}


static
void window(double *array1, int gridID1, double *array2,
	    int lat1, int lat2, int lon11, int lon12, int lon21, int lon22)
{
  int nlon1;
  int ilat, ilon;

  nlon1 = gridInqXsize(gridID1);

  for ( ilat = lat1; ilat <= lat2; ilat++ )
    {
      for ( ilon = lon21; ilon <= lon22; ilon++ )
	*array2++ = array1[ilat*nlon1 + ilon];
      for ( ilon = lon11; ilon <= lon12; ilon++ )
	*array2++ = array1[ilat*nlon1 + ilon];
    }
}

static
void window_cell(double *array1, int gridID1, double *array2, int gridsize2, int *cellidx)
{
  int i;

  for ( i = 0; i < gridsize2; ++i )
    array2[i] = array1[cellidx[i]];
}


void *Selbox(void *argument)
{
  static char func[] = "Selbox";
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
  int *vars;
  int i;
  double missval;
  double *array1 = NULL, *array2 = NULL;
  int taxisID1, taxisID2;
  typedef struct {
    int gridID1, gridID2;
    int *cellidx, nvals;
    int lat1, lat2, lon11, lon12, lon21, lon22; 
  } SBOX;
  SBOX *sbox = NULL;

  cdoInitialize(argument);

  SELLONLATBOX = cdoOperatorAdd("sellonlatbox", 0, 0, "western and eastern longitude and southern and northern latitude");
  SELINDEXBOX  = cdoOperatorAdd("selindexbox",  0, 0, "index of first and last longitude and index of first and last latitude");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nvars = vlistNvars(vlistID1);
  vars  = (int *) malloc(nvars*sizeof(int));
  for ( varID = 0; varID < nvars; varID++ ) vars[varID] = FALSE;

  ngrids = vlistNgrids(vlistID1);
  sbox = (SBOX *) malloc(ngrids*sizeof(SBOX));

  for ( index = 0; index < ngrids; index++ )
    {
      gridID1  = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);

      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR ||
	   (operatorID == SELINDEXBOX && (gridtype == GRID_GENERIC || gridtype == GRID_SINUSOIDAL) && 
	    gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0 ) ||
	   (operatorID == SELLONLATBOX && gridtype == GRID_CELL) )
	{
	  if ( operatorID == SELLONLATBOX )
	    {
	      if ( gridtype == GRID_CELL )
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
      else
	{
	  cdoPrint("Unsupported grid type: %s", gridNamePtr(gridtype));
	  if ( gridtype == GRID_GAUSSIAN_REDUCED )
	    cdoPrint("Use option -R to convert Gaussian reduced grid to a regular grid!");
	  cdoAbort("Unsupported grid type!");
	}
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  array1 = (double *) malloc(gridsize*sizeof(double));

  gridsize2 = vlistGridsizeMax(vlistID2);
  array2 = (double *) malloc(gridsize2*sizeof(double));

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
	      streamReadRecord(streamID1, array1, &nmiss);

	      gridID1 = vlistInqVarGrid(vlistID1, varID);
	      for ( index = 0; index < ngrids; index++ )
		if ( gridID1 == sbox[index].gridID1 ) break;

	      if ( index == ngrids ) cdoAbort("Internal problem, grid not found!");

	      gridsize2 = gridInqSize(sbox[index].gridID2);

	      if ( operatorID == SELLONLATBOX  && gridtype == GRID_CELL )
		window_cell(array1, gridID1, array2, gridsize2, sbox[index].cellidx);
	      else
		window(array1, gridID1, array2, sbox[index].lat1, sbox[index].lat2, sbox[index].lon11, 
		       sbox[index].lon12, sbox[index].lon21, sbox[index].lon22);

	      if ( nmiss )
		{
		  nmiss = 0;
		  missval = vlistInqVarMissval(vlistID2, varID);
		  for ( i = 0; i < gridsize2; i++ )
		    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
		}

	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, array2, nmiss);
	    }
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( vars   ) free(vars);
  if ( array2 ) free(array2);
  if ( array1 ) free(array1);
  free(sbox);
  if ( operatorID == SELLONLATBOX  && gridtype == GRID_CELL ) free(sbox[index].cellidx);

  cdoFinish();

  return (0);
}
