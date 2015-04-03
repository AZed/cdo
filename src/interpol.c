#include <string.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"


double intlinarr2p(int nxm, int nym, double **fieldm, const double *xm, const double *ym,
		   double x, double y)
{
  int ii, jj;
  double value = 0;

  for ( jj = 1; jj < nym; jj++ )
    {
      if ( y < MIN(ym[jj-1], ym[jj]) || 
	   y > MAX(ym[jj-1], ym[jj]) ) continue;
      for ( ii = 1; ii < nxm; ii++ )
	{
	  if ( x < xm[ii-1] || x > xm[ii] ) continue;
	  value = fieldm[jj-1][ii-1] * (x-xm[ii]) * (y-ym[jj])
	              / ((xm[ii-1]-xm[ii]) * (ym[jj-1]-ym[jj]))
	        + fieldm[jj-1][ii] * (x-xm[ii-1]) * (y-ym[jj])
                      / ((xm[ii]-xm[ii-1]) * (ym[jj-1]-ym[jj]))
                + fieldm[jj][ii-1] * (x-xm[ii]) * (y-ym[jj-1])
                      / ((xm[ii-1]-xm[ii]) * (ym[jj]-ym[jj-1]))
                + fieldm[jj][ii] * (x-xm[ii-1]) * (y-ym[jj-1])
	              / ((xm[ii]-xm[ii-1]) * (ym[jj]-ym[jj-1]));
	}
    }

  return value;
}


void intlinarr2(double missval,
		int nxm, int nym,  double **fieldm, const double *xm, const double *ym,
		int nx, int ny, double **field, const double *x, const double *y)
{
  int i, ii, j , jj;
  double ymin, ymax;

  for ( j = 0; j < ny; j++ )
    for ( i = 0; i < nx; i++ )
      field[j][i] = missval;

  for ( jj = 1; jj < nym; jj++ )
    {
      ymin = MIN(ym[jj-1], ym[jj]);
      ymax = MAX(ym[jj-1], ym[jj]);
      for ( j = 0; j < ny; j++ )
	{
	  if ( y[j] < ymin || y[j] > ymax ) continue;
	  for ( ii = 1; ii < nxm; ii++ )
	    {
#if defined (SX)
#pragma vdir nodep
#endif
	      for ( i = 0; i < nx; i++ )
		{
		  if ( x[i] < xm[ii-1] || x[i] > xm[ii] ) continue;
		  field[j][i] = fieldm[jj-1][ii-1] * (x[i]-xm[ii]) * (y[j]-ym[jj])
		                          / ((xm[ii-1]-xm[ii]) * (ym[jj-1]-ym[jj]))
		              + fieldm[jj-1][ii] * (x[i]-xm[ii-1]) * (y[j]-ym[jj])
                                          / ((xm[ii]-xm[ii-1]) * (ym[jj-1]-ym[jj]))
                              + fieldm[jj][ii-1] * (x[i]-xm[ii]) * (y[j]-ym[jj-1])
                                          / ((xm[ii-1]-xm[ii]) * (ym[jj]-ym[jj-1]))
                              + fieldm[jj][ii] * (x[i]-xm[ii-1]) * (y[j]-ym[jj-1])
	               	                  / ((xm[ii]-xm[ii-1]) * (ym[jj]-ym[jj-1]));
		}
	    }
	}
    }
}


double intlin(double x, double y1, double x1, double y2, double x2)
{
  /*
    xlin - lineare interpolation

    Uwe Schulzweida  04/05/1995
  */
  double value;
  
  value = (y2*(x-x1)+y1*(x2-x)) / (x2-x1);

  return (value);
}


void intlinarr(int nxm, double *ym, double *xm, int nx, double *y, double *x)
{
  /*
    xlinarr - lineare interpolation over 1D array

    Uwe Schulzweida  04/05/1995
  */
  int j, jj;

  for ( jj = 1; jj < nxm; jj++ )
    for ( j = 0; j < nx; j++ )
      if ( x[j] >= xm[jj-1] && x[j] <= xm[jj] )
	y[j] = intlin(x[j], ym[jj-1], xm[jj-1], ym[jj], xm[jj]);
}


void intgrid(FIELD *field1, FIELD *field2)
{
  static char func[] = "intgrid";
  int nlonIn, nlatIn;
  int nlonOut, nlatOut;
  int ilat, ilon;
  int gridIDin, gridIDout;
  int i, nmiss;
  double *lonIn, *latIn;
  double *lonOut, *latOut;
  double **fieldIn;
  double **field;
  double *array = NULL;
  double *arrayIn, *arrayOut;
  double missval;
  /* static int index = 0; */

  gridIDin  = field1->grid;
  gridIDout = field2->grid;
  arrayIn   = field1->ptr;
  arrayOut  = field2->ptr;
  missval   = field1->missval;

  if ( ! (gridInqXvals(gridIDin, NULL) && gridInqYvals(gridIDin, NULL)) )
    cdoAbort("Source grid has no values");

  nlonIn = gridInqXsize(gridIDin);
  nlatIn = gridInqYsize(gridIDin);
  lonIn = (double *) malloc(nlonIn*sizeof(double));
  latIn = (double *) malloc(nlatIn*sizeof(double));
  gridInqXvals(gridIDin, lonIn);
  gridInqYvals(gridIDin, latIn);

  if ( ! (gridInqXvals(gridIDout, NULL) && gridInqYvals(gridIDout, NULL)) )
    cdoAbort("Target grid has no values");

  nlonOut = gridInqXsize(gridIDout);
  nlatOut = gridInqYsize(gridIDout);
  lonOut = (double *) malloc(nlonOut*sizeof(double));
  latOut = (double *) malloc(nlatOut*sizeof(double));
  gridInqXvals(gridIDout, lonOut);
  gridInqYvals(gridIDout, latOut);

  fieldIn = (double **) malloc(nlatIn*sizeof(double *));

  for ( ilat = 0; ilat < nlatIn; ilat++ )
    fieldIn[ilat] = arrayIn + ilat*nlonIn;

  if ( nlonOut == 1 && nlatOut == 1 )
    {
      if ( lonOut[0] < lonIn[0] ) lonOut[0] += 360;

      if ( lonOut[0] > lonIn[nlonIn-1] )
	{
	  field = fieldIn;
	  fieldIn = (double **) malloc(nlatIn*sizeof(double *));
	  lonIn = (double *) realloc(lonIn, (nlonIn+1)*sizeof(double));
	  array = (double *) malloc(nlatIn*(nlonIn+1)*sizeof(double));

	  for ( ilat = 0; ilat < nlatIn; ilat++ )
	    {
	      fieldIn[ilat] = array + ilat*(nlonIn+1);  
	      memcpy(fieldIn[ilat], field[ilat], nlonIn*sizeof(double));
	      fieldIn[ilat][nlonIn] = fieldIn[ilat][0];
	      lonIn[nlonIn] = lonIn[0] + 360;
	    }
	  nlonIn++;
	  free(field);
	}

      if ( lonOut[0] < lonIn[0] || lonOut[0] > lonIn[nlonIn-1] )
	cdoAbort("Longitude %f out of bounds (%f to %f)!", lonOut[0], lonIn[0], lonIn[nlonIn-1]);

      if ( latOut[0] < MIN(latIn[0], latIn[nlatIn-1]) ||
	   latOut[0] > MAX(latIn[0], latIn[nlatIn-1]) )
	cdoAbort("Latitude %f out of bounds (%f to %f)!", latOut[0], latIn[0], latIn[nlatIn-1]);

      *arrayOut = intlinarr2p(nlonIn, nlatIn, fieldIn, lonIn, latIn, lonOut[0], latOut[0]);
      /*
      printf("%5d %f %f %f\n", index++, lonOut[0], latOut[0], *arrayOut);
      */
    }
  else
    {
      double **fieldOut;

      fieldOut = (double **) malloc(nlatOut * sizeof(double *));

      for ( ilat = 0; ilat < nlatOut; ilat++ )
	fieldOut[ilat] = arrayOut + ilat*nlonOut;

      for ( ilat = 0; ilat < nlatOut; ilat++ )
	for ( ilon = 0; ilon < nlonOut; ilon++ )
	  fieldOut[ilat][ilon] = 0;

      intlinarr2(missval,
		 nlonIn, nlatIn, fieldIn, lonIn, latIn,
		 nlonOut, nlatOut, fieldOut, lonOut, latOut);

      nmiss = 0;
      for ( i = 0; i < nlatOut*nlonOut; i++ )
	if ( DBL_IS_EQUAL(arrayOut[i], missval) ) nmiss++;

      field2->nmiss = nmiss;

      free(fieldOut);
    }

  if (array) free(array);
  free(lonIn);
  free(latIn);
  free(lonOut);
  free(latOut);
  free(fieldIn);
}


void intarea(FIELD *field1, FIELD *field2)
{
  static char func[] = "intarea";
  int gridsize_i, gridsize_o;
  int gridIDi;
  double *arrayIn;
  int gridIDo;
  double *arrayOut;
  double miss_val_8;
  int nmiss = 0;
  int jlat, jlon, i, j;
  int nlon, nlat, out_nlon, out_nlat;
  double xtm, xbm, xlm, xrm, xtd, xbd, xld, xrd;
  double dx, dy, xdx, xdxm, xt, xb, xr, xl;
  double pp;
  double **fieldm;
  double *lon, *lat, *lono, *lato;
  double **field, **plane;
  double *planedat;
  double *lon_array, *lat_array, *lono_array, *lato_array;
  int irun;
  int latfirst = 0, latlast = 0;
  int lonfirst = 0, lonlast;

  gridIDi    = field1->grid;
  gridIDo    = field2->grid;
  arrayIn    = field1->ptr;
  arrayOut   = field2->ptr;
  miss_val_8 = field1->missval;

  gridsize_i = gridInqSize(gridIDi);
  gridsize_o = gridInqSize(gridIDo);

  nlon  = gridInqXsize(gridIDi);
  nlat  = gridInqYsize(gridIDi);
  out_nlon = gridInqXsize(gridIDo);
  out_nlat = gridInqYsize(gridIDo);

  lon_array = (double *) malloc(nlon * sizeof(double));
  lat_array = (double *) malloc(nlat * sizeof(double));
  lon = lon_array;
  lat = lat_array;

  gridInqXvals(gridIDi, lon);
  gridInqYvals(gridIDi, lat);

  lono_array = (double *) malloc(out_nlon * sizeof(double));
  lono = lono_array;
  lato_array = (double *) malloc(out_nlat * sizeof(double));
  lato = lato_array;

  gridInqXvals(gridIDo, lono);
  gridInqYvals(gridIDo, lato);

  xdx  = 360./out_nlon;
  xdxm = 360./(nlon-2);
  /*
  printf("xdx, xdxm %.9g %.9g\n", xdx, xdxm);
  printf("nlon, nlat %d %d\n", nlon, nlat);
  printf("out_nlon, out_nlat %d %d\n", out_nlon, out_nlat);
  */

  fieldm = (double **) malloc(nlat*sizeof(double *));

  for ( j = 0; j < nlat; j++ )
    fieldm[j] = arrayIn + j*nlon;

  field = (double **) malloc(out_nlat*sizeof(double *));
  for ( j = 0; j < out_nlat; j++ )
    field[j] = arrayOut + j*out_nlon;

  planedat = (double *) malloc(out_nlat*out_nlon*sizeof(double));
  plane = (double **) malloc(out_nlat*sizeof(double *));
  for ( j = 0; j < out_nlat; j++ )
    plane[j] = planedat + j*out_nlon;

  for ( jlat = 0; jlat < out_nlat; jlat++ )
    for ( jlon = 0; jlon < out_nlon; jlon++ )
      {
        field[jlat][jlon] = 0;
        plane[jlat][jlon] = 0;
      }
  /*
  for ( i = 0; i < nlat; i++ ) printf("lat %5d %.9g\n", i+1, lat[i]);
  for ( i = 0; i < nlon; i++ ) printf("lon %5d %.9g\n", i+1, lon[i]);
  for ( i = 0; i < out_nlat; i++ ) printf("lato %5d %.9g\n", i+1, lato[i]);
  for ( i = 0; i < out_nlon; i++ ) printf("lono %5d %.9g\n", i+1, lono[i]);
  */
  
  /*
      WRITE(*,*) NLON,NLAT,OUT_NLON,OUT_NLAT
C      WRITE(*,*) XM
C      WRITE(*,*) YM
C      WRITE(*,*) X
C      WRITE(*,*) Y
      WRITE(*,*) XDX,XDXM
      WRITE(*,*)
  */

  for ( jlat = 0; jlat < out_nlat; jlat++ )
    {
      if ( jlat == 0 )
        xtm = 90.;
      else
        xtm = lato[jlat]+(lato[jlat-1]-lato[jlat])/2.;

      if ( jlat == out_nlat-1 )
        xbm = -90.;
      else
        xbm = lato[jlat]-(lato[jlat]-lato[jlat+1])/2.;
      /*
C      WRITE(*,*) JLAT,Y(JLAT),XBM,XTM
      */
      for ( j = 1; j < nlat; j++ )
	{
	  if ( lat[j-1] > xbm && lat[j] < xtm )
	    {
	      latfirst = j - 1;
	      break;
	    }
	}
      if ( j == nlat ) latfirst = nlat;

      for ( j = latfirst+1; j < nlat; j++ )
	{
	  if ( lat[j] < xbm )
	    {
	      latlast = j + 1;
	      break;
	    }
	}
      if ( j >= nlat ) latlast = nlat;
      
      for ( jlon = 0; jlon < out_nlon; jlon++ )
	{
	  irun = 0;
	  /*
C       IF (JLAT.EQ.45) WRITE(*,*) 'LON=',JLON
	  */
	  xlm = lono[jlon] - xdx/2;
	  xrm = lono[jlon] + xdx/2;
	/*
C        WRITE(*,*) X(JLON),XLM,XRM
	*/
	  for ( j = 1; j < nlon; j++ )
	    {
	      if ( lon[j-1] < xrm && lon[j] > xlm )
		{
		  lonfirst = j - 1;
		  break;
		}
	    }
	  if ( j == nlon ) lonfirst = nlon;

	  for ( j = lonfirst+1; j < nlon; j++ )
	    {
	      if ( lon[j] > xrm )
		{
		  lonlast = j + 1;
		  break;
		}
	    }
	  if ( j >= nlon ) lonlast = nlon;
	  /*  printf("%d %d %d %d\n", latfirst, latlast, lonfirst, lonlast); */

	  for ( j = latfirst; j < latlast; j++ )
	    {
	      if ( j == 0 )
		xtd = lat[0];
	      else
		xtd = lat[j] + (lat[j-1]-lat[j])/2.;

	      if ( j == nlat-1 )
		xbd = lat[nlat-1];
	      else
		xbd = lat[j] - (lat[j]-lat[j+1])/2.;

	      for ( i = 0; i < nlon; i++ )
		{
		  if ( !DBL_IS_EQUAL(fieldm[j][i], miss_val_8) )
		    {
		      xld = lon[i] - xdxm/2;
		      xrd = lon[i] + xdxm/2;
		/*
C        WRITE(*,*) XM(I),XLD,XRD
		*/
		      xt = MIN(xtd, xtm);
		      xb = MAX(xbd, xbm);
		      dy = MAX(0., xt-xb);

		      xr = MIN(xrd, xrm);
		      xl = MAX(xld, xlm);
		      dx = MAX(0., xr-xl);

		  /*
		  printf("xt, xb, dy, xr, xl, dx %8d %g %g %g %g %g %g\n", iii++, xt, xb, dy, xr, xl, dx);
		  */
		/*
C          IF (JLON.LT.5.AND.JLAT.LT.4.AND.I.LT.5.AND.J.LT.5) THEN
C          IF(DX.GT.0..AND.DY.GT.0.) THEN
C            WRITE(*,*)'X',JLON,I,XLM,XRM,XLD,XRD,DX,XL,XR
C            WRITE(*,*)'Y',JLAT,J,XBM,XTM,XBD,XTD,DY,XB,XT
C          ENDIF
C          ENDIF
C
		*/
		      field[jlat][jlon] = field[jlat][jlon] + fieldm[j][i]*dx*dy;
		      plane[jlat][jlon] = plane[jlat][jlon] + dx*dy;
		      irun++;
		    }
		}
	    }
	  /* printf("runs %d %d %d\n", jlat, jlon, irun); */
	}
    }

  for ( jlat = 0; jlat < out_nlat; jlat++ )
    for ( jlon = 0; jlon < out_nlon; jlon++ )
      {
        pp = plane[jlat][jlon];
        if ( pp > 0 )
          field[jlat][jlon] = field[jlat][jlon]/pp;
        else
	  {
	    field[jlat][jlon] = miss_val_8;
	    nmiss++;
	  }
      }

  field2->nmiss = nmiss;

  free(lon_array);
  free(lat_array);
  free(lono_array);
  free(lato_array);
  free(fieldm);
  free(field);
  free(planedat);
  free(plane);
}


/* source code from pingo */
void interpolate(FIELD *field1, FIELD *field2)
{
  static char func[] = "interpolate";
  int i;
  double *lono_array, *lato_array, *lono, *lato;
  double *lon_array, *lat_array, *lon, *lat;
  int gridsize_i, gridsize_o;
  int gridIDi;
  double *arrayIn;
  int gridIDo;
  double *arrayOut;
  double missval;
  int nmiss;
  long ilon, ilat, nxlon, nxlat, olon, olat;
  long l11, l12, l21, l22, l1, l2;
  double volon1, volon2, volat1, volat2;
  double *volon11, *volon12;
  double *volon21, *volon22;
  double vilon1, vilon2, vilat1, vilat2;
  double vlon1, vlon2, vlat1, vlat2;
  long *ilon11, *ilon12;
  long *ilon21, *ilon22;
  long *ilat1, *ilat2;
  long ilon1, ilon2;
  double sum, wsum;
  int k, n;
  double *xin_array, *xlon, *xlat;
  double **in0, **xin, **xout;
  int wrap_around, xlat_is_ascending;
  double a11, a12, a21, a22, b11, b12, b21, b22, t;
  double faclon1, faclon2, faclat1, faclat2;
  int nlon, nlat, out_nlon, out_nlat;

  gridIDi  = field1->grid;
  gridIDo  = field2->grid;
  arrayIn  = field1->ptr;
  arrayOut = field2->ptr;
  missval  = field1->missval;

  gridsize_i = gridInqSize(gridIDi);
  gridsize_o = gridInqSize(gridIDo);

  nlon  = gridInqXsize(gridIDi);
  nlat  = gridInqYsize(gridIDi);
  out_nlon = gridInqXsize(gridIDo);
  out_nlat = gridInqYsize(gridIDo);

  lon_array = (double *) malloc((nlon + 2) * sizeof(double));
  lat_array = (double *) malloc((nlat + 2) * sizeof(double));
  lon = lon_array + 1;
  lat = lat_array + 1;

  if ( ! (gridInqXvals(gridIDi, NULL) && gridInqYvals(gridIDi, NULL)) )
    cdoAbort("Source grid has no values");

  if ( ! (gridInqXvals(gridIDo, NULL) && gridInqYvals(gridIDo, NULL)) )
    cdoAbort("Target grid has no values");

  gridInqXvals(gridIDi, lon);
  gridInqYvals(gridIDi, lat);

  if ( nlon > 1 )
    {
      lon[-1] = lon[nlon - 1] - 360 > 2*lon[0] - lon[1] ?
	        lon[nlon - 1] - 360 : 2*lon[0] - lon[1];
      lon[nlon] = lon[0] + 360 < 2*lon[nlon-1] - lon[nlon-2] ?
 	          lon[0] + 360 : 2*lon[nlon-1] - lon[nlon-2];
    }
  else
    {
      lon[-1] = lon[0] - 360;
      lon[ 1] = lon[0] + 360;
    }

  if ( nlat > 1 )
    {
      lat[-1]   = 2*lat[0] - lat[1];
      lat[nlat] = 2*lat[nlat-1] - lat[nlat-2];
    }
  else
    {
      lat[-1] = lat[0] - 10;
      lat[ 1] = lat[nlat-1] + 10;
    }

  if ( lat[-1]   < -90 ) lat[-1] = -99;
  if ( lat[-1]   >  90 ) lat[-1] =  99;
  if ( lat[nlat] < -90 ) lat[nlat] = -99;
  if ( lat[nlat] >  90 ) lat[nlat] =  99;

  lono_array = (double *) malloc((out_nlon < 2 ? 4 : out_nlon + 2) * sizeof(double));
  lono = lono_array + 1;
  lato_array = (double *) malloc((out_nlat < 2 ? 4 : out_nlat + 2) * sizeof(double));
  lato = lato_array + 1;

  gridInqXvals(gridIDo, lono);
  gridInqYvals(gridIDo, lato);

  for ( i = 0; i < out_nlon - 1; i++ )
    if (lono[i + 1] <= lono[i]) break;

  for ( i++; i < out_nlon; i++ )
    {
      lono[i] += 360;
      if ( i < out_nlon - 1 && lono[i + 1] + 360 <= lono[i] )
	cdoAbort("Longitudes of output grid are not in ascending order!");
    }

  if ( lono[out_nlon - 1] - lono[0] >= 360 )
    cdoAbort("The area covered by the longitudes of output grid must not overlap!");

  if ( lato[0] >  90.001 || lato[out_nlat - 1] >  90.001 ||
       lato[0] < -90.001 || lato[out_nlat - 1] < -90.001 )
    {
      cdoAbort("Latitudes of output grid must be between 90 and -90!");
    }

  for ( i = 0; i < out_nlat - 1; i++ )
    if ( IS_EQUAL(lato[i + 1], lato[i]) || (i < out_nlat - 2 &&
	((lato[i + 1] > lato[i]) != (lato[i + 2] > lato[i + 1]))) )
      {
	cdoAbort("Latitudes of output grid must be in descending or ascending order!");
      }

  if ( out_nlon > 1 )
    {
      lono[-1] = lono[out_nlon - 1] - 360 > 2 * lono[0] - lono[1] ?
	            lono[out_nlon - 1] - 360 : 2 * lono[0] - lono[1];
      lono[out_nlon] = lono[0] + 360 < 2 * lono[out_nlon - 1] - lono[out_nlon - 2] ?
	                  lono[0] + 360 : 2 * lono[out_nlon - 1] - lono[out_nlon - 2];
    }
  else
    {
      lono[-1] = lono[0] - 360;
      lono[ 1] = lono[0] + 360;
    }

  if ( out_nlat > 1 )
    {
      lato[-1]   = 2*lato[0] - lato[1];
      lato[out_nlat] = 2*lato[out_nlat-1] - lato[out_nlat-2];
    }
  else
    {
      lato[-1] = lato[0] - 10;
      lato[ 1] = lato[out_nlat-1] + 10;
    }

  if ( lato[-1]   < -90 ) lato[-1] = -99;
  if ( lato[-1]   >  90 ) lato[-1] =  99;
  if ( lato[out_nlat] < -90 ) lato[out_nlat] = -99;
  if ( lato[out_nlat] >  90 ) lato[out_nlat] =  99;

  nxlon = 2*nlon + 1;
  nxlat = 2*nlat + 1;
  xin_array = (double *) malloc(nxlon * nxlat * sizeof(double));
  xin = (double **) malloc(nxlat * sizeof(double *));

  for (ilat = 0; ilat < nxlat; ilat++)
    xin[ilat] = xin_array + ilat * nxlon;

  xlon = (double *) malloc(nxlon * sizeof (double));
  for ( ilon = 0; ilon < nlon; ilon++ )
    {
      xlon[2*ilon + 1] = lon[ilon];
      xlon[2*ilon] = (lon[ilon - 1] + lon[ilon]) / 2;
    }
  xlon[2 * nlon] = (lon[nlon - 1] + lon[nlon]) / 2;

  xlat = (double *) malloc((2 * nlat + 1) * sizeof (double));
  for ( ilat = 0; ilat < nlat; ilat++ )
    {
      xlat[2*ilat + 1] = lat[ilat];
      xlat[2*ilat] = (lat[ilat - 1] + lat[ilat]) / 2;
    }
  xlat[2 * nlat] = (lat[nlat - 1] + lat[nlat]) / 2;

  in0 = (double **) malloc(nlat * sizeof (double *));
  for (ilat = 0; ilat < nlat; ilat++)
    in0[ilat] = arrayIn + ilat * nlon;

  ilon11 = (long *) malloc(out_nlon * sizeof(long));
  ilon12 = (long *) malloc(out_nlon * sizeof(long));
  ilon21 = (long *) malloc(out_nlon * sizeof(long));
  ilon22 = (long *) malloc(out_nlon * sizeof(long));
  volon11 = (double *) malloc(out_nlon * sizeof(double));
  volon12 = (double *) malloc(out_nlon * sizeof(double));
  volon21 = (double *) malloc(out_nlon * sizeof(double));
  volon22 = (double *) malloc(out_nlon * sizeof(double));

  for (olon = 0; olon < out_nlon; olon++)
    {
      volon1 = (lono[olon - 1] + lono[olon]) / 2;
      volon2 = (lono[olon] + lono[olon + 1]) / 2;
      if ( IS_EQUAL(volon1, volon2) ) volon2 += 360;
      volon2 -= 360 * floor((volon1 - xlon[0]) / 360);
      volon1 -= 360 * floor((volon1 - xlon[0]) / 360);
      volon21[olon] = volon1;
      volon22[olon] = volon2;
      for (l21 = 0; l21 < nxlon && xlon[l21] < volon1; l21++);
      for (l22 = l21; l22 < nxlon && xlon[l22] < volon2; l22++);
      volon1 -= 360;
      volon2 -= 360;
      volon11[olon] = volon1;
      volon12[olon] = volon2;
      for (l11 = 0; xlon[l11] < volon1; l11++);
      for (l12 = l11; l12 < nxlon && xlon[l12] < volon2; l12++);
      ilon11[olon] = l11;
      ilon12[olon] = l12;
      ilon21[olon] = l21;
      ilon22[olon] = l22;
    }

  ilat1 = (long *) malloc(out_nlat * sizeof(long));
  ilat2 = (long *) malloc(out_nlat * sizeof(long));

  xlat_is_ascending = xlat[0] <= xlat[nxlat - 1];
  for ( olat = 0; olat < out_nlat; olat++ )
    {
      volat1 = (lato[olat - 1] + lato[olat]) / 2;
      volat2 = (lato[olat] + lato[olat + 1]) / 2;
      if (!xlat_is_ascending)
	{
	  if (volat1 > volat2)
	    {
	      for (l1 =  0; l1 < nxlat && xlat[l1] > volat1; l1++);
	      for (l2 = l1; l2 < nxlat && xlat[l2] > volat2; l2++);
	    }
	  else
	    {
	      for (l1 =  0; l1 < nxlat && xlat[l1] > volat2; l1++);
	      for (l2 = l1; l2 < nxlat && xlat[l2] > volat1; l2++);
	    }
	}
      else
	{
	  if (volat1 < volat2)
	    {
	      for (l1 =  0; l1 < nxlat && xlat[l1] < volat1; l1++);
	      for (l2 = l1; l2 < nxlat && xlat[l2] < volat2; l2++);
	    }
	  else
	    {
	      for (l1 =  0; l1 < nxlat && xlat[l1] < volat2; l1++);
	      for (l2 = l1; l2 < nxlat && xlat[l2] < volat1; l2++);
	    }
	}

      ilat1[olat] = l1;
      ilat2[olat] = l2;
    }

  xout = (double **) malloc(out_nlat * sizeof (double *));
  for (olat = 0; olat < out_nlat; olat++)
    xout[olat] = arrayOut + olat * out_nlon;

  wrap_around = nlon > 1 && (lon[nlon - 1] >= lon[-1] + 360 - 0.001
			      || lon[nlon] >= lon[ 0] + 360 - 0.001);

  for (ilat = 0; ilat < nlat; ilat++)
#if defined (SX)
#pragma vdir nodep
#endif
    for (ilon = 0; ilon < nlon; ilon++)
      xin[2 * ilat + 1][2 * ilon + 1] = in0[ilat][ilon];

  for (ilat = 0; ilat < nxlat; ilat += 2)
#if defined (SX)
#pragma vdir nodep
#endif
    for (ilon = 1; ilon < nxlon; ilon += 2)
      {
	sum = 0;
	n = 0;
	if (ilat > 0 && !DBL_IS_EQUAL(xin[ilat - 1][ilon], missval))
	  {
	    sum += xin[ilat - 1][ilon];
	    n++;
	  }
	if (ilat < nxlat - 1 && !DBL_IS_EQUAL(xin[ilat + 1][ilon], missval))
	  {
	    sum += xin[ilat + 1][ilon];
	    n++;
	  }
	xin[ilat][ilon] = n ? sum / n : missval;
      }

  for ( ilat = 1; ilat < nxlat; ilat += 2 )
#if defined (SX)
#pragma vdir nodep
#endif
    for ( ilon = 0; ilon < nxlon; ilon += 2 )
      {
	sum = 0;
	n = 0;
	if (ilon > 0 && !DBL_IS_EQUAL(xin[ilat][ilon - 1], missval))
	  {
	    sum += xin[ilat][ilon - 1];
	    n++;
	  }
	if (ilon == 0 && wrap_around && !DBL_IS_EQUAL(xin[ilat][2 * nlon - 1], missval))
	  {
	    sum += xin[ilat][2 * nlon - 1];
	    n++;
	  }
	if (ilon < nxlon - 1 && !DBL_IS_EQUAL(xin[ilat][ilon + 1], missval))
	  {
	    sum += xin[ilat][ilon + 1];
	    n++;
	  }
	if (ilon == nxlon - 1 && wrap_around && !DBL_IS_EQUAL(xin[ilat][1], missval))
	  {
	    sum += xin[ilat][1];
	    n++;
	  }
	xin[ilat][ilon] = n ? sum / n : missval;
      }

  for ( ilat = 0; ilat < nxlat; ilat += 2 )
#if defined (SX)
#pragma vdir nodep
#endif
    for ( ilon = 0; ilon < nxlon; ilon += 2 )
      {
	sum = 0;
	n = 0;
	if (ilon > 0 && !DBL_IS_EQUAL(xin[ilat][ilon - 1], missval))
	  {
	    sum += xin[ilat][ilon - 1];
	    n++;
	  }
	if (ilon == 0 && wrap_around && !DBL_IS_EQUAL(xin[ilat][2 * nlon - 1], missval))
	  {
	    sum += xin[ilat][2 * nlon - 1];
	    n++;
	  }
	if (ilon < nxlon - 1 && !DBL_IS_EQUAL(xin[ilat][ilon + 1], missval))
	  {
	    sum += xin[ilat][ilon + 1];
	    n++;
	  }
	if (ilon == nxlon - 1 && wrap_around && !DBL_IS_EQUAL(xin[ilat][1], missval))
	  {
	    sum += xin[ilat][1];
	    n++;
	  }
	if (ilat > 0 && !DBL_IS_EQUAL(xin[ilat - 1][ilon], missval))
	  {
	    sum += xin[ilat - 1][ilon];
	    n++;
	  }
	if (ilat < nxlat - 1 && !DBL_IS_EQUAL(xin[ilat + 1][ilon], missval))
	  {
	    sum += xin[ilat + 1][ilon];
	    n++;
	  }
	xin[ilat][ilon] = n ? sum / n : missval;
      }

  for ( olat = 0; olat < out_nlat; olat++ )
    {
      if ( lato[-1] < lato[out_nlat] )
	{
	  volat1 = (lato[olat - 1] + lato[olat]) / 2;
	  volat2 = (lato[olat] + lato[olat + 1]) / 2;
	}
      else
	{
	  volat2 = (lato[olat - 1] + lato[olat]) / 2;
	  volat1 = (lato[olat] + lato[olat + 1]) / 2;
	}

      for ( olon = 0; olon < out_nlon; olon++ )
	{
	  sum = 0;
	  wsum = 0;
	  for (k = 0; k < 2; k++)
	    {
	      if (k == 0)
		{
		  ilon1 = ilon11[olon];
		  ilon2 = ilon12[olon];
		  volon1 = volon11[olon];
		  volon2 = volon12[olon];
		}
	      else
		{
		  ilon1 = ilon21[olon];
		  ilon2 = ilon22[olon];
		  volon1 = volon21[olon];
		  volon2 = volon22[olon];
		}

	      for ( ilon = ilon1; ilon <= ilon2; ilon++ )
		{
		  if ( ilon == 0 || ilon == nxlon ) continue;
		  vilon1 = xlon[ilon - 1];
		  vilon2 = xlon[ilon];
		  for ( ilat = ilat1[olat]; ilat <= ilat2[olat]; ilat++ )
		    {
		      if ( ilat == 0 || ilat == nxlat ) continue;
		      if ( xlat_is_ascending )
			{
			  vilat1 = xlat[ilat - 1];
			  vilat2 = xlat[ilat];
			  a11 = xin[ilat - 1][ilon - 1];
			  a12 = xin[ilat - 1][ilon];
			  a21 = xin[ilat][ilon - 1];
			  a22 = xin[ilat][ilon];
			}
		      else
			{
			  vilat1 = xlat[ilat];
			  vilat2 = xlat[ilat - 1];
			  a11 = xin[ilat][ilon - 1];
			  a12 = xin[ilat][ilon];
			  a21 = xin[ilat - 1][ilon - 1];
			  a22 = xin[ilat - 1][ilon];
			}
		      if ( DBL_IS_EQUAL(a11, missval) || DBL_IS_EQUAL(a12, missval) ||
			   DBL_IS_EQUAL(a21, missval) || DBL_IS_EQUAL(a22, missval) )
			{
			  continue;
			}
		      if ( volon1 <= vilon1 && vilon2 <= volon2 &&
			   volat1 <= vilat1 && vilat2 <= volat2 )
			{
			  vlon1 = vilon1 * M_PI / 180;
			  vlon2 = vilon2 * M_PI / 180;
			  vlat1 = vilat1 * M_PI / 180;
			  vlat2 = vilat2 * M_PI / 180;
			  b11 = a11;
			  b12 = a12;
			  b21 = a21;
			  b22 = a22;
			}
		      else
			{
			  vlon1 = (volon1 <= vilon1 ? vilon1 : volon1);
			  vlon2 = (vilon2 <= volon2 ? vilon2 : volon2);
			  vlat1 = (volat1 <= vilat1 ? vilat1 : volat1);
			  vlat2 = (vilat2 <= volat2 ? vilat2 : volat2);
			  if ( vlon1 >= vlon2 - (volon2 - volon1) * 1e-5 ||
			       vlat1 >= vlat2 - (volat2 - volat1) * 1e-5)
			    {
			      continue;
			    }
			  faclon1 = (vlon1 - vilon1) / (vilon2 - vilon1);
			  faclon2 = (vlon2 - vilon1) / (vilon2 - vilon1);
			  faclat1 = (vlat1 - vilat1) / (vilat2 - vilat1);
			  faclat2 = (vlat2 - vilat1) / (vilat2 - vilat1);
			  vlon1 *= M_PI / 180;
			  vlon2 *= M_PI / 180;
			  vlat1 *= M_PI / 180;
			  vlat2 *= M_PI / 180;
			  b11 = a11 + (a12 - a11)*faclon1 + (a21 - a11)*faclat1
			      + (a22 - a12 - a21 + a11)*faclon1*faclat1;
			  b12 = a11 + (a12 - a11)*faclon2 + (a21 - a11)*faclat1
			      + (a22 - a12 - a21 + a11)*faclon2*faclat1;
			  b21 = a11 + (a12 - a11)*faclon1 + (a21 - a11)*faclat2
			      + (a22 - a12 - a21 + a11)*faclon1*faclat2;
			  b22 = a11 + (a12 - a11)*faclon2 + (a21 - a11)*faclat2
			      + (a22 - a12 - a21 + a11)*faclon2*faclat2;
			}
		      wsum += (vlon2 - vlon1) * (sin(vlat2) - sin(vlat1));
		      t = 2 * sin((vlat2 + vlat1) / 2) *
		  	      sin((vlat2 - vlat1) / 2) / (vlat2 - vlat1);
		      sum += (vlon2 - vlon1) / 2 * ((b11 + b12) * (t - sin(vlat1)) +
				                    (b21 + b22) * (sin(vlat2) - t));
		    }
		}
	    }
	  xout[olat][olon] = IS_NOT_EQUAL(wsum, 0) ? sum / wsum : missval;
	}
    }

  nmiss = 0;
  for ( i = 0; i < gridsize_o; i++ )
    if ( DBL_IS_EQUAL(arrayOut[i], missval) ) nmiss++;

  field2->nmiss = nmiss;

  free(lon_array);
  free(lat_array);
  free(lono_array);
  free(lato_array);
  free(xin);
  free(xin_array);
  free(xlon);
  free(xlat);
  free(in0);
  free(ilon11);
  free(ilon12);
  free(ilon21);
  free(ilon22);
  free(volon11);
  free(volon12);
  free(volon21);
  free(volon22);
  free(ilat1);
  free(ilat2);
  free(xout);
}


/* source code from pingo */
void contrast(void)
{
  static char func[] = "contrast";

  int rec = 1;
  int nlat, nlon;
  int i, j, size = 0;
  double missval;
  double **work;
  double **in;
  double **out;
  double *lon;
  double *lat;

  static double *xin_array, **xin, **xout, **xwork[17];
  static double **r, **r_bar, **r_new, **r_bar_new;
  static double **p, **p_bar, **p_new, **p_bar_new, **swap;
  static double a0, a1, a2, a3, a4, a5, a6, a7, a8;
  static double dlon0, dlon1, dslat0, dslat1;
  static double w00, w01, w10, w11, wsum;
  static double lon0, lon1, lon2, lat0, lat1, lat2;
  static double flat00, flat01, flat10, flat11;
  static long ilon, ilat;
  static int determine_matrix, wrap_around, stop_iteration;
  static double table[2][2][2][4]
     = { {{{8 / 32., 8 / 32., 8 / 32., 8 / 32.},
	  {12 / 32., 8 / 32., 0, 12 / 32.}},
	 {{12 / 32., 0, 8 / 32., 12 / 32.},
	  {16 / 32., 0, 0, 16 / 32.}}},
	 {{{0, 12 / 32., 12 / 32., 8 / 32.},
	   {0, 16 / 32., 0, 16 / 32.}},
	  {{0, 0, 16 / 32., 16 / 32.},
	   {0, 0, 0, 32 / 32.}}}
     };
  static double *case_table;
  static double f;
  static double nom, denom, a, b;
  static const double eps = 1e-20;
  static double max;
  static int iter;
  static float iter_sum = 0, iter_n = 0;

  nlat = 0;
  nlon = 0;
  missval = 0;
  work = 0;
  in = 0;
  out = 0;
  lon = 0;
  lat = 0;

  if (rec == 1)
    {
      xin_array = (double *) malloc((nlat + 2) * (nlon + 2) * sizeof(double));
      xin = (double **) malloc((nlat + 2) * sizeof(double *));
      *xin = *(xin + 1);
      for (ilat = -1; ilat <= nlat; ilat++)
	xin[ilat] = xin_array + (ilat + 1) * (nlon + 2) + 1;
      xout = (double **) malloc(nlat * sizeof(double *));
      for (ilat = 0; ilat < nlat; ilat++)
	xout[ilat] = out[0] + ilat * nlon;
      for (j = 0; j < 17; j++)
	{
	  xwork[j] = (double **) malloc(nlat * sizeof(double *));
	  for (ilat = 0; ilat < nlat; ilat++)
	    xwork[j][ilat] = work[j] + ilat * nlon;
	}
      wrap_around = nlon > 1 && (lon[nlon - 1] >= lon[-1] + 360 - 0.001
				 || lon[nlon] >= lon[0] + 360 - 0.001);
      determine_matrix = TRUE;
    }
  else
    {
      determine_matrix = FALSE;
      for (i = 0; i < size; i++)
	if ( DBL_IS_EQUAL(in[0][i], missval) != DBL_IS_EQUAL(out[0][i], missval) )
	  {
	    determine_matrix = TRUE;
	    break;
	  }
    }
  for (ilon = 0; ilon < nlon; ilon++)
    {
      for (ilat = 0; ilat < nlat; ilat++)
	xin[ilat][ilon] = in[0][ilon + ilat * nlon];
      xin[-1][ilon] = xin[nlat][ilon] = missval;
    }
  if (wrap_around)
    for (ilat = -1; ilat <= nlat; ilat++)
      {
	xin[ilat][-1] = xin[ilat][nlon - 1];
	xin[ilat][nlon] = xin[ilat][0];
      }
  else
    for (ilat = -1; ilat <= nlat; ilat++)
      xin[ilat][-1] = xin[ilat][nlon] = missval;

  if (determine_matrix)
    {
      for (ilon = 0; ilon < nlon; ilon++)
	{
	  lon1 = lon[ilon];
	  lon0 = (lon[ilon - 1] + lon1) / 2;
	  lon2 = (lon[ilon + 1] + lon1) / 2;
	  dlon0 = lon1 - lon0;
	  dlon1 = lon2 - lon1;
	  for (ilat = 0; ilat < nlat; ilat++)
	    {
	      lat1 = lat[ilat] * M_PI / 180;
	      lat0 = (lat[ilat - 1] * M_PI / 180 + lat1) / 2;
	      lat2 = (lat[ilat + 1] * M_PI / 180 + lat1) / 2;
	      dslat0 = sin (lat1) - sin (lat0);
	      dslat1 = sin (lat2) - sin (lat1);
	      flat00 = 2 / (lat1 - lat0) * sin ((lat1 + lat0) / 2) * 
		       sin ((lat1 - lat0) /2) - sin (lat0);
	      flat01 = dslat0 - flat00;
	      flat10 = 2 / (lat2 - lat1) * sin ((lat2 + lat1) / 2) * 
                       sin ((lat2 - lat1) / 2) - sin (lat1);
	      flat11 = dslat1 - flat10;
	      flat00 /= 2;
	      flat01 /= 2;
	      flat10 /= 2;
	      flat11 /= 2;
	      if ( DBL_IS_EQUAL(xin[ilat][ilon], missval) )
		{
		  xwork[4][ilat][ilon] = 1;
		  xwork[0][ilat][ilon] = xwork[1][ilat][ilon]
		    = xwork[2][ilat][ilon] = xwork[3][ilat][ilon]
		    = xwork[5][ilat][ilon] = xwork[6][ilat][ilon]
		    = xwork[7][ilat][ilon] = xwork[8][ilat][ilon] = 0;
		}
	      else
		{
		  w00 = dslat0 * dlon0;
		  w01 = dslat0 * dlon1;
		  w10 = dslat1 * dlon0;
		  w11 = dslat1 * dlon1;
		  wsum = w00 + w01 + w10 + w11;

		  a4 = dlon0 * flat01;
		  if ( DBL_IS_EQUAL(xin[ilat - 1][ilon], missval) )
		    {
		      a1 = 0;
		      a4 += dlon0 * flat00;
		    }
		  else
		    {
		      a1 = dlon0 * flat00 / 2;
		      a4 += a1;
		    }
		  if ( DBL_IS_EQUAL(xin[ilat][ilon - 1], missval) )
		    {
		      a3 = 0;
		      a4 += dlon0 * flat01;
		    }
		  else
		    {
		      a3 = dlon0 * flat01 / 2;
		      a4 += a3;
		    }
		  case_table = table[DBL_IS_EQUAL(xin[ilat - 1][ilon - 1], missval)]
		                    [DBL_IS_EQUAL(xin[ilat - 1][ilon], missval)]
		                    [DBL_IS_EQUAL(xin[ilat][ilon - 1], missval)];
		  f = dlon0 * flat00;
		  a0 = case_table[0] * f;
		  a1 += case_table[1] * f;
		  a3 += case_table[2] * f;
		  a4 += case_table[3] * f;
		  xwork[0][ilat][ilon] = a0;
		  xwork[1][ilat][ilon] = a1;
		  xwork[3][ilat][ilon] = a3;
		  xwork[4][ilat][ilon] = a4;

		  a4 = dlon1 * flat01;
		  if ( DBL_IS_EQUAL(xin[ilat - 1][ilon], missval) )
		    {
		      a1 = 0;
		      a4 += dlon1 * flat00;
		    }
		  else
		    {
		      a1 = dlon1 * flat00 / 2;
		      a4 += a1;
		    }

		  if ( DBL_IS_EQUAL(xin[ilat][ilon + 1], missval) )
		    {
		      a5 = 0;
		      a4 += dlon1 * flat01;
		    }
		  else
		    {
		      a5 = dlon1 * flat01 / 2;
		      a4 += a5;
		    }
		  case_table = table[DBL_IS_EQUAL(xin[ilat - 1][ilon + 1], missval)]
		                    [DBL_IS_EQUAL(xin[ilat - 1][ilon], missval)]
		                    [DBL_IS_EQUAL(xin[ilat][ilon + 1], missval)];
		  f = dlon1 * flat00;
		  a2 = case_table[0] * f;
		  a1 += case_table[1] * f;
		  a5 += case_table[2] * f;
		  a4 += case_table[3] * f;
		  xwork[2][ilat][ilon] += a2;
		  xwork[1][ilat][ilon] += a1;
		  xwork[5][ilat][ilon] += a5;
		  xwork[4][ilat][ilon] += a4;

		  a4 = dlon0 * flat10;
		  if ( DBL_IS_EQUAL(xin[ilat + 1][ilon], missval) )
		    {
		      a7 = 0;
		      a4 += dlon0 * flat11;
		    }
		  else
		    {
		      a7 = dlon0 * flat11 / 2;
		      a4 += a7;
		    }

		  if ( DBL_IS_EQUAL(xin[ilat][ilon - 1], missval))
		    {
		      a3 = 0;
		      a4 += dlon0 * flat10;
		    }
		  else
		    {
		      a3 = dlon0 * flat10 / 2;
		      a4 += a3;
		    }

		  case_table = table[DBL_IS_EQUAL(xin[ilat + 1][ilon - 1], missval)]
		                    [DBL_IS_EQUAL(xin[ilat + 1][ilon], missval)]
		                    [DBL_IS_EQUAL(xin[ilat][ilon - 1], missval)];
		  f = dlon0 * flat11;
		  a6 = case_table[0] * f;
		  a7 += case_table[1] * f;
		  a3 += case_table[2] * f;
		  a4 += case_table[3] * f;
		  xwork[6][ilat][ilon] += a6;
		  xwork[7][ilat][ilon] += a7;
		  xwork[3][ilat][ilon] += a3;
		  xwork[4][ilat][ilon] += a4;

		  a4 = dlon1 * flat10;
		  if ( DBL_IS_EQUAL(xin[ilat + 1][ilon], missval) )
		    {
		      a7 = 0;
		      a4 += dlon1 * flat11;
		    }
		  else
		    {
		      a7 = dlon1 * flat11 / 2;
		      a4 += a7;
		    }

		  if ( DBL_IS_EQUAL(xin[ilat][ilon + 1], missval) )
		    {
		      a5 = 0;
		      a4 += dlon1 * flat10;
		    }
		  else
		    {
		      a5 = dlon1 * flat10 / 2;
		      a4 += a5;
		    }

		  case_table = table[DBL_IS_EQUAL(xin[ilat + 1][ilon + 1], missval)]
		                    [DBL_IS_EQUAL(xin[ilat + 1][ilon], missval)]
		                    [DBL_IS_EQUAL(xin[ilat][ilon + 1], missval)];
		  f = dlon1 * flat11;
		  a8 = case_table[0] * f;
		  a7 += case_table[1] * f;
		  a5 += case_table[2] * f;
		  a4 += case_table[3] * f;
		  xwork[8][ilat][ilon] += a8;
		  xwork[7][ilat][ilon] += a7;
		  xwork[5][ilat][ilon] += a5;
		  xwork[4][ilat][ilon] += a4;

		  xwork[0][ilat][ilon] /= wsum;
		  xwork[1][ilat][ilon] /= wsum;
		  xwork[2][ilat][ilon] /= wsum;
		  xwork[3][ilat][ilon] /= wsum;
		  xwork[4][ilat][ilon] /= wsum;
		  xwork[5][ilat][ilon] /= wsum;
		  xwork[6][ilat][ilon] /= wsum;
		  xwork[7][ilat][ilon] /= wsum;
		  xwork[8][ilat][ilon] /= wsum;
		}
	    }
	}
    }

  /* Solve sparse linear equation system
     xin[ilat][ilon] = xwork[0][ilat][ilon]*xout[ilat-1][ilon-1] 
                     + xwork[1][ilat][ilon]*xout[ilat-1][ilon  ] 
	             + xwork[2][ilat][ilon]*xout[ilat-1][ilon+1] 
	             + xwork[3][ilat][ilon]*xout[ilat  ][ilon-1] 
	             + xwork[4][ilat][ilon]*xout[ilat  ][ilon  ] 
	             + xwork[5][ilat][ilon]*xout[ilat  ][ilon+1] 
	             + xwork[6][ilat][ilon]*xout[ilat+1][ilon-1] 
	             + xwork[7][ilat][ilon]*xout[ilat+1][ilon  ] 
	             + xwork[8][ilat][ilon]*xout[ilat+1][ilon+1]
     using the biconjugate gradient method */

  max = 0;
  for (ilat = 0; ilat < nlat; ilat++)
    for (ilon = 0; ilon < nlon; ilon++)
      {
	f = xin[ilat][ilon];
	if (!DBL_IS_EQUAL(f, missval) && f > max)
	  max = f;
      }

  r = xwork[9];
  r_bar = xwork[10];
  r_new = xwork[11];
  r_bar_new = xwork[12];
  p = xwork[13];
  p_bar = xwork[14];
  p_new = xwork[15];
  p_bar_new = xwork[16];

  for (ilat = 0; ilat < nlat; ilat++)
    for (ilon = 0; ilon < nlon; ilon++)
      xout[ilat][ilon] = xin[ilat][ilon];

  for (ilat = 0; ilat < nlat; ilat++)
    for (ilon = 0; ilon < nlon; ilon++)
      r[ilat][ilon] = r_bar[ilat][ilon]
	            = p[ilat][ilon]
	            = p_bar[ilat][ilon]
	            = xin[ilat][ilon] - (xwork[0][ilat][ilon] * xin[ilat - 1][ilon - 1]
	                              + xwork[1][ilat][ilon] * xin[ilat - 1][ilon]
	                              + xwork[2][ilat][ilon] * xin[ilat - 1][ilon + 1]
	                              + xwork[3][ilat][ilon] * xin[ilat][ilon - 1]
	                              + xwork[4][ilat][ilon] * xin[ilat][ilon]
	                              + xwork[5][ilat][ilon] * xin[ilat][ilon + 1]
	                              + xwork[6][ilat][ilon] * xin[ilat + 1][ilon - 1]
	                              + xwork[7][ilat][ilon] * xin[ilat + 1][ilon]
	                              + xwork[8][ilat][ilon] * xin[ilat + 1][ilon + 1]);

  for (iter = 1;; iter++)
    {
      stop_iteration = TRUE;
      for (ilat = 0; ilat < nlat; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  if (fabs (r[ilat][ilon]) > eps * max)
	    {
	      stop_iteration = FALSE;
	      break;
	    }
      if (stop_iteration)
	break;
      /*
      if (user_asked)
	{
	  lock ();
	  fprintf (stderr,
		   "%s: Status: Raising contrast of record %d"
		   " iteration step %d", prompt, rec, iter);
	  if (iter_n)
	    fprintf (stderr, " of approximately %d.\n",
		     (int) (iter_sum / iter_n));
	  else
	    fputs (".\n", stderr);
	  fflush (stderr);
	  unlock ();
	  user_asked = FALSE;
	}
      */
      if (iter == 1)
	{
	  nom = 0;
	  for (ilat = 0; ilat < nlat; ilat++)
	    for (ilon = 0; ilon < nlon; ilon++)
	      nom += r[ilat][ilon] * r_bar[ilat][ilon];
	}
      for (ilat = 0; ilat < nlat; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  r_new[ilat][ilon] = xwork[4][ilat][ilon] * p[ilat][ilon];
      for (ilat = 1; ilat < nlat; ilat++)
	{
	  for (ilon = 1; ilon < nlon; ilon++)
	    r_new[ilat][ilon] +=
	      xwork[0][ilat][ilon] * p[ilat - 1][ilon - 1];
	  r_new[ilat][0] += xwork[0][ilat][0] * p[ilat - 1][nlon - 1];
	}
      for (ilat = 1; ilat < nlat; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  r_new[ilat][ilon] += xwork[1][ilat][ilon] * p[ilat - 1][ilon];
      for (ilat = 1; ilat < nlat; ilat++)
	{
	  for (ilon = 0; ilon < nlon - 1; ilon++)
	    r_new[ilat][ilon] +=
	      xwork[2][ilat][ilon] * p[ilat - 1][ilon + 1];
	  r_new[ilat][nlon - 1] +=
	    xwork[2][ilat][nlon - 1] * p[ilat - 1][0];
	}
      for (ilat = 0; ilat < nlat; ilat++)
	{
	  for (ilon = 1; ilon < nlon; ilon++)
	    r_new[ilat][ilon] += xwork[3][ilat][ilon] * p[ilat][ilon - 1];
	  r_new[ilat][0] += xwork[3][ilat][0] * p[ilat][nlon - 1];
	}
      for (ilat = 0; ilat < nlat; ilat++)
	{
	  for (ilon = 0; ilon < nlon - 1; ilon++)
	    r_new[ilat][ilon] +=
	      xwork[5][ilat][ilon] * p[ilat][ilon + 1];
	  r_new[ilat][nlon - 1] +=
	    xwork[5][ilat][nlon - 1] * p[ilat][0];
	}
      for (ilat = 0; ilat < nlat - 1; ilat++)
	{
	  for (ilon = 1; ilon < nlon; ilon++)
	    r_new[ilat][ilon] +=
	      xwork[6][ilat][ilon] * p[ilat + 1][ilon - 1];
	  r_new[ilat][0] += xwork[6][ilat][0] * p[ilat + 1][nlon - 1];
	}
      for (ilat = 0; ilat < nlat - 1; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  r_new[ilat][ilon] += xwork[7][ilat][ilon] * p[ilat + 1][ilon];
      for (ilat = 0; ilat < nlat - 1; ilat++)
	{
	  for (ilon = 0; ilon < nlon - 1; ilon++)
	    r_new[ilat][ilon] +=
	      xwork[8][ilat][ilon] * p[ilat + 1][ilon + 1];
	  r_new[ilat][nlon - 1] +=
	    xwork[8][ilat][nlon - 1] * p[ilat + 1][0];
	}
      denom = 0;
      for (ilat = 0; ilat < nlat; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  denom += p_bar[ilat][ilon] * r_new[ilat][ilon];

      if ( IS_EQUAL(denom, 0) ) break;

      a = nom / denom;

      for (ilat = 0; ilat < nlat; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  r_new[ilat][ilon] = r[ilat][ilon] - a * r_new[ilat][ilon];
      for (ilat = 0; ilat < nlat; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  r_bar_new[ilat][ilon] =
	    xwork[4][ilat][ilon] * p_bar[ilat][ilon];
      for (ilat = 1; ilat < nlat; ilat++)
	{
	  for (ilon = 1; ilon < nlon; ilon++)
	    r_bar_new[ilat][ilon] += xwork[8][ilat - 1][ilon - 1] * p_bar[ilat - 1][ilon - 1];

	  r_bar_new[ilat][0] += xwork[8][ilat - 1][nlon - 1] * p_bar[ilat - 1][nlon - 1];
	}

      for (ilat = 1; ilat < nlat; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  r_bar_new[ilat][ilon] += xwork[7][ilat - 1][ilon] * p_bar[ilat - 1][ilon];

      for (ilat = 1; ilat < nlat; ilat++)
	{
	  for (ilon = 0; ilon < nlon - 1; ilon++)
	    r_bar_new[ilat][ilon] += xwork[6][ilat - 1][ilon + 1] * p_bar[ilat - 1][ilon + 1];

	  r_bar_new[ilat][nlon - 1] += xwork[6][ilat - 1][0] * p_bar[ilat - 1][0];
	}
      for (ilat = 0; ilat < nlat; ilat++)
	{
	  for (ilon = 1; ilon < nlon; ilon++)
	    r_bar_new[ilat][ilon] += xwork[5][ilat][ilon - 1] * p_bar[ilat][ilon - 1];
	  r_bar_new[ilat][0] += xwork[5][ilat][nlon - 1] * p_bar[ilat][nlon - 1];
	}

      for (ilat = 0; ilat < nlat; ilat++)
	{
	  for (ilon = 0; ilon < nlon - 1; ilon++)
	    r_bar_new[ilat][ilon] += xwork[3][ilat][ilon + 1] * p_bar[ilat][ilon + 1];

	  r_bar_new[ilat][nlon - 1] += xwork[3][ilat][0] * p_bar[ilat][0];
	}
      for (ilat = 0; ilat < nlat - 1; ilat++)
	{
	  for (ilon = 1; ilon < nlon; ilon++)
	    r_bar_new[ilat][ilon] += xwork[2][ilat + 1][ilon - 1] * p_bar[ilat + 1][ilon - 1];
	  r_bar_new[ilat][0] += xwork[2][ilat + 1][nlon - 1] * p_bar[ilat + 1][nlon - 1];
	}
      for (ilat = 0; ilat < nlat - 1; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  r_bar_new[ilat][ilon] += xwork[1][ilat + 1][ilon] * p_bar[ilat + 1][ilon];
      for (ilat = 0; ilat < nlat - 1; ilat++)
	{
	  for (ilon = 0; ilon < nlon - 1; ilon++)
	    r_bar_new[ilat][ilon] += xwork[0][ilat + 1][ilon + 1] * p_bar[ilat + 1][ilon + 1];
	  r_bar_new[ilat][nlon - 1] += xwork[0][ilat + 1][0] * p_bar[ilat + 1][0];
	}
      for (ilat = 0; ilat < nlat; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  r_bar_new[ilat][ilon] = r_bar[ilat][ilon] - a * r_bar_new[ilat][ilon];
      denom = nom;
      nom = 0;
      for (ilat = 0; ilat < nlat; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  nom += r_bar_new[ilat][ilon] * r_new[ilat][ilon];
      if ( IS_EQUAL(denom, 0) )
	break;
      b = nom / denom;
      for (ilat = 0; ilat < nlat; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  {
	    p_new[ilat][ilon] = r_new[ilat][ilon] + b * p[ilat][ilon];
	    p_bar_new[ilat][ilon] = r_bar_new[ilat][ilon] + b * p_bar[ilat][ilon];
	  }
      for (ilat = 0; ilat < nlat; ilat++)
	for (ilon = 0; ilon < nlon; ilon++)
	  if ( !DBL_IS_EQUAL(xout[ilat][ilon], missval) )
	    xout[ilat][ilon] += a * p[ilat][ilon];
      swap = r_new;
      r_new = r;
      r = swap;
      swap = r_bar_new;
      r_bar_new = r_bar;
      r_bar = swap;
      swap = p_new;
      p_new = p;
      p = swap;
      swap = p_bar_new;
      p_bar_new = p_bar;
      p_bar = swap;
    }
  iter_sum = iter_sum * 0.9 + iter;
  iter_n = iter_n * 0.9 + 1;

  free(xin_array);
  free(xin);
  free(xout);
  for (j = 0; j < 17; j++)
    free(xwork[j]);
}
