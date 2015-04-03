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

#include "cdo.h"
#include "cdo_int.h"
#include "cdi.h"
/* RQ */
#include "nth_element.h"
/* QR */

double _FADD_(double x, double y, double missval1, double missval2)
{
  return FADD(x,y);
}

double _FSUB_(double x, double y, double missval1, double missval2)
{
  return FSUB(x, y);
}

double _FMUL_(double x, double y, double missval1, double missval2)
{
  return FMUL(x, y);
}

double _FDIV_(double x, double y, double missval1, double missval2)
{
  return FDIV(x, y);
}

double _FROOT_(double x, double missval1)
{
  return FROOT(x);
}


double fldfun(FIELD field, int function)
{
  double rval = 0;

  if      ( function == func_min  )  rval = fldmin(field);
  else if ( function == func_max  )  rval = fldmax(field);  
  else if ( function == func_sum  )  rval = fldsum(field);  
  else if ( function == func_mean )  rval = fldmean(field);  
  else if ( function == func_avg  )  rval = fldavg(field);  
  else if ( function == func_std  )  rval = fldstd(field);  
  else if ( function == func_var  )  rval = fldvar(field);
  else cdoAbort("function %d not implemented!", function);

  return rval;
}

double fldmin(FIELD field)
{
  int i;
  int    len     = field.size;
  int    nmiss   = field.nmiss;
  double missval = field.missval;
  double *array  = field.ptr;
  double rmin = 0;

  if ( nmiss > 0 )
    {
      rmin = DBL_MAX;
      for ( i = 0; i < len; i++ ) 
	if ( !DBL_IS_EQUAL(array[i], missval) )
	  if ( array[i] < rmin ) rmin = array[i];

      if ( IS_EQUAL(rmin, DBL_MAX) )
	rmin = missval;
    }
  else
    {
      rmin = array[0];
      for ( i = 1; i < len; i++ ) 
	if ( array[i] < rmin )  rmin = array[i];
    }

  return (rmin);
}


double fldmax(FIELD field)
{
  int i;
  int    len     = field.size;
  int    nmiss   = field.nmiss;
  double missval = field.missval;
  double *array  = field.ptr;
  double rmax = 0;

  if ( nmiss > 0 )
    {
      rmax = -DBL_MAX;
      for ( i = 0; i < len; i++ )
        if ( !DBL_IS_EQUAL(array[i], missval) )
          if ( array[i] > rmax ) rmax = array[i];
      
      if ( IS_EQUAL(rmax, -DBL_MAX) )
        rmax = missval;
    }
  else
    {
      rmax = array[0];
      for ( i = 1; i < len; i++ ) 
        if ( array[i] > rmax )  rmax = array[i];
    }

  return (rmax);
}


double fldsum(FIELD field)
{
  int i;
  int    len     = field.size;
  int    nmiss   = field.nmiss;
  double missval = field.missval;
  double *array  = field.ptr;
  double rsum = 0;

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	if ( !DBL_IS_EQUAL(array[i], missval) )
	  rsum += array[i];
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	rsum += array[i];
    }

  return (rsum);
}


double fldmean(FIELD field)
{
  int i;
  int    len      = field.size;
  int    nmiss    = field.nmiss;
  double missval1 = field.missval;
  double missval2 = field.missval;
  double *array   = field.ptr;
  double *w       = field.weight;
  double rsum = 0, rsumw = 0, ravg = 0;

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	if ( !DBL_IS_EQUAL(array[i], missval1) && !DBL_IS_EQUAL(w[i], missval1) )
	  {
	    rsum  += w[i] * array[i];
	    rsumw += w[i];
	  }
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	{
	  rsum  += w[i] * array[i];
	  rsumw += w[i];
	}
    }

  ravg = DIV(rsum, rsumw);

  return (ravg);
}


double fldavg(FIELD field)
{
  int i;
  int    len      = field.size;
  int    nmiss    = field.nmiss;
  double missval1 = field.missval;
  double missval2 = field.missval;
  double *array   = field.ptr;
  double *w       = field.weight;
  double rsum = 0, rsumw = 0, ravg = 0;

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	if ( !DBL_IS_EQUAL(w[i], missval1) )
	  {
	    rsum  = ADD(rsum, MUL(w[i], array[i]));
	    rsumw = ADD(rsumw, w[i]);
	  }
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	{
	  rsum  += w[i] * array[i];
	  rsumw += w[i];
	}
    }

  ravg = DIV(rsum, rsumw);

  return (ravg);
}


double fldvar(FIELD field)
{
  int i;
  int    len     = field.size;
  int    nmiss   = field.nmiss;
  double missval = field.missval;
  double *array  = field.ptr;
  double *w      = field.weight;
  double rsum = 0, rsumw = 0, rvar = 0;
  double rsumq = 0, rsumwq = 0;

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ ) 
        if ( !DBL_IS_EQUAL(array[i], missval) && !DBL_IS_EQUAL(w[i], missval) )
          {
            rsum   += w[i] * array[i];
            rsumq  += w[i] * array[i] * array[i];
            rsumw  += w[i];
            rsumwq += w[i] * w[i];
          }
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
        {
          rsum   += w[i] * array[i];
          rsumq  += w[i] * array[i] * array[i];
          rsumw  += w[i];
          rsumwq += w[i] * w[i];
        }
    }

  rvar = IS_NOT_EQUAL(rsumw, 0) ? (rsumq*rsumw - rsum*rsum) / (rsumw*rsumw) : missval;

  return (rvar);
}


double fldstd(FIELD field)
{
  double missval = field.missval;
  double rvar, rstd;

  rvar = fldvar(field);

  rstd = (IS_NOT_EQUAL(rvar, 0) && !DBL_IS_EQUAL(rvar, missval)) ? sqrt(rvar) : missval;

  return (rstd);
}


void fldrms(FIELD field, FIELD field2, FIELD *field3)
{
  int i, len, rnmiss = 0;
  int    grid1    = field.grid;
  int    nmiss1   = field.nmiss;
  double *array1  = field.ptr;
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double *array2  = field2.ptr;
  double missval1 = field.missval;
  double missval2 = field2.missval;
  double *w       = field.weight;
  double rsum = 0, rsumw = 0, ravg = 0;

  len    = gridInqSize(grid1);
  if ( len != gridInqSize(grid2) )
    cdoAbort("fields have different size!");

  /*
  if ( nmiss1 > 0 )
  */
    {
      for ( i = 0; i < len; i++ ) 
	if ( !DBL_IS_EQUAL(w[i], missval1) )
	  {
	    rsum  = ADD(rsum, MUL(w[i], MUL(SUB(array2[i], array1[i]),
                                            SUB(array2[i], array1[i]))));
	    rsumw = ADD(rsumw, w[i]);
	  }
    }
    /*
  else
    {
      for ( i = 0; i < len; i++ ) 
	{
	  rsum  += w[i] * array1[i];
	  rsumw += w[i];
	}
    }
    */

  ravg = ROOT(DIV(rsum, rsumw));

  if ( DBL_IS_EQUAL(ravg, missval1) ) rnmiss++;

  field3->ptr[0] = ravg;
  field3->nmiss  = rnmiss;
}


void varrms(FIELD field, FIELD field2, FIELD *field3)
{
  int i, k, nlev, len, rnmiss = 0;
  int    zaxis    = field.zaxis;
  int    grid1    = field.grid;
  int    nmiss1   = field.nmiss;
  double *array1  = field.ptr;
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double *array2  = field2.ptr;
  double missval1 = field.missval;
  double missval2 = field2.missval;
  double *w       = field.weight;
  double rsum = 0, rsumw = 0, ravg = 0;

  nlev   = zaxisInqSize(zaxis);
  len    = gridInqSize(grid1);
  if ( len != gridInqSize(grid2) )
    cdoAbort("fields have different size!");

  /*
  if ( nmiss1 > 0 )
  */
    {
      for ( k = 0; k < nlev; k++ )
	for ( i = 0; i < len; i++ )
	  /*	  if ( !DBL_IS_EQUAL(w[i], missval1) ) */
	    {
	      rsum  = ADD(rsum, MUL(w[i], MUL(SUB(array2[k*len+i], array1[k*len+i]),
                                              SUB(array2[k*len+i], array1[k*len+i]))));
	      rsumw = ADD(rsumw, w[i]);
	    }
    }
    /*
  else
    {
      for ( i = 0; i < len; i++ )
	{
	  rsum  += w[i] * array1[i];
	  rsumw += w[i];
	}
    }
    */

  ravg = ROOT(DIV(rsum, rsumw));

  if ( DBL_IS_EQUAL(ravg, missval1) ) rnmiss++;

  field3->ptr[0] = ravg;
  field3->nmiss  = rnmiss;
}

/* RQ */
double fldpctl(FIELD field, int p)
{
  static const char func[] = "fldpctl";
  
  int    len     = field.size;
  int    nmiss   = field.nmiss;
  double missval = field.missval;
  double *array  = field.ptr;
  double *array2;
	
  int i, j;
  double pctl = missval;
  
  if ( len - nmiss > 0 )
    {
      if ( nmiss > 0 )
        {
          array2 = (double *) malloc((len - nmiss)*sizeof(double));
          
          for ( i = 0, j = 0; i < len; i++ ) 
            if ( !DBL_IS_EQUAL(array[i], missval) )
              array2[j++] = array[i];

          pctl = nth_element(array2, j, (int)ceil(j*(p/100.0))-1);
          
          free(array2);	  
        }
      else
        {
          pctl = nth_element(array, len, (int)ceil(len*(p/100.0))-1);
        }
    }

  return pctl;
}
/* QR */
