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


void farfun(FIELD *field1, FIELD field2, int function)
{
  if      ( function == func_add   ) faradd(field1, field2);
  else if ( function == func_min   ) farmin(field1, field2);
  else if ( function == func_max   ) farmax(field1, field2);
  else if ( function == func_sum   ) farsum(field1, field2);
  else if ( function == func_mean  ) farsum(field1, field2);
  else if ( function == func_avg   ) faradd(field1, field2);
  else if ( function == func_sub   ) farsub(field1, field2);
  else if ( function == func_mul   ) farmul(field1, field2);
  else if ( function == func_div   ) fardiv(field1, field2);
  else if ( function == func_atan2 ) faratan2(field1, field2);
  else cdoAbort("function %d not implemented!", function);
}


void faradd(FIELD *field1, FIELD field2)
{
  static char func[] = "faradd";
  long   i, len;
  int    grid1    = field1->grid;
  int    nmiss1   = field1->nmiss;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  /*  double *weight1 = field1->weight; */
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] = ADD(array1[i], array2[i]);

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      /*
#if defined (_OPENMP)
#pragma omp parallel for default(shared)
#endif
      */
      for ( i = 0; i < len; i++ ) 
	array1[i] += array2[i];
    }
}


void farsum(FIELD *field1, FIELD field2)
{
  static char func[] = "farsum";
  long   i, len;
  int    grid1    = field1->grid;
  int    nmiss1   = field1->nmiss;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  /*  double *weight1 = field1->weight; */
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += array2[i];
	    else
	      array1[i] = array2[i];
	  }

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] += array2[i];
    }
}


void farsumq(FIELD *field1, FIELD field2)
{
  static char func[] = "farsumq";
  long   i, len;
  int    grid1    = field1->grid;
  int    nmiss1   = field1->nmiss;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  /*  double *weight1 = field1->weight; */
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += array2[i]*array2[i];
	    else
	      array1[i] = array2[i]*array2[i];
	  }

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] += array2[i]*array2[i];
    }
}


void farsub(FIELD *field1, FIELD field2)
{
  static char func[] = "farsub";
  long   i, len;
  int    grid1    = field1->grid;
  int    nmiss1   = field1->nmiss;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] = SUB(array1[i], array2[i]);

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] -= array2[i];
    }
}


void farmul(FIELD *field1, FIELD field2)
{
  static char func[] = "farmul";
  long   i, len;
  int    grid1    = field1->grid;
  int    nmiss1   = field1->nmiss;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] = MUL(array1[i], array2[i]);

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] *= array2[i];
    }
}


void fardiv(FIELD *field1, FIELD field2)
{
  static char func[] = "fardiv";
  long   i, len;
  int    grid1    = field1->grid;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  int    grid2    = field2.grid;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  for ( i = 0; i < len; i++ ) 
    array1[i] = DIV(array1[i], array2[i]);

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
}


void faratan2(FIELD *field1, FIELD field2)
{
  static char func[] = "fardiv";
  long   i, len;
  int    grid1    = field1->grid;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  int    grid2    = field2.grid;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  for ( i = 0; i < len; i++ ) 
    array1[i] = DBL_IS_EQUAL(array1[i],missval1) || DBL_IS_EQUAL(array2[i],missval2) ? missval1 : atan2(array1[i], array2[i]);

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
}


void farmin(FIELD *field1, FIELD field2)
{
  static char func[] = "farmin";
  long   i, len;
  int    grid1    = field1->grid;
  int    nmiss1   = field1->nmiss;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	{
	  array1[i] = DBL_IS_EQUAL(array2[i], missval2) ? array1[i] :
	              DBL_IS_EQUAL(array1[i], missval1) ? array2[i] :
		      MIN(array1[i], array2[i]);
	}

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ )
	array1[i] = MIN(array1[i], array2[i]);
    }
}


void farmax(FIELD *field1, FIELD field2)
{
  static char func[] = "farmax";
  long   i, len;
  int    grid1    = field1->grid;
  int    nmiss1   = field1->nmiss;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	{
	  array1[i] = DBL_IS_EQUAL(array2[i], missval2) ? array1[i] :
	              DBL_IS_EQUAL(array1[i], missval1) ? array2[i] :
		      MAX(array1[i], array2[i]);
	}

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ )
	array1[i] = MAX(array1[i], array2[i]);
    }
}


void farvar(FIELD *field1, FIELD field2, FIELD field3)
{
  static char func[] = "farstd";
  long   i, len;
  int    grid1    = field1->grid;
  int    nmiss1   = field1->nmiss;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;
  double *array3  = field3.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array1[i], missval1) && !DBL_IS_EQUAL(array2[i], missval2) )
	  array1[i] = array2[i]*array3[i] - (array1[i]*array3[i])*(array1[i]*array3[i]);
	else
	  array1[i] = missval1;
    }
  else
    {
      for ( i = 0; i < len; i++ )
	{
	  array1[i] = array2[i]*array3[i] - (array1[i]*array3[i])*(array1[i]*array3[i]);
	}
    }

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0.0 )
      {
	array1[i] = missval1;
	field1->nmiss++;
      }
}


void farstd(FIELD *field1, FIELD field2, FIELD field3)
{
  static char func[] = "farstd";
  long   i, len;
  int    grid1    = field1->grid;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  int    grid2    = field2.grid;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  farvar(field1, field2, field3);

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0.0 )
      {
	array1[i] = missval1;
	field1->nmiss++;
      }
    else
      {
	array1[i] = sqrt(array1[i]);
      }
}


void farcvar(FIELD *field1, FIELD field2, double rconst1)
{
  static char func[] = "farcvar";
  long   i, len;
  int    grid1    = field1->grid;
  int    nmiss1   = field1->nmiss;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array1[i], missval1) && !DBL_IS_EQUAL(array2[i], missval2) )
	  array1[i] = array2[i]*rconst1 - (array1[i]*rconst1)*(array1[i]*rconst1);
	else
	  array1[i] = missval1;
    }
  else
    {
      for ( i = 0; i < len; i++ )
	{
	  array1[i] = array2[i]*rconst1 - (array1[i]*rconst1)*(array1[i]*rconst1);
	}
    }

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0.0 )
      {
	array1[i] = missval1;
	field1->nmiss++;
      }
}


void farcstd(FIELD *field1, FIELD field2, double rconst1)
{
  static char func[] = "farcstd";
  long   i, len;
  int    grid1    = field1->grid;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  int    grid2    = field2.grid;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  farcvar(field1, field2, rconst1);

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0.0 )
      {
	array1[i] = missval1;
	field1->nmiss++;
      }
    else
      {
	array1[i] = sqrt(array1[i]);
      }
}


void farmoq(FIELD *field1, FIELD field2)
{
  static char func[] = "farmoq";
  long   i, len;
  int    grid1    = field1->grid;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  array1[i] = array2[i]*array2[i];
	else
	  array1[i] = missval1;

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] = array2[i]*array2[i];
    }
}


/* RQ */
/**
 * Counts the number of nonmissing values. The result of the operation
 * is computed according to the following rules:
 * 
 * field1  field2  result
 * a       b       a + 1
 * a       miss    a
 * miss    b       1
 * miss    miss    miss
 * 
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field
 */  
void farcount(FIELD *field1, FIELD field2)
{
  static char func[] = "farcount";
  long   i, len;
  int    grid1    = field1->grid;
  int    nmiss1   = field1->nmiss;
  double missval1 = field1->missval;
  double *array1  = field1->ptr;
  /*  double *weight1 = field1->weight; */
  int    grid2    = field2.grid;
  int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  len    = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += 1.0;
	    else
	      array1[i] = 1.0;
	  }

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] += 1.0;
    }
}
/* QR */
