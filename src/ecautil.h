/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef ECAUTIL_H_
#define ECAUTIL_H_

#include "field.h"

/**
 * Computes the day-of-year correspnding a given Gregorian date.
 * 
 * @param date a Gregorian date in the form YYYYMMDD
 * 
 * @return the day-of-year
 */
unsigned long day_of_year(int date);

/**
 * Counts the number of nonmissing values in a field. The result
 * of the operation is computed according to the following rules:
 * 
 * field1  field2  result
 * a       b       a + 1
 * a       miss    a
 * miss    b       1
 * miss    miss    0
 * 
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field  
 */  
void farnum(FIELD *field1, FIELD field2);

/**
 * Counts the number of consecutive nonmissing values in a field.
 * The result of the operation is computed according to the following
 * rules:
 * 
 * field1  field2  result
 * a       b       a + 1
 * a       miss    0
 * miss    b       1
 * miss    miss    0
 * 
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field  
 */  
void farnum2(FIELD *field1, FIELD field2);

/**
 * Counts the number of values in series of at least n consecutive
 * nonmissing values. The result of the operation is computed according
 * to the following rules:
 * 
 * field1  field2  result
 * a       b       b < n ? a : b > n ? a + 1 : a + n
 * a       miss    a
 * miss    b       b < n ? 0 : b
 * miss    miss    0    
 * 
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field
 * @param n      the number of consecutive values, must be an exact
 *               mathematical integer
 */  
void farnum3(FIELD *field1, FIELD field2, double n);

/**
 * Selects field elements according to a given mask. The result of
 * the operation is computed according to the following rules:
 * 
 * field1  field2  result
 * a       b       b != 0.0 ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss    
 * 
 * @param field1  the input field, also holds the result
 * @param field2  the mask
 */  
void farsel(FIELD *field1, FIELD field2);

/**
 * Selects all field elements that are less than or equal to the
 * corresponding element of a reference field. The result of the
 * operation is computed according to the following rules:
 * 
 * field1  field2  result
 * a       b       a <= b ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss    
 * 
 * @param field1 the input field, also holds the result
 * @param field2 the reference field
 */  
void farselle(FIELD *field1, FIELD field2);

/**
 * Selects all field elements that are less than the
 * corresponding element of a reference field. The result of the
 * operation is computed according to the following rules:
 * 
 * field1  field2  result
 * a       b       a < b ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss    
 * 
 * @param field1 the input field, also holds the result
 * @param field2 the reference field
 */  
void farsellt(FIELD *field1, FIELD field2);

/**
 * Selects all field elements that are greater than or equal to
 * the corresponding element of a reference field. The result of
 * the operation is computed according to the following rules:
 * 
 * field1  field2  result
 * a       b       a >= b ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss    
 * 
 * @param field1 the input field, also holds the result
 * @param field2 the reference field
 */  
void farselge(FIELD *field1, FIELD field2);

/**
 * Selects all field elements that are greater than the
 * corresponding element of a reference field. The result of the
 * operation is computed according to the following rules:
 * 
 * field1  field2  result
 * a       b       a > b ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss    
 * 
 * @param field1 the input field, also holds the result
 * @param field2 the reference field
 */  
void farselgt(FIELD *field1, FIELD field2);

/**
 * Selects all field elements that are equal to the
 * corresponding element of a reference field. The result of the
 * operation is computed according to the following rules:
 * 
 * field1  field2  result
 * a       b       a == b ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss    
 * 
 * @param field1 the input field, also holds the result
 * @param field2 the reference field
 */  
void farseleq(FIELD *field1, FIELD field2);

/**
 * Selects all field elements that are not equal to the
 * corresponding element of a reference field. The result of the
 * operation is computed according to the following rules:
 * 
 * field1  field2  result
 * a       b       a != b ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss    
 * 
 * @param field1 the input field, also holds the result
 * @param field2 the reference field
 */  
void farselne(FIELD *field1, FIELD field2);

/**
 * Selects all field elements that are less than or equal to a
 * certain reference value. The result of the operation is computed
 * according to the following rules:
 * 
 * field  c      result
 * a      c      a <= c ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss    
 * 
 * @param field the input field, also holds the result
 * @param c     the reference value
 */  
void farsellec(FIELD *field, double c);

/**
 * Selects all field elements that are less a
 * certain reference value. The result of the operation is computed
 * according to the following rules:
 * 
 * field  c      result
 * a      c      a < c ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss    
 * 
 * @param field the input field, also holds the result
 * @param c     the reference value
 */  
void farselltc(FIELD *field, double c);

/**
 * Selects all field elements that are greater than or equal to a
 * certain reference value. The result of the operation is computed
 * according to the following rules:
 * 
 * field  c      result
 * a      c      a >= c ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss    
 * 
 * @param field the input field, also holds the result
 * @param c     the reference value
 */  
void farselgec(FIELD *field, double c);

/**
 * Selects all field elements that are greater than a
 * certain reference value. The result of the operation is computed
 * according to the following rules:
 * 
 * field  c      result
 * a      c      a > c ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss    
 * 
 * @param field the input field, also holds the result
 * @param c     the reference value
 */  
void farselgtc(FIELD *field, double c);

/**
 * Selects all field elements that are equal to a
 * certain reference value. The result of the operation is computed
 * according to the following rules:
 * 
 * field  c      result
 * a      c      a == c ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss    
 * 
 * @param field the input field, also holds the result
 * @param c     the reference value
 */  
void farseleqc(FIELD *field, double c);

/**
 * Selects all field elements that are not equal to a
 * certain reference value. The result of the operation is computed
 * according to the following rules:
 * 
 * field  c      result
 * a      c      a != c ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss    
 * 
 * @param field the input field, also holds the result
 * @param c     the reference value
 */  
void farselnec(FIELD *field, double c);

/**
 * reset the fields real values to the missval for all levels
 *
 * @param field the input field, also holds the result
 * @param nlevels   number of available levels
 * @param gridsize  number of grid points
 */
void updateHist(FIELD *field[2], int nlevels, int gridsize, double *yvals, int onlyNorth);

/*
 * Compute the Gsl and its starting day
 *
 * @param int nlevels
 * @param int gridsize
 * @param double *yvals = array of latitudes
 * @param int ysize = number of gridpoints in lat-direction
 * @param double missval
 * @param int ovdate = the last year, which has been fully processed
 * @param FIELD *startDate
 * @param FIELD *endDate
 * @param FIELD *startDateWithHist[2]
 * @param FIELD *endDateWithHist[2]
 * @param FIELD *gslDuration
 * @param FIELD *gslFirstDay
 * @param int useCurrentYear = if TRUE, only measurements of the current year
 *                             (index 0) are used for computation, i.e. that
 *                             gsl can only be computed for the northern
 *                             hemisphere (see definition of GSL: EcaGsl()
 */
void computeGsl(int nlevels, int gridsize, double *yvals, double missval, 
                FIELD *startDateWithHist[2], FIELD *endDateWithHist[2],
                FIELD *gslDuration, FIELD *gslFirstDay,
                int useCurrentYear);

/*
 * Adjust the endDates found in the current year:
 * if a starting date for gsl could be found, but no ending date, the end
 * should be the last day of the corresponding year for norther and June, 30th
 * for southern hemisphere
 */
void adjustEndDate(int nlevels, int gridsize, double *yvals, double missval, int ovdate,
                   FIELD *startDateWithHist[2], FIELD *endDateWithHist[2]);
/*
 * Write GSL related fields to an output stream
 */
void writeGslStream(int ostreamID, int otaxisID, int otsID, 
                    int ovarID1, int ovarID2, int ivlistID1,
                    int first_var_id,
                    FIELD *gslDuration, FIELD *gslFirstDay,
                    int vdate, int vtime,  int nlevels);

#endif /*ECAUTIL_H_*/
