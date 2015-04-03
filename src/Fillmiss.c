/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

*/

#include <cdi.h>
#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void fillmiss(field_t *field1, field_t *field2, int nfill)
{
  int gridID, nx, ny, i, j;
  int nmiss1, nmiss2 = 0;
  int kr, ku, kl, ko;
  int ir, iu, il, io;
  int kh, kv, k1, k2, kk;
  int globgrid = FALSE,gridtype;
  double s1, s2;
  double xr, xu, xl, xo;
  double missval;
  double *array1, *array2;
  double **matrix1, **matrix2;

  gridID   = field1->grid;
  nmiss1   = field1->nmiss;
  missval  = field1->missval;
  array1   = field1->ptr;
  array2   = field2->ptr;

  nx       = gridInqXsize(gridID);
  ny       = gridInqYsize(gridID);
  globgrid = gridIsCircular(gridID);

  gridtype = gridInqType(gridID);
  if ( !(gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN ) )
    cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

  matrix1 = (double **) malloc(ny * sizeof(double *));
  matrix2 = (double **) malloc(ny * sizeof(double *));

  for ( j = 0; j < ny; j++ )
    {
      matrix1[j] = array1 + j*nx;
      matrix2[j] = array2 + j*nx;
    }

  for ( j = 0; j < ny; j++ )
    for ( i = 0; i < nx; i++ )
      {
	if ( DBL_IS_EQUAL(matrix1[j][i], missval) )
	  {
	    nmiss2++;

	    kr = ku = kl = ko = 0;
	    xr = xu = xl = xo = 0.;

	    for ( ir = i + 1; ir < nx; ir++ )
	      if ( !DBL_IS_EQUAL(matrix1[j][ir], missval) )
		{ kr = ir - i; xr = matrix1[j][ir]; break; }

	    if ( globgrid && ir == nx )
	      {
		for ( ir = 0; ir < i; ir++ )
		  if ( !DBL_IS_EQUAL(matrix1[j][ir], missval) )
		    { kr = nx + ir - i; xr = matrix1[j][ir]; break; }
	      }

	    for ( il = i-1; il >= 0; il-- )
	      if ( !DBL_IS_EQUAL(matrix1[j][il], missval) )
		{ kl = i - il; xl = matrix1[j][il]; break; }

	    if ( globgrid && il == -1 )
	      {
		for ( il = nx-1; il > i; il-- )
		  if ( !DBL_IS_EQUAL(matrix1[j][il], missval) )
		    { kl = nx + i - il; xl = matrix1[j][il]; break; }
	      }

	    for ( iu = j + 1; iu < ny; iu++ )
	      if ( !DBL_IS_EQUAL(matrix1[iu][i], missval) )
		{ ku = iu - j; xu = matrix1[iu][i]; break; }
	    
	    for ( io = j - 1; io >= 0; io-- )
	      if ( !DBL_IS_EQUAL(matrix1[io][i], missval) )
		{ ko = j - io; xo = matrix1[io][i]; break; }
	    
	    /*  printf("%d %d %d %d %d %d %g %g %g %g\n", j,i,kr,kl,ku,ko,xr,xl,xu,xo);*/

	    kh = kl + kr;
	    kv = ko + ku;
	    if      ( kh == 0 ) { s1 = 0.; k1 = 0; }
	    else if ( kl == 0 ) { s1 = xr; k1 = 1; }
	    else if ( kr == 0 ) { s1 = xl; k1 = 1; }
	    else { s1 = xr*kl/kh + xl*kr/kh; k1 = 2; }

	    if      ( kv == 0 ) { s2 = 0.; k2 = 0; }
	    else if ( ku == 0 ) { s2 = xo; k2 = 1; }
	    else if ( ko == 0 ) { s2 = xu; k2 = 1; }
	    else { s2 = xu*ko/kv + xo*ku/kv; k2 = 2; }

	    kk = k1 + k2;
	    if ( kk >= nfill )
	      {
		if      ( kk == 0 ) cdoAbort("no point found!");
		else if ( k1 == 0 ) matrix2[j][i] = s2;
		else if ( k2 == 0 ) matrix2[j][i] = s1;
		else  matrix2[j][i] = s1*k2/kk + s2*k1/kk;
	      }
	    else
	      matrix2[j][i] = matrix1[j][i];

	    /* matrix1[j][i] = matrix2[j][i]; */
	  }
	else
	  {
	    matrix2[j][i] = matrix1[j][i];
	  }
      }

  if ( nmiss1 != nmiss2 ) cdoAbort("found only %d of %d missing values!", nmiss2, nmiss1);

  free(matrix2);
  free(matrix1);
}

void fillmiss_one_step(field_t *field1, field_t *field2, int maxfill)
{
  int gridID, nx, ny, i, j;
  int nmiss2 = 0;
  int kr, ku, kl, ko;
  int ir, iu, il, io;
  int kh, kv, k1, k2, kk;
  double s1, s2;
  double xr, xu, xl, xo;
  double missval;
  double *array1, *array2;
  double **matrix1, **matrix2;

  gridID  = field1->grid;
  missval = field1->missval;
  array1  = field1->ptr;
  array2  = field2->ptr;

  nx  = gridInqXsize(gridID);
  ny  = gridInqYsize(gridID);

  matrix1 = (double **) malloc(ny * sizeof(double *));
  matrix2 = (double **) malloc(ny * sizeof(double *));

  for ( j = 0; j < ny; j++ ) { matrix1[j] = array1 + j*nx; matrix2[j] = array2 + j*nx; }

  for (int fill_iterations=0; fill_iterations < maxfill; fill_iterations++) {
  for ( j = 0; j < ny; j++ )
    for ( i = 0; i < nx; i++ )
      {
        if ( DBL_IS_EQUAL(matrix1[j][i], missval) )
          {
            nmiss2++;

            kr = ku = kl = ko = 0;
            xr = xu = xl = xo = 0.;

            for ( ir = i + 1; ir < nx; ir++ )
              if ( !DBL_IS_EQUAL(matrix1[j][ir], missval) )
                { kr = ir - i; xr = matrix1[j][ir]; break; }

            for ( il = i-1; il >= 0; il-- )
              if ( !DBL_IS_EQUAL(matrix1[j][il], missval) )
                { kl = i - il; xl = matrix1[j][il]; break; }


            for ( iu = j + 1; iu < ny; iu++ )
              if ( !DBL_IS_EQUAL(matrix1[iu][i], missval) )
                { ku = iu - j; xu = matrix1[iu][i]; break; }

            for ( io = j - 1; io >= 0; io-- )
              if ( !DBL_IS_EQUAL(matrix1[io][i], missval) )
                { ko = j - io; xo = matrix1[io][i]; break; }


            kh = kl + kr;
            kv = ko + ku;
            if      ( kh == 0 ) { s1 = 0.; k1 = 0; }
            else if ( kl == 0 ) { s1 = xr; k1 = kr; }
            else if ( kr == 0 ) { s1 = xl; k1 = kl; }
            else
              {
                if ( kl < kr )
                {
                  s1 = xl;
                  k1 = kl;
                }
              else
                {
                  s1 = xr;
                  k1 = kr;
                }
              }

            if      ( kv == 0 ) { s2 = 0.; k2 = 0; }
            else if ( ku == 0 ) { s2 = xo; k2 = ko; }
            else if ( ko == 0 ) { s2 = xu; k2 = ku; }
            else
              {
                if ( ku < ko )
                  {
                    s2 = xu;
                    k2 = ku;
                  }
                else
                  {
                    s2 = xo;
                    k2 = ko;
                  }
              }

            kk = k1 + k2;
            if      ( kk == 0 ) matrix2[j][i] = matrix1[j][i];
            else if ( k1 == 0 ) matrix2[j][i] = s2;
            else if ( k2 == 0 ) matrix2[j][i] = s1;
            else
              {
                if ( k1 <= k2 )
                {
                  matrix2[j][i] = s1;
                }
                else
                {
                  matrix2[j][i] = s2;
                }

              }

            //printf("%d %d %2d %2d %2d %2d %2g %2g %2g %2g %2g %2g %2g\n", j,i,kr,kl,ku,ko,xr,xl,xu,xo,s1,s2,matrix2[j][i]);
            /* matrix1[j][i] = matrix2[j][i]; */
          }
        else
          {
            matrix2[j][i] = matrix1[j][i];
          }
      }
  for ( j = 0; j < ny; j++ ) for ( i = 0; i < nx; i++ ) matrix1[j][i] = matrix2[j][i];
  }

  free(matrix2);
  free(matrix1);
}


void *Fillmiss(void *argument)
{
  int FILLMISS,FILLMISSONESTEP;
  int operatorID;
  int streamID1, streamID2;
  int nrecs, ngrids;
  int index;
  int tsID, recID, varID, levelID;
  int gridsize;
  int vlistID1, vlistID2;
  int gridID1 = -1;
  field_t field1, field2;
  int taxisID1, taxisID2;
  int nmiss, i, nfill = 1;
  void (*fill_method) (field_t *fin , field_t *fout , int);

  cdoInitialize(argument);

  FILLMISS        = cdoOperatorAdd("fillmiss"  ,  0, 0, "nfill");
  FILLMISSONESTEP = cdoOperatorAdd("fillmiss2" ,  0, 0, "nfill");
  operatorID      = cdoOperatorID();
  if ( operatorID == FILLMISS )
     {
       fill_method = &fillmiss;
     }
  else
     {
       fill_method = &fillmiss_one_step;
     }

  /* Argument handling */
  {
    int oargc = operatorArgc();
    char **oargv = operatorArgv();

    if ( oargc == 1 )
      {
        nfill = atoi(oargv[0]);
        if ( operatorID != FILLMISSONESTEP ) 
          {
            if ( nfill < 1 || nfill > 4 ) cdoAbort("nfill out of range!");
          }
      }
    else if ( oargc > 1 )
      cdoAbort("Too many arguments!");
  }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);

      if ( gridInqType(gridID1) == GRID_GME ||
	   gridInqType(gridID1) == GRID_UNSTRUCTURED )
	cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridInqType(gridID1)) );
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);

  field_init(&field1);
  field_init(&field2);
  field1.ptr   = (double*) malloc(gridsize*sizeof(double));
  field2.ptr   = (double*) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field1.ptr, &field1.nmiss);

	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  field1.missval = vlistInqVarMissval(vlistID1, varID);

	  field2.grid    = field1.grid;
	  field2.nmiss   = 0;
	  field2.missval = field1.missval;

	  fill_method(&field1, &field2, nfill);

	  gridsize = gridInqSize(field2.grid);
	  nmiss = 0;
	  for ( i = 0; i < gridsize; ++i )
	    if ( DBL_IS_EQUAL(field2.ptr[i], field2.missval) ) nmiss++;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, field2.ptr, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field2.ptr ) free(field2.ptr);
  if ( field1.ptr ) free(field1.ptr);

  cdoFinish();

  return (0);
}
