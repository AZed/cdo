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

      Info       info            Dataset information
      Info       map             Dataset information and simple map
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void printMap(int nlon, int nlat, double *array, double missval, double min, double max)
{
  /* source code from PINGO */
  int ilon, ilat, i;
  double x, a, b;
  double step;
  double level[10];
  int min_n, max_n;
  char c;

  step = (max - min) / 10;

  if ( IS_NOT_EQUAL(step, 0) )
    {
      a = pow(10, floor(log(step) / M_LN10));
      b = step / a;

      if ( b > 5 )
	b = 0.5 * ceil (b / 0.5);
      else if ( b > 2 )
	b = 0.2 * ceil (b / 0.2);
      else if ( b > 1 )
	b = 0.1 * ceil (b / 0.1);
      else
	b = 1;

      step = b * a;

      if ( min < 0 && max > 0 )
	{
	  min_n = (int) floor (10 * (-min) / (max - min) - 0.5);
	  max_n = (int) ceil (10 * (-min) / (max - min) - 0.5);
	  level[min_n] = 0;
	  for (i = min_n - 1; i >= 0; i--)
	    level[i] = level[i + 1] - step;
	  for (i = max_n; i < 9; i++)
	    level[i] = level[i - 1] + step;
	}
      else
	{
	  level[0] = step * ceil (min / step + 0.5);
	  for ( i = 1; i < 9; i++ )
	    level[i] = level[i - 1] + step;
	}
    }
  else
    for ( i = 0; i < 9; i++ )
      level[i] = min;

  fputc ('\n', stdout);
  fflush (stdout);

  if ( nlon >= 1000 )
    {
      printf ("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
      for ( ilon = 0; ilon < nlon; ilon++ )
	printf ("%d", ((ilon + 1) / 1000) % 10);
      putchar ('\n');
      fflush (stdout);
    }

  if ( nlon >= 100 )
    {
      printf ("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
      for (ilon = 0; ilon < nlon; ilon++)
	printf ("%d", ((ilon + 1) / 100) % 10);
      putchar ('\n');
      fflush (stdout);
    }

  if ( nlon >= 10 )
    {
      printf ("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
      for (ilon = 0; ilon < nlon; ilon++)
	printf ("%d", ((ilon + 1) / 10) % 10);
      putchar ('\n');
      fflush (stdout);
    }

  printf ("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");

  for ( ilon = 0; ilon < nlon; ilon++ )
    printf ("%d", (ilon + 1) % 10);
  putchar ('\n');
  fflush (stdout);
  putchar ('\n');
  fflush (stdout);

  for ( ilat = 0; ilat < nlat; ilat++ )
    {
      printf ("%0*d ", nlat < 10 ? 1 : nlat < 100 ? 2 : nlat < 1000 ? 3 : 4, ilat + 1);
      for ( ilon = 0; ilon < nlon; ilon++ )
	{
	  x = array[ilat * nlon + ilon];
	  if ( DBL_IS_EQUAL(x, missval) )
	    c = '.';
	  else if ( DBL_IS_EQUAL(x, min) && !DBL_IS_EQUAL(min, max) )
	    c = 'm';
	  else if ( DBL_IS_EQUAL(x, max) && !DBL_IS_EQUAL(min, max) )
	    c = 'M';
	  else if ( DBL_IS_EQUAL(x, 0.) )
	    c = '*';
	  else if ( x < 0 )
	    {
	      c = '9';
	      for ( i = 0; i < 9; i++ )
		if ( level[i] > x )
		  {
		    c = i + '0';
		    break;
		  }
	    }
	  else
	    {
	      c = '0';
	      for ( i = 8; i >= 0; i-- )
		if ( level[i] < x )
		  {
		    c = i + 1 + '0';
		    break;
		  }
	    }
	  putchar (c);
	}
      printf (" %0*d\n", nlat < 10 ? 1 : nlat < 100 ? 2 : nlat < 1000 ? 3 : 4, ilat + 1);
      fflush (stdout);
    }
  putchar ('\n');
  fflush (stdout);

  if ( nlon >= 1000 )
    {
      printf ("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
      for ( ilon = 0; ilon < nlon; ilon++ )
	printf ("%d", ((ilon + 1) / 1000) % 10);
      putchar ('\n');
      fflush (stdout);
    }

  if ( nlon >= 100 )
    {
      printf ("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
      for ( ilon = 0; ilon < nlon; ilon++ )
	printf ("%d", ((ilon + 1) / 100) % 10);
      putchar ('\n');
      fflush (stdout);
    }

  if ( nlon >= 10 )
    {
      printf ("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
      for ( ilon = 0; ilon < nlon; ilon++ )
	printf ("%d", ((ilon + 1) / 10) % 10);
      putchar ('\n');
      fflush (stdout);
    }

  printf ("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
  for ( ilon = 0; ilon < nlon; ilon++ )
    printf ("%d", (ilon + 1) % 10);
  putchar ('\n');
  fflush (stdout);
  putchar ('\n');
  fflush (stdout);

  for ( i = 0; i < 10; i++ )
    {
      printf ("%d=%c%+9.3e,%+9.3e%c%s", (int) i,
	      i == 0 || level[i - 1] >= 0 ? '[' : '[',
	      i == 0 ? min : level[i - 1],
	      i == 9 ? max : level[i],
	      i == 9 || level[i] <= 0 ? ']' : ']',
	      i != 2 && i != 5 && i != 8 ? "  " : "");

      if ( i == 2 || i == 5 || i == 8 )
	{
	  fputc ('\n', stdout);
	  fflush (stdout);
	}
    }

  printf ("*=0  .=miss  m=min=%+9.3e  M=max=%+9.3e\n", min, max);
  fflush (stdout);
  putchar ('\n');
  fflush (stdout);
}


void *Info(void *argument)
{
  int INFO, INFOP, INFON, INFOC, MAP;
  int operatorID;
  int i;
  int indf, indg;
  int varID, recID;
  int gridsize = 0;
  int gridID, zaxisID;
  int code, param;
  int vdate, vtime;
  int nrecs;
  int levelID;
  int tsID, taxisID;
  int streamID = 0;
  int vlistID;
  int nmiss;
  int number;
  int ivals = 0, nvals = 0;
  int imiss = 0;
  char varname[CDI_MAX_NAME];
  char paramstr[32];
  char vdatestr[32], vtimestr[32];
  double missval;
  double *array = NULL;
  double level;
  double arrmin = 0, arrmax = 0, arrmean = 0, arrvar = 0;

  cdoInitialize(argument);

  INFO  = cdoOperatorAdd("info",  0, 0, NULL);
  INFOP = cdoOperatorAdd("infop", 0, 0, NULL);
  INFON = cdoOperatorAdd("infon", 0, 0, NULL);
  INFOC = cdoOperatorAdd("infoc", 0, 0, NULL);
  MAP   = cdoOperatorAdd("map",   0, 0, NULL);

  operatorID = cdoOperatorID();

  for ( indf = 0; indf < cdoStreamCnt(); indf++ )
    {
      streamID = streamOpenRead(cdoStreamName(indf));

      vlistID = streamInqVlist(streamID);

      if ( vlistNvars(vlistID) == 0 ) continue;

      gridsize = vlistGridsizeMax(vlistID);
      if ( vlistNumber(vlistID) != CDI_REAL ) gridsize *= 2;

      array = malloc(gridsize*sizeof(double));

      indg = 0;
      tsID = 0;
      taxisID = vlistInqTaxis(vlistID);
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{
	  vdate = taxisInqVdate(taxisID);
	  vtime = taxisInqVtime(taxisID);

	  date2str(vdate, vdatestr, sizeof(vdatestr));
	  time2str(vtime, vtimestr, sizeof(vtimestr));

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      if ( (tsID == 0 && recID == 0) || operatorID == MAP )
		{
		  fprintf(stdout, "%6d :       Date     Time   Level Gridsize    Miss :"
			  "     Minimum        Mean     Maximum : ",  -(indf+1));

		  if      ( operatorID == INFON ) fprintf(stdout, "Parameter name");
		  else if ( operatorID == INFOC ) fprintf(stdout, "Code number");
		  else                            fprintf(stdout, "Parameter ID");

		  if ( cdoVerbose ) fprintf(stdout, " : Extra" );              
		  fprintf(stdout, "\n" );              
		}

	      streamInqRecord(streamID, &varID, &levelID);
	      streamReadRecord(streamID, array, &nmiss);

	      indg += 1;
	      param    = vlistInqVarParam(vlistID, varID);
	      code     = vlistInqVarCode(vlistID, varID);
	      gridID   = vlistInqVarGrid(vlistID, varID);
	      zaxisID  = vlistInqVarZaxis(vlistID, varID);
	      missval  = vlistInqVarMissval(vlistID, varID);
	      gridsize = gridInqSize(gridID);
	      number   = vlistInqVarNumber(vlistID, varID);

	      cdiParamToString(param, paramstr, sizeof(paramstr));

	      if ( operatorID == INFON ) vlistInqVarName(vlistID, varID, varname);

	      fprintf(stdout, "%6d :%s %s ", indg, vdatestr, vtimestr);

	      level = zaxisInqLevel(zaxisID, levelID);
	      fprintf(stdout, "%7g ", level);

	      fprintf(stdout, "%8d %7d :", gridsize, nmiss);

	      if ( /* gridInqType(gridID) == GRID_SPECTRAL || */
		   (gridsize == 1 && nmiss == 0 && number == CDI_REAL) )
		{
		  fprintf(stdout, "            %#12.5g            ", array[0]);
		}
	      else
		{
		  if ( number == CDI_REAL )
		    {
		      if ( nmiss > 0 )
			{
			  ivals   = 0;
			  arrmean = 0;
			  arrvar  = 0;
			  arrmin  =  1.e300;
			  arrmax  = -1.e300;
			  for ( i = 0; i < gridsize; ++i )
			    {
			      if ( !DBL_IS_EQUAL(array[i], missval) )
				{
				  if ( array[i] < arrmin ) arrmin = array[i];
				  if ( array[i] > arrmax ) arrmax = array[i];
				  arrmean += array[i];
				  arrvar  += array[i]*array[i];
				  ivals++;
				}
			    }
			  imiss = gridsize - ivals;
			  nvals = ivals;
			}
		      else
			{
			  arrmean = array[0];
			  arrvar  = array[0];
			  arrmin  = array[0];
			  arrmax  = array[0];
			  /*
#pragma omp parallel for default(none) shared(arrmin, arrmax, array, gridsize)	\
                                       reduction(+:arrmean, arrvar)
			  */
			  for ( i = 1; i < gridsize; i++ )
			    {
			      /* #pragma omp critical */
			      if ( array[i] < arrmin ) arrmin = array[i];
			      /* #pragma omp critical */
			      if ( array[i] > arrmax ) arrmax = array[i];
			      arrmean += array[i];
			      arrvar  += array[i]*array[i];
			    }
			  nvals = gridsize;
			}

		      if ( nvals )
			{
			  arrmean = arrmean/nvals;
			  arrvar  = arrvar/nvals - arrmean*arrmean;
			  fprintf(stdout, "%#12.5g%#12.5g%#12.5g", arrmin, arrmean, arrmax);
			}
		      else
			{
			  fprintf(stdout, "                     nan            ");
			}
		    }
		  else
		    {
		      int nvals_r = 0, nvals_i = 0;
		      double arrsum_r, arrsum_i, arrmean_r = 0, arrmean_i = 0;
		      arrsum_r = 0;
		      arrsum_i = 0;
		      
		      for ( i = 0; i < gridsize; i++ )
			{
			  if ( !DBL_IS_EQUAL(array[i*2],   missval) && 
			       !DBL_IS_EQUAL(array[i*2+1], missval) )
			    {
			      arrsum_r += array[i*2];
			      arrsum_i += array[i*2+1];
			      nvals_r++;
			      nvals_i++;
			    }
			}

		      imiss = gridsize - nvals_r;

		      if ( nvals_r > 0 ) arrmean_r = arrsum_r / nvals_r;
		      if ( nvals_i > 0 ) arrmean_i = arrsum_i / nvals_i;
		      fprintf(stdout, "   -  (%#12.5g,%#12.5g)  -", arrmean_r, arrmean_i);
		    }
		}

	      if ( operatorID == INFON )
		fprintf(stdout, " : %-14s", varname);
	      else if ( operatorID == INFOC )
		fprintf(stdout, " : %4d   ", code);
	      else
		fprintf(stdout, " : %-14s", paramstr);

	      if ( cdoVerbose )
		{
		  char varextra[CDI_MAX_NAME];
		  vlistInqVarExtra(vlistID, varID, varextra);
		  fprintf(stdout, " : %s", varextra );              
		}

	      fprintf(stdout, "\n");

	      if ( imiss != nmiss && nmiss > 0 )
		fprintf(stdout, "Found %d of %d missing values!\n", imiss, nmiss);

	      if ( operatorID == MAP )
		{
		  int nlon, nlat;
		  
		  nlon = gridInqXsize(gridID);
		  nlat = gridInqYsize(gridID);

		  if ( gridInqType(gridID) == GRID_GAUSSIAN    ||
		       gridInqType(gridID) == GRID_LONLAT      ||
		       gridInqType(gridID) == GRID_CURVILINEAR ||
		       (gridInqType(gridID) == GRID_GENERIC && 
			nlon*nlat == gridInqSize(gridID) && nlon < 1024) )
		    {
		      printMap(nlon, nlat, array, missval, arrmin, arrmax);
		    }
		}
	    }
	  tsID++;
	}
      streamClose(streamID);

      if ( array ) free(array);
    }

  cdoFinish();

  return (0);
}
