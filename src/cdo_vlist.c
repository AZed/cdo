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

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "error.h"
#include "functs.h"

void vlistCompare(int vlistID1, int vlistID2, int function)
{
  static char func[] = "vlistCompare";
  int varID, nvars;

  if ( vlistNvars(vlistID1) != vlistNvars(vlistID2) )
    cdoAbort("Input streams have different number of variables per timestep!");

  if ( vlistNrecs(vlistID1) != vlistNrecs(vlistID2) )
    cdoAbort("Input streams have different number of records per timestep!");

  nvars = vlistNvars(vlistID1);

  if ( function == func_hrd )
    for ( varID = 0; varID < nvars; varID++ )
      {
	if ( vlistInqVarCode(vlistID1, varID) != vlistInqVarCode(vlistID2, varID) )
	  cdoAbort("Input streams have different structure!");

	if ( gridInqSize(vlistInqVarGrid(vlistID1, varID)) !=
	     gridInqSize(vlistInqVarGrid(vlistID2, varID)) )
	  cdoAbort("Grid size of the input fields do not match!");
      }
  else if ( function == func_code )
    for ( varID = 0; varID < nvars; varID++ )
      {
	if ( vlistInqVarCode(vlistID1, varID) != vlistInqVarCode(vlistID2, varID) )
	  cdoAbort("Input streams have different structure!");

	if ( zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) !=
	     zaxisInqSize(vlistInqVarZaxis(vlistID2, varID)) )
	  cdoAbort("Number of levels of the input fields do not match!");
      }
  else if ( function == func_sft || function == func_sftn /* || function == func_sftc */ )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( function == func_sftn )
	    {
	      char name1[256], name2[256];
	      vlistInqVarName(vlistID1, varID, name1);
	      vlistInqVarName(vlistID2, varID, name2);
	      if ( strcmp(name1, name2) != 0 )
		{
		  cdoWarning("Input streams have different variable names!");
		  break;
		}
	    }
	  /*
	  if ( function == func_sftc )
	    {
	      if ( vlistInqVarCode(vlistID1, varID) != vlistInqVarCode(vlistID2, varID) )
		{
		  cdoWarning("Input streams have different code numbers!");
		  break;
		}
	    }
	  */
	  if ( gridInqSize(vlistInqVarGrid(vlistID1, varID)) !=
	       gridInqSize(vlistInqVarGrid(vlistID2, varID)) )
	    cdoAbort("Grid size of the input fields do not match!");

	  if ( zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) !=
	       zaxisInqSize(vlistInqVarZaxis(vlistID2, varID)) )
	    cdoAbort("Number of levels of the input fields do not match!");
	}
      
      /* compare grids of first variable */
      {
	int gridID1, gridID2;
	int xsize, ysize;
	int i;

	gridID1 = vlistInqVarGrid(vlistID1, 0);
	gridID2 = vlistInqVarGrid(vlistID2, 0);

	if ( gridInqType(gridID1) == gridInqType(gridID2) )
	  {
	    if ( gridInqType(gridID1) == GRID_GAUSSIAN || gridInqType(gridID1) == GRID_LONLAT )
	      {
		xsize = gridInqXsize(gridID1);
		ysize = gridInqYsize(gridID1);
		
		if ( ysize == gridInqYsize(gridID2) )
		  {
		    if ( ysize > 1 )
		      {
			double *yvals1, *yvals2;

			yvals1 = (double *) malloc(ysize*sizeof(double));
			yvals2 = (double *) malloc(ysize*sizeof(double));

			gridInqYvals(gridID1, yvals1);
			gridInqYvals(gridID2, yvals2);
		
			if ( IS_EQUAL(yvals1[0], yvals2[ysize-1]) &&
			     IS_EQUAL(yvals1[ysize-1], yvals2[0]) )
			  {
			    if ( yvals1[0] > yvals2[0] )
			      cdoWarning("Grid orientation differ! First grid: N->S; second grid: S->N");
			    else
			      cdoWarning("Grid orientation differ! First grid: S->N; second grid: N->S");
			  }
			else
			  {
			    for ( i = 0; i < ysize; ++i )
			      if ( fabs(yvals1[i] - yvals2[i]) > 1.e-10 )
				{
				  cdoWarning("Grid latitudes differ!");
				  break;
				}
			  }

			free(yvals1);
			free(yvals2);
		      }
		  }
		else
		  cdoWarning("ysize of input grids differ!");
		
		if ( xsize == gridInqXsize(gridID2) )
		  {
		    if ( xsize > 1 )
		      {
			double *xvals1, *xvals2;

			xvals1 = (double *) malloc(xsize*sizeof(double));
			xvals2 = (double *) malloc(xsize*sizeof(double));

			gridInqXvals(gridID1, xvals1);
			gridInqXvals(gridID2, xvals2);
		
			for ( i = 0; i < xsize; ++i )
			  if ( fabs(xvals1[i] - xvals2[i]) > 1.e-10 )
			    {
			      cdoWarning("Grid longitudes differ!");
			      break;
			    }

			free(xvals1);
			free(xvals2);
		      }
		  }
		else
		  cdoWarning("ysize of input grids differ!");
	      }
	  }
	else
	  {
	    cdoWarning("Grids have different types! First grid: %s; second grid: %s",
		       gridNamePtr(gridInqType(gridID1)), gridNamePtr(gridInqType(gridID2)));
	  }
      }
    }
  else
    cdoAbort("Internal problem! Invalid function %d", function);
}


int vlistIsSzipped(int vlistID)
{
  int lszip = FALSE;
  int nvars, varID, ztype;

  nvars = vlistNvars(vlistID);

  for ( varID = 0; varID < nvars; varID++ )
    {						
      ztype = vlistInqVarZtype(vlistID, varID);
      if ( ztype == COMPRESS_SZIP )
	{
	  lszip = TRUE;
	  break;
	}
    }      

  return (lszip);
}


