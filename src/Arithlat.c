/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Arithlat   mulcoslat       Multiply with cos(lat)
      Arithlat   divcoslat       Divide by cos(lat)
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#include <math.h>

#ifndef  M_PI
#define  M_PI		3.14159265358979323846	/* pi */
#endif

#ifndef  DEG2RAD
#define  DEG2RAD  (M_PI/180.)   /* conversion for deg to rad */
#endif

void *Arithlat(void *argument)
{
  static char func[] = "Arithlat";
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int gridsize, gridtype;
  int gridID, gridID0 = -1;
  int nlon = 0, nlat = 0, i, j;
  int nrecs, recID;
  int tsID;
  int varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int nmiss;
  double *scale = NULL;
  double *array = NULL;

  cdoInitialize(argument);

  cdoOperatorAdd("mulcoslat", func_mul, 0, NULL);
  cdoOperatorAdd("divcoslat", func_div, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operfunc == func_mul || operfunc == func_div )
    nospec(vlistID1);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);

  array = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array, &nmiss);
	  
	  gridID = vlistInqVarGrid(vlistID1, varID);

	  if ( gridID != gridID0 )
	    {
	      gridtype = gridInqType(gridID);
	      if ( gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN )
		{
		  if ( gridInqType(gridID) == GRID_GAUSSIAN_REDUCED )
		    cdoAbort("Gaussian reduced grid found. Use option -R to convert it to a regular grid!");
		  else
		    cdoAbort("LONLAT or GAUSSIAN grid not found!");
		}

	      gridsize = gridInqSize(gridID);
	      nlon = gridInqXsize(gridID);
	      nlat = gridInqYsize(gridID);

	      scale = (double *) realloc(scale, nlat*sizeof(double));
	      gridInqYvals(gridID, scale);

	      if ( operfunc == func_mul )
		for ( j = 0; j < nlat; j++ ) scale[j] = cos(scale[j]*DEG2RAD);
	      else
		for ( j = 0; j < nlat; j++ ) scale[j] = 1./cos(scale[j]*DEG2RAD);

	      if ( cdoVerbose ) for ( j = 0; j < nlat; j++ ) cdoPrint("coslat  %3d  %g", j+1, scale[j]);
		  
	      gridID0 = gridID;
	    }

	  for ( j = 0; j < nlat; j++ )
	    for ( i = 0; i < nlon; i++ )
	      array[i+j*nlon] *= scale[j];

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);
  if ( scale ) free(scale);

  cdoFinish();

  return (0);
}
