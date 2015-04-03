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

      Vargen     const           Create a constant field
      Vargen     random          Field with random values
*/

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#if defined (__GNUC__)
#if __GNUC__ > 2
#  define WITH_ETOPO 1
#endif
#endif


void *Vargen(void *argument)
{
  static char func[] = "Vargen";
  int RANDOM, CONST, TOPO, FOR;
  int operatorID;
  int streamID;
  int nrecs, ntimesteps;
  int tsID, recID, varID, levelID;
  int gridsize, i;
  int vlistID;
  int gridID = -1, zaxisID, taxisID;
  int vdate, vtime, julday;
  const char *gridfile;
  double rval, rstart = 0, rstop = 0, rinc = 0;
  double rconst = 0;
  double *array;
#if defined(WITH_ETOPO)
  double etopo_scale = 3;
  static const short etopo[] = {
#include "etopo.h"
  };
#endif

  cdoInitialize(argument);

  RANDOM = cdoOperatorAdd("random", 0, 0, "grid description file or name, <seed>");
  CONST  = cdoOperatorAdd("const",  0, 0, "constant value, grid description file or name");
  TOPO   = cdoOperatorAdd("topo",   0, 0, "");
  FOR    = cdoOperatorAdd("for",    0, 0, "start, end<, increment>");

  operatorID = cdoOperatorID();

  if ( operatorID == RANDOM )
    {
      unsigned int seed = 1;
      operatorInputArg(cdoOperatorEnter(operatorID));
      if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");
      if ( operatorArgc() > 2 ) cdoAbort("Too many arguments!");
      gridfile = operatorArgv()[0];
      gridID   = cdoDefineGrid(gridfile);
      if ( operatorArgc() == 2 )
	{
	  long idum;
	  idum = atol(operatorArgv()[1]);
	  if ( idum >= 0 && idum < 0x7FFFFFFF )
	    seed = idum;
	}
      srand(seed);
    }
  else if ( operatorID == CONST )
    {
      operatorInputArg(cdoOperatorEnter(operatorID));
      operatorCheckArgc(2);
      rconst   = atof(operatorArgv()[0]);
      gridfile = operatorArgv()[1];
      gridID   = cdoDefineGrid(gridfile);
    }
  else if ( operatorID == TOPO )
    {
      int nlon, nlat, i;
      double lon[720], lat[360];
      nlon = 720;
      nlat = 360;
      gridID = gridCreate(GRID_LONLAT, nlon*nlat);
      gridDefXsize(gridID, nlon);
      gridDefYsize(gridID, nlat);

      for ( i = 0; i < nlon; i++ ) lon[i] = -179.75 + i*0.5;
      for ( i = 0; i < nlat; i++ ) lat[i] = -89.75 + i*0.5;

      gridDefXvals(gridID, lon);
      gridDefYvals(gridID, lat);
    }
  else if ( operatorID == FOR )
    {
      operatorInputArg(cdoOperatorEnter(operatorID));
      if ( operatorArgc() < 2 ) cdoAbort("Too few arguments!");
      if ( operatorArgc() > 3 ) cdoAbort("Too many arguments!");

      rstart = atof(operatorArgv()[0]);
      rstop  = atof(operatorArgv()[1]);
      if ( operatorArgc() == 3 )
	rinc = atof(operatorArgv()[2]);
      else
	rinc = 1;

      if ( DBL_IS_EQUAL(rinc, 0.0) ) cdoAbort("Increment is zero!");

      gridID = gridCreate(GRID_GENERIC, 1);
      gridDefXsize(gridID, 1);
      gridDefYsize(gridID, 1);
    }


  zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  vlistID = vlistCreate();

  if ( operatorID == FOR )
    varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
  else
    varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_CONSTANT);

  taxisID = taxisCreate(TAXIS_RELATIVE);
  vlistDefTaxis(vlistID, taxisID);

  if ( operatorID == RANDOM || operatorID == CONST || operatorID == TOPO )
    vlistDefNtsteps(vlistID, 1);

  streamID = streamOpenWrite(cdoStreamName(0), cdoFiletype());
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));

  streamDefVlist(streamID, vlistID);

  gridsize = gridInqSize(gridID);
  array = (double *) malloc(gridsize*sizeof(double));

  if ( operatorID == FOR )
    ntimesteps = 1.001 + ((rstop-rstart)/rinc);
  else
    ntimesteps = 1;

  julday = date_to_julday(CALENDAR_PROLEPTIC, 10101);

  for ( tsID = 0; tsID < ntimesteps; tsID++ )
    {
      rval  = rstart + rinc*tsID;
      vdate = julday_to_date(CALENDAR_PROLEPTIC, julday + tsID);
      vtime = 0;
      taxisDefVdate(taxisID, vdate);
      taxisDefVtime(taxisID, vtime);
      streamDefTimestep(streamID, tsID);

      nrecs = 1;
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  levelID = 0;
	  streamDefRecord(streamID, varID, levelID);

	  if ( operatorID == RANDOM )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array[i] = rand()/(RAND_MAX+1.0);
	    }
	  else if ( operatorID == CONST )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array[i] = rconst;
	    }
	  else if ( operatorID == TOPO )
	    {
#if defined(WITH_ETOPO)
	      for ( i = 0; i < gridsize; i++ )
		array[i] = (double)etopo[i]/etopo_scale;
#else
	      cdoAbort("Operator support disabled!");
#endif
	    }
	  else if ( operatorID == FOR )
	    {
	      array[0] = rval;
	    }

	  streamWriteRecord(streamID, array, 0);
	}
    }

  streamClose(streamID);

  vlistDestroy(vlistID);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}
