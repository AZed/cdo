/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2013 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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
      Vargen     stdatm          Field values for pressure and temperature for
                                 the standard atmosphere
*/


#if defined(HAVE_CONFIG_H)
#  include "config.h" // ENABLE_DATA
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "list.h"


#if defined(ENABLE_DATA)
  static double etopo_scale  = 3;
  static double etopo_offset = 11000;
  static const unsigned short etopo[] = {
#include "etopo.h"
  };

  static double temp_scale  =  500;
  static double temp_offset = -220;
  static const unsigned short temp[] = {
#include "temp.h"
  };

  static double mask_scale  =  1;
  static double mask_offset =  0;
  static const unsigned short mask[] = {
#include "mask.h"
  };
#endif

/*  some Constants for creating temperatur and pressure for the standard atmosphere */
#define T_ZERO          (213.0)
#define T_DELTA          (75.0)
#define SCALEHEIGHT   (10000.0)   /* [m] */
#define P_ZERO         (1013.25)  /* surface pressure [hPa] */
#define C_EARTH_GRAV      (9.80665)
#define C_R             (287.05)  /*  specific gas constant for air */
static double TMP4PRESSURE = (C_EARTH_GRAV*SCALEHEIGHT)/(C_R*T_ZERO);

static double
std_atm_temperatur(double height)
{
  /*
    Compute the temperatur for the given height (in meters) according to the
    solution of the hydrostatic atmosphere
   */
   return (T_ZERO + T_DELTA * exp((-1)*(height/SCALEHEIGHT)));
}

static double
std_atm_pressure(double height)
{
  /*
    Compute the pressure for the given height (in meters) according to the
    solution of the hydrostatic atmosphere
   */
  return (P_ZERO * exp((-1)*TMP4PRESSURE*log((exp(height/SCALEHEIGHT)*T_ZERO + T_DELTA)/(T_ZERO + T_DELTA))));
}

void *Vargen(void *argument)
{
  int RANDOM, SINCOS, CONST, FOR, TOPO, TEMP, MASK, STDATM;
  int operatorID;
  int streamID;
  int nvars, ntimesteps, nlevels = 1;
  int tsID, varID, varID2 = -1, levelID;
  int gridsize, i;
  int vlistID;
  int gridID = -1, zaxisID, taxisID;
  int vdate, vtime, julday;
  const char *gridfile;
  double rval, rstart = 0, rstop = 0, rinc = 0;
  double rconst = 0;
  double *array, *levels = NULL;

  cdoInitialize(argument);

  RANDOM = cdoOperatorAdd("random", 0, 0, "grid description file or name, <seed>");
  SINCOS = cdoOperatorAdd("sincos", 0, 0, "grid description file or name");
  CONST  = cdoOperatorAdd("const",  0, 0, "constant value, grid description file or name");
  FOR    = cdoOperatorAdd("for",    0, 0, "start, end, <increment>");
  TOPO   = cdoOperatorAdd("topo",   0, 0, NULL);
  TEMP   = cdoOperatorAdd("temp",   0, 0, NULL);
  MASK   = cdoOperatorAdd("mask",   0, 0, NULL);
  STDATM = cdoOperatorAdd("stdatm", 0, 0, "levels");

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
  else if ( operatorID == SINCOS )
    {
      operatorInputArg(cdoOperatorEnter(operatorID));
      operatorCheckArgc(1);
      gridfile = operatorArgv()[0];
      gridID   = cdoDefineGrid(gridfile);
    }
  else if ( operatorID == CONST )
    {
      operatorInputArg(cdoOperatorEnter(operatorID));
      operatorCheckArgc(2);
      rconst   = atof(operatorArgv()[0]);
      gridfile = operatorArgv()[1];
      gridID   = cdoDefineGrid(gridfile);
    }
  else if ( operatorID == TOPO || operatorID == TEMP || operatorID == MASK )
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
      double lon = 0, lat = 0;
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

      gridID = gridCreate(GRID_LONLAT, 1);
      gridDefXsize(gridID, 1);
      gridDefYsize(gridID, 1);
      gridDefXvals(gridID, &lon);
      gridDefYvals(gridID, &lat);
    }
  else if ( operatorID == STDATM )
    {
      double lon = 0, lat = 0;
      LIST *flist = listNew(FLT_LIST);

      operatorInputArg("levels");
      nlevels = args2fltlist(operatorArgc(), operatorArgv(), flist);
      levels  = (double *) listArrayPtr(flist);
      //listDelete(flist);

      if ( cdoVerbose ) for ( i = 0; i < nlevels; ++i ) printf("levels %d: %g\n", i, levels[i]);

      gridID = gridCreate(GRID_LONLAT, 1);
      gridDefXsize(gridID, 1);
      gridDefYsize(gridID, 1);
      gridDefXvals(gridID, &lon);
      gridDefYvals(gridID, &lat);
    }

  if ( operatorID == STDATM )
    {
      zaxisID = zaxisCreate(ZAXIS_HEIGHT, nlevels);
      zaxisDefLevels(zaxisID  , levels);
      zaxisDefName(zaxisID    , "level");
      zaxisDefLongname(zaxisID, "Level");
      zaxisDefUnits(zaxisID   , "m");
    }
  else
    {
      zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
      nlevels = 1;
    }

  vlistID = vlistCreate();

  if ( operatorID == FOR )
    varID = vlistDefVar(vlistID, gridID, zaxisID, TSTEP_INSTANT);
  else
    varID = vlistDefVar(vlistID, gridID, zaxisID, TSTEP_CONSTANT);
  /*
     For the standard atmosphere two output variables are generated: pressure and
     temperatur. The first (varID) is pressure, second (varID2) is temperatur.
     Add an additional variable for the standard atmosphere.
   */
  if ( operatorID == STDATM )
    varID2 = vlistDefVar(vlistID, gridID, zaxisID, TSTEP_CONSTANT);

  if ( operatorID == MASK )
    vlistDefVarDatatype(vlistID, varID, DATATYPE_INT8);

  if ( operatorID == STDATM )
    {
      vlistDefVarName(vlistID    , varID , "P");
      vlistDefVarCode(vlistID    , varID , 1);
      vlistDefVarStdname(vlistID , varID , "air_pressure");
      vlistDefVarLongname(vlistID, varID , "pressure");
      vlistDefVarUnits(vlistID   , varID , "hPa");
      vlistDefVarName(vlistID    , varID2, "T");
      vlistDefVarCode(vlistID    , varID2, 130);
      vlistDefVarStdname(vlistID , varID2, "air_temperature");
      vlistDefVarLongname(vlistID, varID2, "temperature");
      vlistDefVarUnits(vlistID   , varID2, "K");
    }
  else
    {
      vlistDefVarName(vlistID, varID, cdoOperatorName(operatorID));
      if ( operatorID == TOPO )
	vlistDefVarUnits(vlistID, varID , "m");	
      if ( operatorID == TEMP )
	vlistDefVarUnits(vlistID, varID , "K");	
    }

  taxisID = taxisCreate(TAXIS_RELATIVE);
  vlistDefTaxis(vlistID, taxisID);

  if ( operatorID == RANDOM || operatorID == SINCOS || operatorID == CONST || operatorID == TOPO ||
       operatorID == TEMP || operatorID == MASK || operatorID == STDATM )
    vlistDefNtsteps(vlistID, 1);

  streamID = streamOpenWrite(cdoStreamName(0), cdoFiletype());

  streamDefVlist(streamID, vlistID);

  gridsize = gridInqSize(gridID);
  array = (double *) malloc(gridsize*sizeof(double));

  if ( operatorID == FOR )
    ntimesteps = 1.001 + ((rstop-rstart)/rinc);
  else
    {
      vlistDefNtsteps(vlistID, 0);
      ntimesteps = 1;
    }

  julday = date_to_julday(CALENDAR_PROLEPTIC, 10101);

  nvars = vlistNvars(vlistID);

  for ( tsID = 0; tsID < ntimesteps; tsID++ )
    {
      rval  = rstart + rinc*tsID;
      vdate = julday_to_date(CALENDAR_PROLEPTIC, julday + tsID);
      vtime = 0;
      taxisDefVdate(taxisID, vdate);
      taxisDefVtime(taxisID, vtime);
      streamDefTimestep(streamID, tsID);

      for ( varID = 0; varID < nvars; varID++ )
        {
          nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
          for ( levelID = 0; levelID < nlevels; levelID++ )
            {
              streamDefRecord(streamID, varID, levelID);

              if ( operatorID == RANDOM )
                {
                  for ( i = 0; i < gridsize; i++ )
                    array[i] = rand()/(RAND_MAX+1.0);
                }
              else if ( operatorID == SINCOS )
                {
		  int nlon = gridInqXsize(gridID);
		  int nlat = gridInqYsize(gridID);
		  double dlon = 360./nlon;
		  double dlat = 180./nlat;
		  double lon0 = 0;
		  double lat0 = -90 + dlat/2;

                  for ( i = 0; i < gridsize; i++ )
		    {
		      int ilat = (i%gridsize)/ nlon;
		      int ilon = i%nlon;
		      array[i] = cos(2.0 * M_PI * (lon0 + ilon*dlon)/360)
		   	       * sin(2.0 * M_PI * (lat0 + ilat*dlat)/180);
		    }
		}
              else if ( operatorID == CONST )
                {
                  for ( i = 0; i < gridsize; i++ )
                    array[i] = rconst;
                }
              else if ( operatorID == TOPO )
                {
#if defined(ENABLE_DATA)
                  for ( i = 0; i < gridsize; i++ )
                    array[i] = etopo[i]/etopo_scale - etopo_offset;
#else
                  cdoAbort("Operator support disabled!");
#endif
                }
              else if ( operatorID == TEMP )
                {
#if defined(ENABLE_DATA)
                  for ( i = 0; i < gridsize; i++ )
                    array[i] = temp[i]/temp_scale - temp_offset;
#else
                  cdoAbort("Operator support disabled!");
#endif
                }
              else if ( operatorID == MASK )
                {
#if defined(ENABLE_DATA)
                  for ( i = 0; i < gridsize; i++ )
                    array[i] = mask[i]/mask_scale - mask_offset;
#else
                  cdoAbort("Operator support disabled!");
#endif
                }
              else if ( operatorID == FOR )
                {
                  array[0] = rval;
                }
              else if ( operatorID == STDATM )
                {
                  array[0] = (varID == varID2) ? std_atm_temperatur(levels[levelID]) : std_atm_pressure(levels[levelID]);
                }

              streamWriteRecord(streamID, array, 0);
            }
        }
    }

  streamClose(streamID);

  vlistDestroy(vlistID);

  if ( array ) free(array);
  if ( levels ) free(levels); 

  cdoFinish();

  return (0);
}
