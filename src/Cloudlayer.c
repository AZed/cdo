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

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "vinterp.h"


#define  SCALESLP        (101325.0)


/* ================================================= */
/* LayerCloud calculates random overlap cloud cover */
/* ================================================= */

static
void layer_cloud(const double *cc, double *ll, long MaxLev, long MinLev, long dimgp)
{
  long i, k;
  double maxval, minval;
  double ZEPSEC;

  ZEPSEC = 1.-1.0e-12;

  for ( i = 0; i < dimgp; i++ ) ll[i] = 1. - cc[i+MaxLev*dimgp];

  //  printf("maxlev %d minlev %d\n", MaxLev, MinLev);

  for ( k = MaxLev + 1; k <= MinLev; k++ )
    {
      for ( i = 0; i < dimgp; i++ )
	{
	  maxval = MAX(cc[i+(k-1)*dimgp], cc[i+k*dimgp]);
	  minval = MIN(cc[i+(k-1)*dimgp], ZEPSEC);
	  ll[i] *= (1. - maxval) / (1. - minval);
	}
    }

  for ( i = 0; i < dimgp; i++ ) ll[i] = 1. - ll[i];
}

static
void vct2plev(const double *vct, double *plevs, long nlevs)
{
  long k;

  for ( k = 0; k < nlevs; k++ )
    plevs[k] = vct[k] + vct[k+nlevs] * SCALESLP;
  /*
  for ( k = 0; k < nlevs; k++ )
    printf("plevs %ld %g\n", k, plevs[k]);
  
  for ( k = 1; k < nlevs; k++ )
    printf("plevs %ld %g\n", k-1, (plevs[k]+plevs[k-1])/2);
  */
}

static
void hl_index(int *kmax, int *kmin, double pmax, double pmin, long nhlevs, double *pph)
{
  long k;
  long MaxLev, MinLev;
   
  for ( k = 0; k < nhlevs; k++ )
    if ( pph[k] > pmax ) break;
  
  MaxLev = k - 1;
  
  for ( k  = nhlevs - 1; k >= 0; k-- )
    if ( pph[k] < pmin ) break;
   
  MinLev = k;

  *kmax = MaxLev;
  *kmin = MinLev;
}

static
void pl_index(int *kmax, int *kmin, double pmax, double pmin, long nlevs, double *plevs)
{
  long k;
  long MaxLev = -1, MinLev = -1;
   
  for ( k = 0; k < nlevs; k++ )
    if ( plevs[k] >= pmax )
      {
	MaxLev = k;
	break;
      }
  
  for ( k  = nlevs - 1; k >= 0; k-- )
    if ( plevs[k] < pmin )
      {
	MinLev = k;
	break;
      }
   
  *kmax = MaxLev;
  *kmin = MinLev;
}


#define NVARS  3

void *Cloudlayer(void *argument)
{
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int gridID, zaxisID, tsID;
  int surfaceID;
  int nlevel, nhlev, nlevs, nrecs, recID, code;
  int varID, levelID;
  int zrev = FALSE;
  int nvars;
  int gridsize, i;
  int offset;
  int nmiss;
  int ngp = 0, ngrids;
  int aclcac_code;
  int aclcacID = -1;
  int nvars2 = 0;
  int aclcac_code_found = 0;
  int kmin[NVARS], kmax[NVARS];
  char varname[CDI_MAX_NAME];
  double sfclevel = 0;
  double *plevs = NULL;
  double *aclcac = NULL;
  double *cloud[NVARS];
  double missval;
  double pmin = 0, pmax = 0;

  cdoInitialize(argument);

  if ( operatorArgc() > 0 )
    {
      operatorCheckArgc(2);
      nvars2 = 1;
      pmin = atof(operatorArgv()[0]);
      pmax = atof(operatorArgv()[1]);
    }
  else
    {
      nvars2 = NVARS;
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  ngrids  = vlistNgrids(vlistID1);
  for ( i = 0; i < ngrids; i++ )
    {
      gridID = vlistGrid(vlistID1, i);
      if ( gridInqType(gridID) != GRID_SPECTRAL )
	{
	  ngp = gridInqSize(gridID);
	  break;
	}
    }

  /* check gridsize */
  for ( i = 0; i < ngrids; i++ )
    {
      gridID = vlistGrid(vlistID1, i);
      if ( gridInqType(gridID) != GRID_SPECTRAL )
	{
	  if ( ngp != gridInqSize(gridID) )
	    cdoAbort("Grids have different size!");
	}
    }

  aclcac_code = 223;

  nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; ++varID )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);

      code = vlistInqVarCode(vlistID1, varID);

      if ( code <= 0 )
	{
	  vlistInqVarName(vlistID1, varID, varname);
	  strtolower(varname);
	  if ( strcmp(varname, "aclcac") == 0 ) code = 223;
	}

      if  ( code == aclcac_code )
	{
	  aclcac_code_found = 1;
	  if ( zaxisInqType(zaxisID) == ZAXIS_PRESSURE || zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	    {
	      aclcacID  = varID;
	      break;
	    }
	}
    }

  if ( aclcacID == -1 )
    {
      if ( aclcac_code_found )
	cdoAbort("Cloud cover (parameter 223) not found on pressure or hybrid levels!");
      else
	cdoAbort("Cloud cover (parameter 223) not found!");
    }

  missval = vlistInqVarMissval(vlistID1, aclcacID);
  gridID  = vlistInqVarGrid(vlistID1, aclcacID);
  zaxisID = vlistInqVarZaxis(vlistID1, aclcacID);

  gridsize = gridInqSize(gridID);
  nlevel = zaxisInqSize(zaxisID);
  nhlev  = nlevel+1;

  aclcac = malloc(gridsize*nlevel*sizeof(double));
  for ( varID = 0; varID < nvars2; ++varID )
    cloud[varID] = malloc(gridsize*sizeof(double));

  if ( zaxisInqType(zaxisID) == ZAXIS_PRESSURE )
    {
      plevs = malloc(nlevel*sizeof(double));
      zaxisInqLevels(zaxisID, plevs);
      if ( plevs[0] > plevs[nlevel-1] )
	{
	  double ptmp;
	  zrev = TRUE;
	  for ( levelID = 0; levelID < nlevel/2; ++levelID )
	    {
	      ptmp = plevs[levelID];
	      plevs[levelID] = plevs[nlevel-1-levelID];
	      plevs[nlevel-1-levelID] = ptmp;
	    }
	}
      /*
      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  printf("level %d %g\n", levelID, plevs[levelID]);
	}
      */
      if ( nvars2 == 1 )
	{
	  pl_index(&kmax[0], &kmin[0], pmin, pmax, nlevel, plevs);
	}
      else
	{
	  pl_index(&kmax[2], &kmin[2],  5000., 44000., nlevel, plevs);
	  pl_index(&kmax[1], &kmin[1], 46000., 73000., nlevel, plevs);
	  pl_index(&kmax[0], &kmin[0], 75000.,101300., nlevel, plevs);
	}

      free(plevs);
    }
  else if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
    {
      int nvct;

      nvct = zaxisInqVctSize(zaxisID);
      if ( nlevel == (nvct/2 - 1) )
	{
	  double *vct;

	  vct = malloc(nvct*sizeof(double));
	  zaxisInqVct(zaxisID, vct);

	  nlevs = nlevel + 1;
	  plevs = malloc(nlevs*sizeof(double));
	  vct2plev(vct, plevs, nlevs);
	  free(vct);

	  if ( nvars2 == 1 )
	    {
	      hl_index(&kmax[0], &kmin[0], pmin, pmax, nhlev, plevs);
	    }
	  else
	    {
	      hl_index(&kmax[2], &kmin[2],  5000., 44000., nhlev, plevs);
	      hl_index(&kmax[1], &kmin[1], 46000., 73000., nhlev, plevs);
	      hl_index(&kmax[0], &kmin[0], 75000.,101300., nhlev, plevs);
	    }

	  free(plevs);
 	}
      else
	cdoAbort("Unsupported vertical coordinate table format!");
   }
  else
    cdoAbort("Unsupported Z-Axis type!");


  surfaceID = zaxisCreate(ZAXIS_SURFACE, 1);
  zaxisDefLevels(surfaceID, &sfclevel);

  vlistID2 = vlistCreate();

  if ( nvars2 == 1 )
    {
      varID = vlistDefVar(vlistID2, gridID, surfaceID, TSTEP_INSTANT);
      vlistDefVarCode(vlistID2, varID, 33);
      vlistDefVarName(vlistID2, varID, "cld_lay");
      vlistDefVarLongname(vlistID2, varID, "cloud layer");
      vlistDefVarMissval(vlistID2, varID, missval);
    }
  else
    {
      varID = vlistDefVar(vlistID2, gridID, surfaceID, TSTEP_INSTANT);
      vlistDefVarCode(vlistID2, varID, 34);
      vlistDefVarName(vlistID2, varID, "low_cld");
      vlistDefVarLongname(vlistID2, varID, "low cloud");
      vlistDefVarMissval(vlistID2, varID, missval);

      varID = vlistDefVar(vlistID2, gridID, surfaceID, TSTEP_INSTANT);
      vlistDefVarCode(vlistID2, varID, 35);
      vlistDefVarName(vlistID2, varID, "mid_cld");
      vlistDefVarLongname(vlistID2, varID, "mid cloud");
      vlistDefVarMissval(vlistID2, varID, missval);

      varID = vlistDefVar(vlistID2, gridID, surfaceID, TSTEP_INSTANT);
      vlistDefVarCode(vlistID2, varID, 36);
      vlistDefVarName(vlistID2, varID, "hih_cld");
      vlistDefVarLongname(vlistID2, varID, "high cloud");
      vlistDefVarMissval(vlistID2, varID, missval);
    }

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
     
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( zrev )
	    offset = (nlevel-1-levelID)*gridsize;
	  else
	    offset = levelID*gridsize;

	  if ( varID == aclcacID )
	    {
	      streamReadRecord(streamID1, aclcac+offset, &nmiss);
	      if ( nmiss != 0 ) cdoAbort("Missing values unsupported!");
	    }
	}

      for ( varID = 0; varID < nvars2; ++varID )
	{
	  for ( i = 0; i < gridsize; i++ ) cloud[varID][i] = missval;
	}

      for ( varID = 0; varID < nvars2; ++varID )
	{
	  if ( kmax[varID] != -1 && kmin[varID] != -1 )
	    layer_cloud(aclcac, cloud[varID], kmax[varID], kmin[varID], gridsize);
	}

      for ( varID = 0; varID < nvars2; ++varID )
	{
	  nmiss = 0;
	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(cloud[varID][i], missval) ) nmiss++;

	  streamDefRecord(streamID2, varID, 0);
	  streamWriteRecord(streamID2, cloud[varID], nmiss);
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);
 
  vlistDestroy(vlistID2);

  free(aclcac);
  for ( varID = 0; varID < nvars2; ++varID )
    free(cloud[varID]);

  cdoFinish();

  return (0);
}
