/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2007-2011 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Derivepar     geopotheight          geopotential height
*/

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "vinterp.h"

#define  C_RKBOL         (1.380658e-23)     /* Boltzmann constant in J/K   */
#define  C_RNAVO         (6.0221367e+23)    /* Avogadro constant in 1/mol  */
#define  C_RMD           (28.9644)          /* molecular weight of dry air */
#define  C_RMV           (18.0153)          /* molecular weight of water vapor */
#define  C_R             (C_RKBOL * C_RNAVO)
#define  C_RV            (1000. * C_R / C_RMV)

#define  C_EARTH_GRAV    (9.80665)
#define  C_RKBOL         (1.380658e-23)     /* Boltzmann constant in J/K   */
#define  C_RNAVO         (6.0221367e+23)    /* Avogadro constant in 1/mol  */
#define  C_RMD           (28.9644)          /* molecular weight of dry air */
#define  C_R             (C_RKBOL * C_RNAVO)
#define  C_EARTH_RD      (1000. * C_R / C_RMD)

static double Grav          = C_EARTH_GRAV;
static double RD            = C_EARTH_RD;

void MakeGeopotHeight(double *geop, double* gt, double *gq, double *ph, double *pcw, double * pci, int nhor, int nlev)
{
  int i, j;
  double vtmp;
  double zrg;
  double z2log2;
  double *geopl, *gtl, *gql, *phl;

  z2log2 = 2.0 * log(2.0);
  vtmp   = (C_RV / RD) - 1.0;
  zrg    = 1.0 / Grav;

  if ( gq ) /* Humidity is present */
    {
      for ( j = nlev ; j > 1 ; j-- )
        {
          geopl = geop + nhor*(j-1);
          gtl   = gt   + nhor*(j-1);
          gql   = gq   + nhor*(j-1);
          phl   = ph   + nhor*(j-1);
#if defined (SX)
#pragma vdir nodep
#endif
#if defined (_OPENMP)
#pragma omp parallel for
#endif
          for ( i = 0; i < nhor; i++ )
            geopl[i] = geopl[i+nhor] + RD * gtl[i] * (1.0 + vtmp * gql[i])
                     * log(phl[i+nhor] / phl[i]);
        }

#if defined (SX)
#pragma vdir nodep
#endif
#if defined (_OPENMP)
#pragma omp parallel for
#endif
      for ( i = 0; i < nhor; i++ )
        geop[i] = geop[i+nhor] + RD * gt[i] * (1.0 + vtmp * gq[i]) * z2log2;
    }
  else    /* No humidity */
    {
      for ( j = nlev ; j > 1 ; j-- )
#if defined (SX)
#pragma vdir nodep
#endif
        for ( i = nhor * (j-1) ; i < nhor * j ; i++ )
          geop[i] = geop[i+nhor] + RD * gt[i] * log(ph[i+nhor] / ph[i]);

#if defined (SX)
#pragma vdir nodep
#endif
      for ( i = 0; i < nhor; i++ )
        geop[i] = geop[i+nhor] + RD * gt[i] * z2log2;
    }

#if defined (SX)
#pragma vdir nodep
#endif
#if defined (_OPENMP)
#pragma omp parallel for
#endif
  for ( i = 0; i < nhor * (nlev+1); i++ ) geop[i] *= zrg;
}

static
void minmax(int nvals, double *array, int *imiss, double *minval, double *maxval)
{
  long i;
  double xmin =  DBL_MAX;
  double xmax = -DBL_MAX;

  if ( imiss )
    {
      for ( i = 0; i < nvals; i++ )
	{
	  if ( ! imiss[i] )
	    {
	      if      ( array[i] > xmax ) xmax = array[i];
	      else if ( array[i] < xmin ) xmin = array[i];
	    }
	}
    }
  else
    {
      for ( i = 0; i < nvals; i++ )
	{
	  if      ( array[i] > xmax ) xmax = array[i];
	  else if ( array[i] < xmin ) xmin = array[i];
	}
    }

  *minval = xmin;
  *maxval = xmax;
}


void *Derivepar(void *argument)
{
  int GEOPOTHEIGHT;
  int operatorID;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridsize, ngp = 0;
  int recID, nrecs;
  int i, offset, iv;
  int tsID, varID, levelID;
  int nvars;
  int zaxisID2, zaxisIDh = -1, nzaxis, surfaceID;
  int ngrids, gridID, zaxisID;
  int nlevel;
  int nvct;
  int geopID = -1, tempID = -1, humID = -1, psID = -1, lnpsID = -1, presID = -1, clwcID = -1, ciwcID = -1;
  int code;
  char varname[CDI_MAX_NAME];
  double *single2;
  int taxisID1, taxisID2;
  int lhavevct;
  int nhlevf = 0;
  double *lev2;
  double *vct = NULL;
  double *geop = NULL, *ps = NULL, *temp = NULL, *hum = NULL, *lwater = NULL, *iwater = NULL;
  double *geopotheight = NULL;
  int nmiss, nmissout = 0;
  int ltq = FALSE;
  int *imiss = NULL;
  double *array = NULL;
  double *half_press = NULL;
  double minval, maxval;
  double missval = 0;
  double ps_min =  20000, ps_max = 120000;
  double fis_min = -100000, fis_max = 100000;
  double t_min = 150, t_max = 400;
  double q_min = 0, q_max = 0.1;
  double cconst = 1.E-6;
  const char *fname;


  cdoInitialize(argument);

  GEOPOTHEIGHT = cdoOperatorAdd("geopotheight",   0, 0, NULL);

  operatorID = cdoOperatorID();


  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  ngrids  = vlistNgrids(vlistID1);
  for ( i = 0; i < ngrids; i++ )
    {
      gridID = vlistGrid(vlistID1, i);
      if ( gridInqType(gridID) == GRID_SPECTRAL )
	{
	  cdoAbort("Spectral data unsupported!");
	}
      else
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


  nzaxis  = vlistNzaxis(vlistID1);
  lhavevct = FALSE;

  if ( cdoVerbose )
    cdoPrint("nzaxis: %d", nzaxis);

  for ( i = 0; i < nzaxis; i++ )
    {
      zaxisID = vlistZaxis(vlistID1, i);
      nlevel  = zaxisInqSize(zaxisID);
      if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	{
	  if ( nlevel > 1 )
	    {
	      nvct = zaxisInqVctSize(zaxisID);

              if ( cdoVerbose )
                cdoPrint("i: %d, vct size of zaxisID %d = %d", i, zaxisID, nvct);

	      if ( nlevel == (nvct/2 - 1) )
		{
		  if ( lhavevct == FALSE )
		    {
		      lhavevct = TRUE;
		      zaxisIDh = zaxisID;
		      nhlevf   = nlevel;
	      
                      if ( cdoVerbose )
                        cdoPrint("lhavevct=TRUE  zaxisIDh = %d, nhlevf   = %d", zaxisIDh, nlevel);
 
		      vct = (double *) malloc(nvct*sizeof(double));
		      zaxisInqVct(zaxisID, vct);

		      if ( cdoVerbose )
			for ( i = 0; i < nvct/2; ++i )
			  cdoPrint("vct: %5d %25.17f %25.17f", i, vct[i], vct[nvct/2+i]);
		    }
		}
              else 
                {
		  if ( cdoVerbose )
		    cdoPrint("nlevel /= (nvct/2 - 1): nlevel = %d", nlevel);
                }
	    }
	}
    }

  if ( zaxisIDh == -1 )
    cdoAbort("No data on hybrid model level found!");

  nvars = vlistNvars(vlistID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID  = vlistInqVarGrid(vlistID1, varID);
      zaxisID = vlistInqVarZaxis(vlistID1, varID);
      nlevel  = zaxisInqSize(zaxisID);

      code = vlistInqVarCode(vlistID1, varID);
      /* code = -1; */
      if ( code <= 0 )
	{
	  vlistInqVarName(vlistID1, varID, varname);

	  strtolower(varname);

	  if ( nlevel == 1 )
	    {
	      if      ( strcmp(varname, "geosp")   == 0 ) code = 129;
	      else if ( strcmp(varname, "aps")     == 0 ) code = 134;
	      else if ( strcmp(varname, "ps")      == 0 ) code = 134;
	      else if ( strcmp(varname, "lsp")     == 0 ) code = 152;
	    }

	  if ( nlevel == nhlevf )
	    {
	      if      ( strcmp(varname, "t")       == 0 ) code = 130;
	      else if ( strcmp(varname, "q")       == 0 ) code = 133;
	      else if ( strcmp(varname, "clwc")    == 0 ) code = 246;
	      else if ( strcmp(varname, "ciwc")    == 0 ) code = 247;
	    }
	}

      if      ( code == 129 ) geopID    = varID;
      else if ( code == 130 ) tempID    = varID;
      else if ( code == 133 ) humID     = varID;
      else if ( code == 134 ) psID      = varID;
      else if ( code == 152 ) lnpsID    = varID;
      else if ( code == 246 ) clwcID    = varID;
      else if ( code == 247 ) ciwcID    = varID;

      if ( gridInqType(gridID) == GRID_SPECTRAL && zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	cdoAbort("Spectral data on model level unsupported!");

      if ( gridInqType(gridID) == GRID_SPECTRAL )
	cdoAbort("Spectral data unsupported!");


      if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID && zaxisIDh != -1 && nlevel == nhlevf )
	{
	}
      else
	{
	  if ( code == 130 ) tempID = -1;
	  if ( code == 133 ) humID  = -1;
	  if ( code == 246 ) clwcID  = -1;
	  if ( code == 247 ) ciwcID  = -1;
	}
    }

  if ( tempID == -1 ) cdoAbort("Temperature not found!");

  array  = (double *) malloc(ngp*sizeof(double));

  geop   = (double *) malloc(ngp*sizeof(double));
  ps     = (double *) malloc(ngp*sizeof(double));

  temp   = (double *) malloc(ngp*nhlevf*sizeof(double));
  hum    = (double *) malloc(ngp*nhlevf*sizeof(double));
  lwater = (double *) malloc(ngp*nhlevf*sizeof(double));
  iwater = (double *) malloc(ngp*nhlevf*sizeof(double));

  half_press   = (double *) malloc(ngp*(nhlevf+1)*sizeof(double));
  geopotheight = (double *) malloc(ngp*(nhlevf+1)*sizeof(double));

  if ( zaxisIDh != -1 && geopID == -1 )
    {
      if ( ltq )
	cdoWarning("Orography (surf. geopotential) not found - using zero orography!");

      memset(geop, 0, ngp*sizeof(double));
    }

  presID = lnpsID;
  if ( zaxisIDh != -1 && lnpsID == -1 )
    {
      presID = psID;
      if ( psID != -1 )
	cdoWarning("LOG surface pressure (lsp) not found - using surface pressure (asp)!");
      else
	cdoAbort("Surface pressure not found!");
    }


  vlistID2 = vlistCreate();
  varID = vlistDefVar(vlistID2, gridID, zaxisIDh, TIME_VARIABLE);
  vlistDefVarParam(vlistID2, varID, cdiEncodeParam(156, 128, 255));
  vlistDefVarName(vlistID2, varID, "geopotheight");
  vlistDefVarStdname(vlistID2, varID, "geopotental_height");
  vlistDefVarUnits(vlistID2, varID, "m");

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
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	  nlevel   = zaxisInqSize(zaxisID);
	  offset   = gridsize*levelID;
	  streamReadRecord(streamID1, array, &nmiss);

	  if ( zaxisIDh != -1 )
	    {
	      if ( varID == geopID )
		{
		  memcpy(geop, array, ngp*sizeof(double));
		}
	      else if ( varID == presID )
		{
		  if ( lnpsID != -1 )
		    for ( i = 0; i < ngp; ++i ) ps[i] = exp(array[i]);
		  else if ( psID != -1 )
		    memcpy(ps, array, ngp*sizeof(double));
		}
	      else if ( varID == tempID )
		memcpy(temp+offset, array, ngp*sizeof(double));
	      else if ( varID == humID )
		memcpy(hum+offset, array, ngp*sizeof(double));
	      else if ( varID == clwcID )
		memcpy(lwater+offset, array, ngp*sizeof(double));
	      else if ( varID == ciwcID )
		memcpy(iwater+offset, array, ngp*sizeof(double));
	    }
	}

      if ( zaxisIDh != -1 )
	{
	  /* check range of ps_prog */

	  minmax(ngp, ps, imiss, &minval, &maxval);

	  if ( minval < ps_min || maxval > ps_max )
	    cdoWarning("Surface pressure out of range (min=%g max=%g)!", minval, maxval);

	  /* check range of geop */

	  minmax(ngp, geop, imiss, &minval, &maxval);

	  if ( minval < fis_min || maxval > fis_max )
	    cdoWarning("Orography out of range (min=%g max=%g)!", minval, maxval);
	}

      varID = tempID;
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single2  = temp + offset;

	  minmax(ngp, single2, imiss, &minval, &maxval);
	  if ( minval < t_min || maxval > t_max )
	    cdoWarning("Input temperature at level %d out of range (min=%g max=%g)!",
		       levelID+1, minval, maxval);
	}

      /*
      if ( humID != -1 )
	{
	  varID = humID;
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	      offset   = gridsize*levelID;
	      single2  = hum + offset;

	      corr_hum(gridsize, single2, q_min);

	      minmax(ngp, single2, imiss, &minval, &maxval);
	      if ( minval < q_min || maxval > q_max )
		cdoWarning("Input humidity at level %d out of range (min=%g max=%g)!",
			   levelID+1, minval, maxval);
	    }
	}
      */
      presh(NULL, half_press, vct, ps, nhlevf, ngp);

      memcpy(geopotheight+ngp*nhlevf, geop, ngp*sizeof(double));
      MakeGeopotHeight(geopotheight, temp, hum, half_press,lwater,iwater, ngp, nhlevf);

      nmissout = 0;
      varID = 0;
      nlevel = nhlevf;
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, geopotheight+levelID*ngp, nmissout);
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( imiss ) free(imiss);

  free(ps);
  free(geop);
  free(temp);
  free(geopotheight);
  if ( hum ) free(hum);

  if ( half_press ) free(half_press);

  free(array);
  if ( vct ) free(vct);

  cdoFinish();

  return (0);
}
