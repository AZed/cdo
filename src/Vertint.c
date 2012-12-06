/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Vertint    ml2pl           Model to pressure level interpolation
      Vertint    ml2hl           Model to height level interpolation
*/


#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "vinterp.h"
#include "list.h"

#define  C_EARTH_GRAV    (9.80665)

void *Vertint(void *argument)
{
  int ML2PL, ML2HL, ML2PLX, ML2HLX;
  int ML2PL_LP, ML2HL_LP, ML2PLX_LP, ML2HLX_LP;
  int operatorID;
  int mode;
  enum {ECHAM_MODE, WMO_MODE};
  enum {func_pl, func_hl};
  enum {type_lin, type_log};
  int geop_code = 0, temp_code = 0, ps_code = 0, lsp_code = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridsize, ngp = 0;
  int recID, nrecs;
  int i, k, offset;
  int tsID, varID, levelID;
  int nvars;
  int zaxisIDp, zaxisIDh = -1, nzaxis;
  int ngrids, gridID, zaxisID;
  int nplev, nhlev = 0, nhlevf = 0, nhlevh = 0, nlevel, maxlev;
  int *vert_index = NULL;
  int nvct;
  int geop_needed = FALSE;
  int geopID = -1, tempID = -1, psID = -1, lnpsID = -1, gheightID = -1;
  int code, param;
  int pnum, pcat, pdis;
  int **varnmiss = NULL, *pnmiss = NULL;
  int *varinterp = NULL;
  char paramstr[32];
  char varname[CDI_MAX_NAME], stdname[CDI_MAX_NAME];
  int *vars = NULL;
  double minval, maxval;
  double missval;
  double *plev = NULL, *phlev = NULL, *vct = NULL;
  double *rvct = NULL; /* reduced VCT for LM */
  double *single1, *single2;
  double **vardata1 = NULL, **vardata2 = NULL;
  double *geop = NULL, *ps_prog = NULL, *full_press = NULL, *half_press = NULL;
  double *hyb_press = NULL;
  int Extrapolate = 0;
  int taxisID1, taxisID2;
  int lhavevct;
  int mono_level;
  int instNum, tableNum;
  int useTable;
  int operfunc, opertype;
  LIST *flist = listNew(FLT_LIST);

  cdoInitialize(argument);

  ML2PL  = cdoOperatorAdd("ml2pl",  func_pl, type_lin, "pressure levels in pascal");
  ML2PLX = cdoOperatorAdd("ml2plx", func_pl, type_lin, "pressure levels in pascal");
  ML2HL  = cdoOperatorAdd("ml2hl",  func_hl, type_lin, "height levels in meter");
  ML2HLX = cdoOperatorAdd("ml2hlx", func_hl, type_lin, "height levels in meter");
  ML2PL_LP  = cdoOperatorAdd("ml2pl_lp",  func_pl, type_log, "pressure levels in pascal");
  ML2PLX_LP = cdoOperatorAdd("ml2plx_lp", func_pl, type_log, "pressure levels in pascal");
  ML2HL_LP  = cdoOperatorAdd("ml2hl_lp",  func_hl, type_log, "height levels in meter");
  ML2HLX_LP = cdoOperatorAdd("ml2hlx_lp", func_hl, type_log, "height levels in meter");

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);
  opertype = cdoOperatorF2(operatorID);

  if ( operatorID == ML2PL || operatorID == ML2HL || operatorID == ML2PL_LP || operatorID == ML2HL_LP )
    {
      char *envstr;
      envstr = getenv("EXTRAPOLATE");

      if ( envstr )
	{
	  if ( isdigit((int) envstr[0]) )
	    {
	      Extrapolate = atoi(envstr);
	      if ( Extrapolate == 1 )
		cdoPrint("Extrapolation of missing values enabled!");
	    }
	}
    }
  else
    {
      Extrapolate = 1;
    }

  operatorInputArg(cdoOperatorEnter(operatorID));

  nplev = args2fltlist(operatorArgc(), operatorArgv(), flist);
  plev  = (double *) listArrayPtr(flist);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

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

  if ( operfunc == func_hl )
    zaxisIDp = zaxisCreate(ZAXIS_HEIGHT, nplev);
  else
    zaxisIDp = zaxisCreate(ZAXIS_PRESSURE, nplev);

  zaxisDefLevels(zaxisIDp, plev);
  nzaxis  = vlistNzaxis(vlistID1);
  lhavevct = FALSE;
  for ( i = 0; i < nzaxis; i++ )
    {
      /* mono_level = FALSE; */
      mono_level = TRUE;
      zaxisID = vlistZaxis(vlistID1, i);
      nlevel  = zaxisInqSize(zaxisID);

      if ( (zaxisInqType(zaxisID) == ZAXIS_HYBRID || zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF) &&
	   nlevel > 1 )
	{
	  double *level;
	  int l;
	  level = (double *) malloc(nlevel*sizeof(double));
	  zaxisInqLevels(zaxisID, level);
	  for ( l = 0; l < nlevel; l++ )
	    {
	      if ( (l+1) != (int) (level[l]+0.5) ) break;
	    }
	  if ( l == nlevel ) mono_level = TRUE; 
	  free(level);
	}

      if ( (zaxisInqType(zaxisID) == ZAXIS_HYBRID || zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF) &&
	   nlevel > 1 && mono_level )
	{
	  nvct = zaxisInqVctSize(zaxisID);
	  if ( nlevel == (nvct/2 - 1) )
	    {
	      if ( lhavevct == FALSE )
		{
		  lhavevct = TRUE;
		  zaxisIDh = zaxisID;
		  nhlev    = nlevel;
		  nhlevf   = nhlev;
		  nhlevh   = nhlevf + 1;
	      
		  vct = (double *) malloc(nvct*sizeof(double));
		  zaxisInqVct(zaxisID, vct);

		  vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	      else
		{
		  if ( memcmp(vct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double)) == 0 )
		    vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	    }
	  else if ( nlevel == (nvct/2) )
	    {
	      if ( lhavevct == FALSE )
		{
		  lhavevct = TRUE;
		  zaxisIDh = zaxisID;
		  nhlev    = nlevel;
		  nhlevf   = nhlev - 1;
		  nhlevh   = nhlev;
	      
		  vct = (double *) malloc(nvct*sizeof(double));
		  zaxisInqVct(zaxisID, vct);

		  vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	      else
		{
		  if ( memcmp(vct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double)) == 0 )
		    vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	    }
	  else if ( nlevel == (nvct - 4 - 1) )
	    {
	      if ( lhavevct == FALSE )
		{
		  int vctsize;
		  int voff = 4;
		  
		  rvct = (double *) malloc(nvct*sizeof(double));
		  zaxisInqVct(zaxisID, rvct);

		  if ( (int)(rvct[0]+0.5) == 100000 && rvct[voff] < rvct[voff+1] )
		    {
		      lhavevct = TRUE;
		      zaxisIDh = zaxisID;
		      nhlev    = nlevel;
		      nhlevf   = nhlev;
		      nhlevh   = nhlev + 1;

		      vctsize = 2*nhlevh;
		      vct = (double *) malloc(vctsize*sizeof(double));

		      vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);

		      /* calculate VCT for LM */

		      for ( i = 0; i < vctsize/2; i++ )
			{
			  if ( rvct[voff+i] >= rvct[voff] && rvct[voff+i] <= rvct[3] )
			    {
			      vct[i] = rvct[0]*rvct[voff+i];
			      vct[vctsize/2+i] = 0;
			    }
			  else
			    {
			      vct[i] = (rvct[0]*rvct[3]*(1-rvct[voff+i]))/(1-rvct[3]);
			      vct[vctsize/2+i] = (rvct[voff+i]-rvct[3])/(1-rvct[3]);
			    }
			}
		      
		      if ( cdoVerbose )
			{
			  for ( i = 0; i < vctsize/2; i++ )
			    fprintf(stdout, "%5d %25.17f %25.17f\n", i, vct[i], vct[vctsize/2+i]);
			}
		    }
		}
	      else
		{
		  if ( memcmp(rvct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double)) == 0 )
		    vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	    }
	}
    }

  nvars = vlistNvars(vlistID1);

  vars      = (int *) malloc(nvars*sizeof(int));
  vardata1  = (double **) malloc(nvars*sizeof(double*));
  vardata2  = (double **) malloc(nvars*sizeof(double*));
  varnmiss  = (int **) malloc(nvars*sizeof(int*));
  varinterp = (int *) malloc(nvars*sizeof(int));

  maxlev   = nhlevh > nplev ? nhlevh : nplev;

  if ( Extrapolate == 0 )
    pnmiss   = (int *) malloc(nplev*sizeof(int));

  if ( zaxisIDh != -1 && ngp > 0 )
    {
      vert_index = (int *) malloc(ngp*nplev*sizeof(int));
      ps_prog    = (double *) malloc(ngp*sizeof(double));
      full_press = (double *) malloc(ngp*nhlevf*sizeof(double));
      half_press = (double *) malloc(ngp*nhlevh*sizeof(double));
    }
  else
    cdoWarning("No data on hybrid model level found!");

  if ( operfunc == func_hl )
    {
      phlev = (double *) malloc(nplev*sizeof(double));
      h2p(phlev, plev, nplev);
      memcpy(plev, phlev, nplev*sizeof(double));
      free(phlev);
    }

  if ( opertype == type_log )
    for ( k = 0; k < nplev; k++ ) plev[k] = log(plev[k]);

  useTable = FALSE;
  for ( varID = 0; varID < nvars; varID++ )
    {
      tableNum = tableInqNum(vlistInqVarTable(vlistID1, varID));
      if ( tableNum > 0 && tableNum != 255 )
	{
	  useTable = TRUE;
	  break;
	}
    }

  if ( cdoVerbose && useTable ) cdoPrint("Using code tables!");

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);
      instNum  = institutInqCenter(vlistInqVarInstitut(vlistID1, varID));
      tableNum = tableInqNum(vlistInqVarTable(vlistID1, varID));

      code     = vlistInqVarCode(vlistID1, varID);
      param    = vlistInqVarParam(vlistID1, varID);

      cdiParamToString(param, paramstr, sizeof(paramstr));
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      if ( pdis >= 0 && pdis < 255 ) code = -1;

      if ( useTable )
	{
	  if ( tableNum == 2 )
	    {
	      mode = WMO_MODE;
	      geop_code  =   6;
	      temp_code  =  11;
	      ps_code    =   1;
	    }
	  else if ( tableNum == 128 || tableNum == 0 )
	    {
	      mode = ECHAM_MODE;
	      geop_code  = 129;
	      temp_code  = 130;
	      ps_code    = 134;
	      lsp_code   = 152;
	    }
	  else
	    mode = -1;
	}
      else
	{
	  mode = ECHAM_MODE;
	  geop_code  = 129;
	  temp_code  = 130;
	  ps_code    = 134;
	  lsp_code   = 152;
	}

      if ( cdoVerbose )
	cdoPrint("Mode = %d  Center = %d  Code = %d  Param = %s", mode, instNum, code, paramstr);

      if ( code <= 0 )
	{
	  vlistInqVarName(vlistID1, varID, varname);
	  strtolower(varname);

	  vlistInqVarStdname(vlistID1, varID, stdname);
	  strtolower(stdname);

	  if ( strcmp(stdname, "surface_air_pressure") == 0 ) code = 134;
	  else
	    {
	      /*                        ECHAM                            ECMWF       */
	      if      ( strcmp(varname, "geosp") == 0 || strcmp(varname, "z")    == 0 ) code = 129;
	      else if ( strcmp(varname, "st")    == 0 || strcmp(varname, "t")    == 0 ) code = 130;
	      else if ( strcmp(varname, "aps")   == 0 || strcmp(varname, "sp"  ) == 0 ) code = 134;
	      else if ( strcmp(varname, "lsp")   == 0 || strcmp(varname, "lnsp") == 0 ) code = 152;
	      /* else if ( strcmp(varname, "geopoth") == 0 ) code = 156; */
	    }
	}

      if ( mode == ECHAM_MODE )
	{
	  if      ( code == geop_code  && nlevel == 1      ) geopID    = varID;
	  else if ( code == temp_code  && nlevel == nhlevf ) tempID    = varID;
	  else if ( code == ps_code    && nlevel == 1      ) psID      = varID;
	  else if ( code == lsp_code   && nlevel == 1      ) lnpsID    = varID;
	  else if ( code == 156        && nlevel == nhlevf ) gheightID = varID;
	}
      else if ( mode == WMO_MODE )
	{
	  if      ( code == geop_code  && nlevel == 1      ) geopID  = varID;
	  else if ( code == temp_code  && nlevel == nhlevf ) tempID  = varID;
	  else if ( code == ps_code    && nlevel == 1      ) psID    = varID;
	}

      if ( gridInqType(gridID) == GRID_SPECTRAL && zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	cdoAbort("Spectral data on model level unsupported!");

      if ( gridInqType(gridID) == GRID_SPECTRAL )
	cdoAbort("Spectral data unsupported!");

      if ( varID == gheightID )
	vardata1[varID] = (double *) malloc(gridsize*(nlevel+1)*sizeof(double));
      else
	vardata1[varID] = (double *) malloc(gridsize*nlevel*sizeof(double));

      /* if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID && zaxisIDh != -1 && nlevel == nhlev ) */
      if ( zaxisID == zaxisIDh ||
	   (zaxisInqType(zaxisID) == ZAXIS_HYBRID && zaxisIDh != -1 && (nlevel == nhlevh || nlevel == nhlevf)) )
	{
	  varinterp[varID] = TRUE;
	  vardata2[varID]  = (double *) malloc(gridsize*nplev*sizeof(double));
	  varnmiss[varID]  = (int *) malloc(maxlev*sizeof(int));
	  memset(varnmiss[varID], 0, maxlev*sizeof(int));
	}
      else
	{
	  if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID && zaxisIDh != -1 && nlevel > 1 )
	    cdoWarning("Parameter %d has wrong number of levels, skipped! (param=%s nlevel=%d)",
		       varID+1, paramstr, nlevel);
	  varinterp[varID] = FALSE;
	  vardata2[varID]  = vardata1[varID];
	  varnmiss[varID]  = (int *) malloc(nlevel*sizeof(int));
	}
    }

  if ( tempID != -1 || gheightID != -1 ) geop_needed = TRUE;

  if ( zaxisIDh != -1 && geop_needed )
    {
      geop = (double *) malloc(ngp*sizeof(double));
      if ( geopID == -1 )
	{
	  cdoWarning("Orography (surf. geopotential) not found - using zero orography!");
	  memset(geop, 0, ngp*sizeof(double));
	}
    }

  if ( zaxisIDh != -1 && gheightID != -1 && tempID == -1 )
    cdoAbort("Temperature not found, needed to compute geopotheight!");

  if ( zaxisIDh != -1 && lnpsID == -1 )
    {
      if ( psID != -1 )
	{
	  param = vlistInqVarParam(vlistID1, psID);
	  cdiParamToString(param, paramstr, sizeof(paramstr));
	  if ( cdoVerbose )
	    cdoPrint("LOG surface pressure not found - using surface pressure (param=%s)!", paramstr);
	}
      else
	cdoAbort("Surface pressure not found!");
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      for ( varID = 0; varID < nvars; ++varID ) vars[varID] = FALSE;

      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single1  = vardata1[varID] + offset;
	  
	  streamReadRecord(streamID1, single1, &varnmiss[varID][levelID]);
	  vars[varID] = TRUE;
	}

      if ( zaxisIDh != -1 )
	{
	  if ( geop_needed && geopID != -1 )
	    {
	      memcpy(geop, vardata1[geopID], ngp*sizeof(double));

	      /* check range of geop */
	      minmaxval(ngp, geop, NULL, &minval, &maxval);
	      if ( minval < MIN_FIS || maxval > MAX_FIS )
		cdoWarning("Surface geopotential out of range (min=%g max=%g)!", minval, maxval);
	      if ( minval >= 0 && maxval <= 1000 )
		cdoWarning("Surface geopotential has an unexpected range (min=%g max=%g)!", minval, maxval);
	    }

	  if ( lnpsID != -1 )
	    for ( i = 0; i < ngp; i++ ) ps_prog[i] = exp(vardata1[lnpsID][i]);
	  else if ( psID != -1 )
	    memcpy(ps_prog, vardata1[psID], ngp*sizeof(double));

	  /* check range of ps_prog */
	  minmaxval(ngp, ps_prog, NULL, &minval, &maxval);
	  if ( minval < MIN_PS || maxval > MAX_PS )
	    cdoWarning("Surface pressure out of range (min=%g max=%g)!", minval, maxval);


	  presh(full_press, half_press, vct, ps_prog, nhlevf, ngp);

	  if ( opertype == type_log )
	    {
	      for ( i = 0; i < ngp; i++ ) ps_prog[i] = log(ps_prog[i]);

	      for ( k = 0; k < nhlevh; k++ )
		for ( i = 0; i < ngp; i++ )
		  half_press[k*ngp+i] = log(half_press[k*ngp+i]);

	      for ( k = 0; k < nhlevf; k++ )
		for ( i = 0; i < ngp; i++ )
		  full_press[k*ngp+i] = log(full_press[k*ngp+i]);
	    }

	  genind(vert_index, plev, full_press, ngp, nplev, nhlevf);

	  if ( Extrapolate == 0 )
	    genindmiss(vert_index, plev, ngp, nplev, ps_prog, pnmiss);
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	      missval  = vlistInqVarMissval(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevel   = zaxisInqSize(zaxisID);
	      if ( varinterp[varID] )
		{
		  /*
		  if ( nlevel == nhlevh )
		    {
		      int i, k;
		      double *vl1, *vl2;

		      for ( k = 1; k < nlevel; k++ )
			{
			  vl1  = vardata1[varID] + gridsize*(k-1);
			  vl2  = vardata1[varID] + gridsize*(k);
			  for ( i = 0; i < gridsize; i++ )
			    vl1[i] = 0.5*(vl1[i] + vl2[i]);
			}
		      
		      nlevel = nhlevf;
		    }
		  */
		  if ( nlevel == nhlevh )
		    {
		      hyb_press = half_press;
		    }
		  else if ( nlevel == nhlevf )
		    {
		      hyb_press = full_press;
		    }
		  else
		    {
		      param = vlistInqVarParam(vlistID1, varID);
		      cdiParamToString(param, paramstr, sizeof(paramstr));
		      cdoAbort("Number of hybrid level differ from full/half level (param=%s)!", paramstr);
		    }

		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    {
		      if ( varnmiss[varID][levelID] )
			cdoAbort("Missing values unsupported for this operator!");
		    }

		  if ( varID == tempID )
		    {
		      if ( nlevel == nhlevh )
			cdoAbort("Temperature on half level unsupported!");

		      if ( opertype == type_log && Extrapolate )
			cdoAbort("Log. extrapolation of temperature unsupported!");

		      interp_T(geop, vardata1[varID], vardata2[varID],
			       full_press, half_press, vert_index,
			       plev, nplev, ngp, nlevel, missval);
		    }
		  else if ( varID == gheightID )
		    {
		      for ( i = 0; i < ngp; ++i )
			vardata1[varID][ngp*nlevel+i] = geop[i]/C_EARTH_GRAV;

		      interp_Z(geop, vardata1[varID], vardata2[varID],
			       full_press, half_press, vert_index, vardata1[tempID],
			       plev, nplev, ngp, nlevel, missval);
		    }
		  else
		    {
		      interp_X(vardata1[varID], vardata2[varID], hyb_press,
			       vert_index, plev, nplev, ngp, nlevel, missval);
		    }
		  
		  if ( Extrapolate == 0 )
		    memcpy(varnmiss[varID], pnmiss, nplev*sizeof(int));
		}
	    }
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
		  offset   = gridsize*levelID;
		  single2  = vardata2[varID] + offset;
		  streamDefRecord(streamID2, varID, levelID);
		  streamWriteRecord(streamID2, single2, varnmiss[varID][levelID]);
		}
	    }
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      free(varnmiss[varID]);
      free(vardata1[varID]);
      if ( varinterp[varID] ) free(vardata2[varID]);
    }

  free(varinterp);
  free(varnmiss);
  free(vardata2);
  free(vardata1);
  free(vars);

  if ( pnmiss     ) free(pnmiss);

  if ( geop       ) free(geop);
  if ( ps_prog    ) free(ps_prog);
  if ( vert_index ) free(vert_index);
  if ( full_press ) free(full_press);
  if ( half_press ) free(half_press);
  if ( vct        ) free(vct);
  if ( rvct       ) free(rvct);

  listDelete(flist);

  cdoFinish();

  return (0);
}
