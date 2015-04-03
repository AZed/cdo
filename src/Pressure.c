/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2009 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Pressure    fpressure          Pressure on full hybrid level
      Pressure    hpressure          Pressure on half hybrid level
*/


#include <ctype.h>
#include <string.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "vinterp.h"
#include "list.h"


void *Pressure(void *argument)
{
  static char func[] = "Pressure";
  int FPRESSURE, HPRESSURE;
  int operatorID;
  int mode;
  enum {ECHAM_MODE, WMO_MODE};
  int geop_code = 0, temp_code = 0, ps_code = 0, lsp_code = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridsize, ngp = 0;
  int recID, nrecs;
  int i, offset;
  int tsID, varID, levelID;
  int nvars;
  int zaxisIDp, zaxisIDh = -1, nzaxis;
  int ngrids, gridID, zaxisID;
  int nhlev = 0, nhlevf = 0, nhlevh = 0, nlevel;
  int nvct;
  int geopID = -1, tempID = -1, psID = -1, lnpsID = -1, pvarID;
  int code;
  char varname[128];
  double *vct = NULL;
  double *rvct = NULL; /* reduced VCT for LM */
  double *ps_prog = NULL, *full_press = NULL, *half_press = NULL;
  double *hyb_press = NULL;
  double *pdata = NULL;
  int taxisID1, taxisID2;
  int lhavevct;
  int nmiss;
  int mono_level;
  int instNum, tableNum;
  int useTable;
  LIST *flist = listNew(FLT_LIST);

  cdoInitialize(argument);

  FPRESSURE = cdoOperatorAdd("fpressure", 0, 0, NULL);
  HPRESSURE = cdoOperatorAdd("hpressure", 0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

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

  nzaxis  = vlistNzaxis(vlistID1);
  lhavevct = FALSE;
  for ( i = 0; i < nzaxis; i++ )
    {
      mono_level = FALSE;
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
		  memcpy(vct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double));
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
		  memcpy(vct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double));
		}
	    }
	  else if ( nlevel == (nvct - 4 - 1) )
	    {
	      if ( lhavevct == FALSE )
		{
		  int vctsize;
		  int voff = 4;
		  const double *pvct = zaxisInqVctPtr(zaxisID);

		  if ( (int)(pvct[0]+0.5) == 100000 && pvct[voff] < pvct[voff+1] )
		    {
		      lhavevct = TRUE;
		      zaxisIDh = zaxisID;
		      nhlev    = nlevel;
		      nhlevf   = nhlev;
		      nhlevh   = nhlev + 1;

		      vctsize = 2*nhlevh;
		      vct = (double *) malloc(vctsize*sizeof(double));
		      rvct = (double *) malloc(nvct*sizeof(double));
		      memcpy(rvct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double));

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
	    }
	}
    }


  nvars = vlistNvars(vlistID1);

  if ( zaxisIDh != -1 && ngp > 0 )
    {
      ps_prog    = (double *) malloc(ngp*sizeof(double));
      full_press = (double *) malloc(ngp*nhlevf*sizeof(double));
      half_press = (double *) malloc(ngp*nhlevh*sizeof(double));
    }
  else
    cdoAbort("No data on hybrid model level found!");

  if ( operatorID == FPRESSURE )
    zaxisIDp = zaxisCreate(ZAXIS_HYBRID, nhlevf);
  else
    zaxisIDp = zaxisCreate(ZAXIS_HYBRID_HALF, nhlevh);

  {
    double *level;
    int l;
    level = (double *) malloc(nhlevh*sizeof(double));
    for ( l = 0; l < nhlevh; l++ ) level[l] = l+1;
    zaxisDefLevels(zaxisIDp, level);
    free(level);
  }

  zaxisDefVct(zaxisIDp, 2*nhlevh, vct);

  useTable = FALSE;
  for ( varID = 0; varID < nvars; varID++ )
    {
      tableNum = tableInqNum(vlistInqVarTable(vlistID1, varID));

      if ( tableNum > 0 )
	{
	  useTable = TRUE;
	}
    }

  if ( cdoVerbose && useTable ) cdoPrint("Use code tables!");

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);
      instNum  = institutInqCenter(vlistInqVarInstitut(vlistID1, varID));
      tableNum = tableInqNum(vlistInqVarTable(vlistID1, varID));

      code = vlistInqVarCode(vlistID1, varID);

      if ( useTable )
	{
	  if ( tableNum == 2 )
	    {
	      mode = WMO_MODE;
	      geop_code  =   6;
	      temp_code  =  11;
	      ps_code    =   1;
	    }
	  else if ( tableNum == 128 )
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
	cdoPrint("Mode = %d  Center = %d  Table = %d  Code = %d", mode, instNum, tableNum, code);

      if ( code <= 0 )
	{
	  vlistInqVarName(vlistID1, varID, varname);

	  strtolower(varname);

	  /*                        ECHAM                            ECMWF       */
	  if      ( strcmp(varname, "geosp") == 0 || strcmp(varname, "z")    == 0 ) code = 129;
	  else if ( strcmp(varname, "st")    == 0 || strcmp(varname, "t")    == 0 ) code = 130;
	  else if ( strcmp(varname, "aps")   == 0 || strcmp(varname, "sp"  ) == 0 ) code = 134;
	  else if ( strcmp(varname, "lsp")   == 0 || strcmp(varname, "lnsp") == 0 ) code = 152;
	  /* else if ( strcmp(varname, "geopoth") == 0 ) code = 156; */
	}

      if ( mode == ECHAM_MODE )
	{
	  if      ( code == geop_code  && nlevel == 1     ) geopID  = varID;
	  else if ( code == temp_code  && nlevel == nhlev ) tempID  = varID;
	  else if ( code == ps_code    && nlevel == 1     ) psID    = varID;
	  else if ( code == lsp_code   && nlevel == 1     ) lnpsID  = varID;
	  /* else if ( code == 156 ) gheightID = varID; */
	}
      else if ( mode == WMO_MODE )
	{
	  if      ( code == geop_code  && nlevel == 1     ) geopID  = varID;
	  else if ( code == temp_code  && nlevel == nhlev ) tempID  = varID;
	  else if ( code == ps_code    && nlevel == 1     ) psID    = varID;
	}
    }

  pvarID = lnpsID;
  /* Log. surface pressure is spectral, use the surface pressure instead */
  lnpsID = -1;
  if ( zaxisIDh != -1 && lnpsID == -1 )
    {
      if ( psID != -1 )
	{
	  pvarID = psID;
	  code = vlistInqVarCode(vlistID1, psID);
	  /* cdoPrint("LOG surface pressure not found - using surface pressure (code %d)!", code); */
	}
      else
	cdoAbort("Surface pressure not found!");
    }

  gridID   = vlistInqVarGrid(vlistID1, pvarID);
  if ( gridInqType(gridID) == GRID_SPECTRAL )
    cdoAbort("Surface pressure on spectral representation not supported!");
    
  gridsize = gridInqSize(gridID);
  pdata = (double *) malloc(gridsize*sizeof(double));


  vlistID2 = vlistCreate();
  varID = vlistDefVar(vlistID2, gridID, zaxisIDp, TIME_VARIABLE);
  vlistDefVarCode(vlistID2, varID, 1);
  vlistDefVarName(vlistID2, varID, "pressure");
  vlistDefVarLongname(vlistID2, varID, "Air pressure");
  vlistDefVarStdname(vlistID2, varID, "air_pressure");
  vlistDefVarUnits(vlistID2, varID, "Pa");

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);


  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( varID == pvarID )
	    {	  
	      streamReadRecord(streamID1, pdata, &nmiss);
	      if ( nmiss > 0 ) cdoAbort("Missing values unsupported!");
	    }
	}

      if ( zaxisIDh != -1 )
	{
	  if ( lnpsID != -1 )
	    for ( i = 0; i < ngp; i++ ) ps_prog[i] = exp(pdata[i]);
	  else if ( psID != -1 )
	    memcpy(ps_prog, pdata, ngp*sizeof(double));

	  /* check range of ps_prog */
	  {
	    double minval = ps_prog[0];
	    double maxval = ps_prog[0];
	    for ( i = 1; i < ngp; i++ )
	      {
		if      ( ps_prog[i] > maxval ) maxval = ps_prog[i];
		else if ( ps_prog[i] < minval ) minval = ps_prog[i];
	      }

	    if ( minval < 20000 || maxval > 150000 )
	      cdoWarning("Surface pressure out of range (min=%g max=%g)!", minval, maxval);
	  }

	  presh(full_press, half_press, vct, ps_prog, nhlevf, ngp);
	}

      if ( operatorID == FPRESSURE )
	{
	  nlevel = nhlevf;
	  hyb_press = full_press;
	}
      else
	{
	  nlevel = nhlevh;
	  hyb_press = half_press;
	}
	  
      varID = 0;
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  streamDefRecord(streamID2, varID, levelID);
	  offset = levelID*gridsize;
	  streamWriteRecord(streamID2, hyb_press+offset, 0);
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( pdata      ) free(pdata);
  if ( ps_prog    ) free(ps_prog);
  if ( full_press ) free(full_press);
  if ( half_press ) free(half_press);
  if ( vct        ) free(vct);
  if ( rvct       ) free(rvct);

  listDelete(flist);

  cdoFinish();

  return (0);
}
