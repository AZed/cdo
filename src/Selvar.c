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

      Selvar     selparam        Select parameters by identifier (format: code.tabnum  or  pnum.cat.dis)
      Selvar     delparam        Delete parameters by identifier (format: code.tabnum  or  pnum.cat.dis)
      Selvar     selcode         Select parameters by code number
      Selvar     delcode         Delete parameters by code number
      Selvar     selname         Select parameters by name
      Selvar     delname         Delete parameters by name
      Selvar     selstdname      Select parameters by CF standard name
      Selvar     sellevel        Select levels
      Selvar     sellevidx       Select levels by index
      Selvar     selgrid         Select grids
      Selvar     selzaxis        Select zaxis
      Selvar     seltabnum       Select parameter table number
      Selvar     selltype        Select GRIB level type 
*/

#include <ctype.h>  /* isdigit */

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "error.h"
#include "util.h"
#include "list.h"


void *Selvar(void *argument)
{
  int SELPARAM, SELCODE, SELNAME, SELLEVEL, SELLEVIDX, SELGRID, SELZAXIS, SELLTYPE; 
  int SELTABNUM, DELPARAM, DELCODE, DELNAME, SELSTDNAME;
  int operatorID;
  int streamID1, streamID2;
  int tsID, nrecs;
  int nvars, nlevs;
  int code, tabnum, param, gridID, zaxisID, levID;
  int grididx, zaxisidx;
  double level;
  int varID2, levelID2;
  int recID, varID, levelID;
  int *intarr = NULL, nsel = 0;
  int *selfound = NULL;
  double *fltarr = NULL;
  char paramstr[32];
  char varname[CDI_MAX_NAME];
  char stdname[CDI_MAX_NAME];
  char gridname[CDI_MAX_NAME];
  char zaxisname[CDI_MAX_NAME];
  char **argnames = NULL;
  int vlistID1 = -1, vlistID2 = -1;
  int isel;
  int i;
  int npar;
  int intlist, byname;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  double *array = NULL;
  int taxisID1, taxisID2;
  int ltype;
  LIST *ilist = listNew(INT_LIST);
  LIST *flist = listNew(FLT_LIST);

  cdoInitialize(argument);

  SELPARAM     = cdoOperatorAdd("selparam",     0, 0, "parameters");
  SELCODE      = cdoOperatorAdd("selcode",      0, 0, "code numbers");
  SELNAME      = cdoOperatorAdd("selname",      0, 0, "variable names");
  SELSTDNAME   = cdoOperatorAdd("selstdname",   0, 0, "standard names");
  SELLEVEL     = cdoOperatorAdd("sellevel",     0, 0, "levels");
  SELLEVIDX    = cdoOperatorAdd("sellevidx",    0, 0, "index of levels");
  SELGRID      = cdoOperatorAdd("selgrid",      0, 0, "list of grid names or numbers");
  SELZAXIS     = cdoOperatorAdd("selzaxis",     0, 0, "list of zaxis names or numbers");
  SELTABNUM    = cdoOperatorAdd("seltabnum",    0, 0, "table numbers");
  DELPARAM     = cdoOperatorAdd("delparam",     0, 0, "parameter");
  DELCODE      = cdoOperatorAdd("delcode",      0, 0, "code numbers");
  DELNAME      = cdoOperatorAdd("delname",      0, 0, "variable names");
  SELLTYPE     = cdoOperatorAdd("selltype",     0, 0, "GRIB level types"); 

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  intlist = FALSE;
  byname  = TRUE;

  if ( operatorID == SELPARAM || operatorID == DELPARAM || operatorID == SELNAME || operatorID == DELNAME || 
       operatorID == SELSTDNAME || operatorID == SELGRID || operatorID == SELZAXIS )
    {
      nsel     = operatorArgc();
      argnames = operatorArgv();

      if ( cdoVerbose )
	for ( i = 0; i < nsel; i++ )
	  fprintf(stderr, "name %d = %s\n", i+1, argnames[i]);

      if ( operatorID == SELGRID || operatorID == SELZAXIS )
	{
	  if ( nsel > 0 && isdigit(*argnames[0]) )
	    {
	      intlist = TRUE;
	      byname  = FALSE;
	    }
	}
    }
  else if ( operatorID == SELLEVEL )
    {
      nsel = args2fltlist(operatorArgc(), operatorArgv(), flist);
      fltarr = (double *) listArrayPtr(flist);

      if ( cdoVerbose )
	for ( i = 0; i < nsel; i++ )
	  printf("flt %d = %g\n", i+1, fltarr[i]);
    }
  else
    {
      intlist = TRUE;
    }

  if ( intlist )
    {
      nsel = args2intlist(operatorArgc(), operatorArgv(), ilist);
      intarr = (int *) listArrayPtr(ilist);

      if ( cdoVerbose )
	for ( i = 0; i < nsel; i++ )
	  printf("int %d = %d\n", i+1, intarr[i]);
    }

  if ( nsel )
    {
      selfound = malloc(nsel*sizeof(int));
      for ( i = 0; i < nsel; i++ ) selfound[i] = FALSE;
    }

  /*
  if ( nsel == 0 )
    cdoAbort("missing code argument!");
  */
  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  vlistClearFlag(vlistID1);
  nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    {
      vlistInqVarName(vlistID1, varID, varname);
      vlistInqVarStdname(vlistID1, varID, stdname);
      param    = vlistInqVarParam(vlistID1, varID);
      code     = vlistInqVarCode(vlistID1, varID);
      tabnum   = tableInqNum(vlistInqVarTable(vlistID1, varID));
      gridID   = vlistInqVarGrid(vlistID1, varID);
      grididx  = vlistGridIndex(vlistID1, gridID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      zaxisidx = vlistZaxisIndex(vlistID1, zaxisID);
      nlevs    = zaxisInqSize(zaxisID);
      gridName(gridInqType(gridID), gridname);
      zaxisName(zaxisInqType(zaxisID), zaxisname);

      cdiParamToString(param, paramstr, sizeof(paramstr));

      for ( levID = 0; levID < nlevs; levID++ )
	{
	  level = zaxisInqLevel(zaxisID, levID);

	  if ( operatorID == DELCODE || operatorID == DELNAME || operatorID == DELPARAM )
	    vlistDefFlag(vlistID1, varID, levID, TRUE);

	  for ( isel = 0; isel < nsel; isel++ )
	    {
	      if ( operatorID == SELCODE )
		{
		  if ( intarr[isel] == code )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELPARAM )
		{
		  if ( strcmp(argnames[isel], paramstr) == 0 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELNAME )
		{
		  if ( strcmp(argnames[isel], varname) == 0 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELSTDNAME )
		{
		  if ( strcmp(argnames[isel], stdname) == 0 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELLEVEL )
		{
		  if ( fabs(fltarr[isel] - level) < 0.0001 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELLEVIDX )
		{
		  if ( intarr[isel] == (levID+1) )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELGRID && byname == FALSE )
		{
		  if ( intarr[isel] == (grididx+1) )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELGRID && byname == TRUE )
		{
		  if ( memcmp(argnames[isel], gridname, strlen(argnames[isel])) == 0 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELZAXIS && byname == FALSE )
		{
		  if ( intarr[isel] == (zaxisidx+1) )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELZAXIS && byname == TRUE )
		{
		  if ( memcmp(argnames[isel], zaxisname, strlen(argnames[isel])) == 0 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELTABNUM )
		{
		  if ( intarr[isel] == tabnum )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == DELCODE )
		{
		  if ( intarr[isel] == code )
		    {
		      vlistDefFlag(vlistID1, varID, levID, FALSE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == DELNAME )
		{
		  if ( strcmp(argnames[isel], varname) == 0 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, FALSE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == DELPARAM )
		{
		  if ( strcmp(argnames[isel], paramstr) == 0 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, FALSE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELLTYPE )
		{
		  ltype = zaxis2ltype(zaxisID);

		  if ( intarr[isel] == ltype )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	    }
	}
    }

  npar = 0;
  for ( varID = 0; varID < nvars; varID++ )
    {
      zaxisID = vlistInqVarZaxis(vlistID1, varID);
      nlevs   = zaxisInqSize(zaxisID);

      for ( levID = 0; levID < nlevs; levID++ )
	if ( vlistInqFlag(vlistID1, varID, levID) == TRUE ) break;
	      
      if ( levID < nlevs ) npar++;
    }

  for ( isel = 0; isel < nsel; isel++ )
    {
      if ( selfound[isel] == FALSE )
	{
	  if ( operatorID == SELCODE || operatorID == DELCODE )
	    {
	      cdoWarning("Code number %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELPARAM || operatorID == DELPARAM )
	    {
	      cdoWarning("Parameter %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELNAME || operatorID == DELNAME )
	    {
	      cdoWarning("Variable name %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELSTDNAME )
	    {
	      cdoWarning("Variable with standard name %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELLEVEL )
	    {
	      cdoWarning("Level %g not found!", fltarr[isel]);
	    }
	  else if ( operatorID == SELLEVIDX )
	    {
	      cdoWarning("Level index %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELGRID && byname == FALSE )
	    {
	      cdoWarning("Grid %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELGRID && byname == TRUE )
	    {
	      cdoWarning("Grid name %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELZAXIS && byname == FALSE )
	    {
	      cdoWarning("Zaxis %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELZAXIS && byname == TRUE )
	    {
	      cdoWarning("Zaxis name %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELTABNUM )
	    {
	      cdoWarning("Table number %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELLTYPE )
	    {
	      cdoWarning("GRIB level type %d not found!", intarr[isel]);
	    }
	}
    }
  
  if ( npar == 0 ) cdoAbort("No variables selected!");

  vlistID2 = vlistCreate();
  vlistCopyFlag(vlistID2, vlistID1);

  nvars = vlistNvars(vlistID2);
  for ( varID = 0; varID < nvars; ++varID )
    if ( vlistInqVarTsteptype(vlistID2, varID) != TSTEP_CONSTANT ) break;
  if ( varID == nvars ) vlistDefNtsteps(vlistID2, 0);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nrecs = vlistNrecs(vlistID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = malloc(gridsize*sizeof(double));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
     
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  if ( vlistInqFlag(vlistID1, varID, levelID) == TRUE )
	    {
	      varID2   = vlistFindVar(vlistID2, varID);
	      levelID2 = vlistFindLevel(vlistID2, varID, levelID);

	      streamDefRecord(streamID2, varID2, levelID2);
	      if ( lcopy )
		{
		  streamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  streamReadRecord(streamID1, array, &nmiss);
		  streamWriteRecord(streamID2, array, nmiss);
		}
     	    }
	}
       
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);
 
  vlistDestroy(vlistID2);

  if ( ! lcopy )
    if ( array ) free(array);

  if ( selfound ) free(selfound);

  listDelete(ilist);
  listDelete(flist);

  cdoFinish();

  return (NULL);
}
