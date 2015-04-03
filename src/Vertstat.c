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

      Vertstat   vertmin         Vertical minimum
      Vertstat   vertmax         Vertical maximum
      Vertstat   vertsum         Vertical sum
      Vertstat   vertmean        Vertical mean
      Vertstat   vertavg         Vertical average
      Vertstat   vertvar         Vertical variance
      Vertstat   vertstd         Vertical standard deviation
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "functs.h"


void *Vertstat(void *argument)
{
  static char func[] = "Vertstat";
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridsize;
  int recID, nrecs;
  int gridID;
  int i;
  int tsID, varID, levelID;
  int nmiss, nvars;
  int zaxisID, nzaxis;
  double missval;
  field_t *vars1 = NULL, *vars2 = NULL, *samp1 = NULL;
  field_t field;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  cdoOperatorAdd("vertmin",  func_min,  0, NULL);
  cdoOperatorAdd("vertmax",  func_max,  0, NULL);
  cdoOperatorAdd("vertsum",  func_sum,  0, NULL);
  cdoOperatorAdd("vertmean", func_mean, 0, NULL);
  cdoOperatorAdd("vertavg",  func_avg,  0, NULL);
  cdoOperatorAdd("vertvar",  func_var,  0, NULL);
  cdoOperatorAdd("vertstd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  vlistClearFlag(vlistID1);
  nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    vlistDefFlag(vlistID1, varID, 0, TRUE);

  vlistID2 = vlistCreate();
  vlistCopyFlag(vlistID2, vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
  nzaxis  = vlistNzaxis(vlistID1);
  for ( i = 0; i < nzaxis; i++ )
    if ( zaxisInqSize(vlistZaxis(vlistID1, i)) > 1 )
      vlistChangeZaxisIndex(vlistID2, i, zaxisID);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);

  field.ptr = (double *) malloc(gridsize*sizeof(double));

  vars1 = (field_t *) malloc(nvars*sizeof(field_t));
  samp1 = (field_t *) malloc(nvars*sizeof(field_t));
  if ( operfunc == func_std || operfunc == func_var )
    vars2 = (field_t *) malloc(nvars*sizeof(field_t));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      missval  = vlistInqVarMissval(vlistID1, varID);

      vars1[varID].grid    = gridID;
      vars1[varID].nsamp   = 0;
      vars1[varID].nmiss   = 0;
      vars1[varID].missval = missval;
      vars1[varID].ptr     = (double *) malloc(gridsize*sizeof(double));
      samp1[varID].grid    = gridID;
      samp1[varID].nmiss   = 0;
      samp1[varID].missval = missval;
      samp1[varID].ptr     = NULL;
      if ( operfunc == func_std || operfunc == func_var )
	{
	  vars2[varID].grid    = gridID;
	  vars2[varID].nmiss   = 0;
	  vars2[varID].missval = missval;
	  vars2[varID].ptr     = (double *) malloc(gridsize*sizeof(double));
	}
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
          vars1[varID].nsamp++;
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  if ( levelID == 0 )
	    {
	      streamReadRecord(streamID1, vars1[varID].ptr, &nmiss);
	      vars1[varID].nmiss = nmiss;

	      if ( operfunc == func_std || operfunc == func_var )
		farmoq(&vars2[varID], vars1[varID]);

	      if ( nmiss > 0 || samp1[varID].ptr )
		{
		  if ( samp1[varID].ptr == NULL )
		    samp1[varID].ptr = (double *) malloc(gridsize*sizeof(double));

		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(vars1[varID].ptr[i], vars1[varID].missval) )
		      samp1[varID].ptr[i] = 0;
		    else
		      samp1[varID].ptr[i] = 1;
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &field.nmiss);
	      field.grid    = vars1[varID].grid;
	      field.missval = vars1[varID].missval;

	      if ( field.nmiss > 0 || samp1[varID].ptr )
		{
		  if ( samp1[varID].ptr == NULL )
		    {
		      samp1[varID].ptr = (double *) malloc(gridsize*sizeof(double));
		      for ( i = 0; i < gridsize; i++ )
			samp1[varID].ptr[i] = vars1[varID].nsamp;
		    }

		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], vars1[varID].missval) )
		      samp1[varID].ptr[i]++;
		}

	      if ( operfunc == func_std || operfunc == func_var )
		{
		  farsumq(&vars2[varID], field);
		  farsum(&vars1[varID], field);
		}
	      else
		{
		  farfun(&vars1[varID], field, operfunc);
		}
	    }
	}

      if ( operfunc == func_mean || operfunc == func_avg )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( samp1[varID].ptr == NULL )
	      farcmul(&vars1[varID], 1.0/vars1[varID].nsamp);
	    else
	      fardiv(&vars1[varID], samp1[varID]);
	  }
      else if ( operfunc == func_std || operfunc == func_var )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( samp1[varID].ptr == NULL )
	      {
		if ( operfunc == func_std )
		  farcstd(&vars1[varID], vars2[varID], 1.0/vars1[varID].nsamp);
		else
		  farcvar(&vars1[varID], vars2[varID], 1.0/vars1[varID].nsamp);
	      }
	    else
	      {
		farinv(&samp1[varID]);
		if ( operfunc == func_std )
		  farstd(&vars1[varID], vars2[varID], samp1[varID]);
		else
		  farvar(&vars1[varID], vars2[varID], samp1[varID]);
	      }
	  }

      for ( varID = 0; varID < nvars; varID++ )
	{
	  streamDefRecord(streamID2, varID, 0);
	  streamWriteRecord(streamID2, vars1[varID].ptr, vars1[varID].nmiss);
	  vars1[varID].nsamp = 0;
	}

      tsID++;
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      free(vars1[varID].ptr);
      if ( samp1[varID].ptr ) free(samp1[varID].ptr);
      if ( operfunc == func_std || operfunc == func_var ) free(vars2[varID].ptr);
    }

  free(vars1);
  free(samp1);
  if ( operfunc == func_std || operfunc == func_var ) free(vars2);

  if ( field.ptr ) free(field.ptr);

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  cdoFinish();

  return (0);
}
