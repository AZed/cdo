/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2011 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Yhourstat   yhourmin         Multi-year hourly minimum
      Yhourstat   yhourmax         Multi-year hourly maximum
      Yhourstat   yhoursum         Multi-year hourly sum
      Yhourstat   yhourmean        Multi-year hourly mean
      Yhourstat   yhouravg         Multi-year hourly average
      Yhourstat   yhourvar         Multi-year hourly variance
      Yhourstat   yhourstd         Multi-year hourly standard deviation
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  NHOUR       8952  /* 31*12*24 */


void *Yhourstat(void *argument)
{
  int operatorID;
  int operfunc;
  int gridsize;
  int i;
  int varID;
  int recID;
  int gridID;
  int vdate, vtime;
  int year, month, day, houroy;
  int hour, minute, second;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  long nsets[NHOUR];
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  int vdates[NHOUR], vtimes[NHOUR];
  double missval;
  field_t **vars1[NHOUR], **vars2[NHOUR], **samp1[NHOUR];
  field_t field;

  cdoInitialize(argument);

  cdoOperatorAdd("yhourmin",  func_min,  0, NULL);
  cdoOperatorAdd("yhourmax",  func_max,  0, NULL);
  cdoOperatorAdd("yhoursum",  func_sum,  0, NULL);
  cdoOperatorAdd("yhourmean", func_mean, 0, NULL);
  cdoOperatorAdd("yhouravg",  func_avg,  0, NULL);
  cdoOperatorAdd("yhourvar",  func_var,  0, NULL);
  cdoOperatorAdd("yhourstd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);

  for ( houroy = 0; houroy < NHOUR; houroy++ )
    {
      vars1[houroy] = NULL;
      vars2[houroy] = NULL;
      samp1[houroy] = NULL;
      nsets[houroy] = 0;
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);
  field.ptr = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  otsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      cdiDecodeDate(vdate, &year, &month, &day);
      cdiDecodeTime(vtime, &hour, &minute, &second);

      if ( month >= 1 && month <= 12 && hour >= 0 && hour < 24 )
	houroy = ((month-1)*31 + day - 1)*24 + hour;
      else
	houroy = 0;

      if ( houroy < 0 || houroy >= NHOUR )
	cdoAbort("hour of year %d out of range (date=%d time=%d)!", houroy, vdate, vtime);

      vdates[houroy] = vdate;
      vtimes[houroy] = vtime;

      if ( vars1[houroy] == NULL )
	{
	  vars1[houroy] = (field_t **) malloc(nvars*sizeof(field_t *));
	  samp1[houroy] = (field_t **) malloc(nvars*sizeof(field_t *));
	  if ( operfunc == func_std || operfunc == func_var )
	    vars2[houroy] = (field_t **) malloc(nvars*sizeof(field_t *));

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      missval  = vlistInqVarMissval(vlistID1, varID);

	      vars1[houroy][varID] = (field_t *)  malloc(nlevel*sizeof(field_t));
	      samp1[houroy][varID] = (field_t *)  malloc(nlevel*sizeof(field_t));
	      if ( operfunc == func_std || operfunc == func_var )
		vars2[houroy][varID] = (field_t *)  malloc(nlevel*sizeof(field_t));
	      
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  vars1[houroy][varID][levelID].grid    = gridID;
		  vars1[houroy][varID][levelID].nmiss   = 0;
		  vars1[houroy][varID][levelID].missval = missval;
		  vars1[houroy][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		  samp1[houroy][varID][levelID].grid    = gridID;
		  samp1[houroy][varID][levelID].nmiss   = 0;
		  samp1[houroy][varID][levelID].missval = missval;
		  samp1[houroy][varID][levelID].ptr     = NULL;
		  if ( operfunc == func_std || operfunc == func_var )
		    {
		      vars2[houroy][varID][levelID].grid    = gridID;
		      vars2[houroy][varID][levelID].nmiss   = 0;
		      vars2[houroy][varID][levelID].missval = missval;
		      vars2[houroy][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		    }
		}
	    }
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;
	    }

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  if ( nsets[houroy] == 0 )
	    {
	      streamReadRecord(streamID1, vars1[houroy][varID][levelID].ptr, &nmiss);
	      vars1[houroy][varID][levelID].nmiss = nmiss;

	      if ( nmiss > 0 || samp1[houroy][varID][levelID].ptr )
		{
		  if ( samp1[houroy][varID][levelID].ptr == NULL )
		    samp1[houroy][varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));

		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(vars1[houroy][varID][levelID].ptr[i],
				      vars1[houroy][varID][levelID].missval) )
		      samp1[houroy][varID][levelID].ptr[i] = 0;
		    else
		      samp1[houroy][varID][levelID].ptr[i] = 1;
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &field.nmiss);
	      field.grid    = vars1[houroy][varID][levelID].grid;
	      field.missval = vars1[houroy][varID][levelID].missval;

	      if ( field.nmiss > 0 || samp1[houroy][varID][levelID].ptr )
		{
		  if ( samp1[houroy][varID][levelID].ptr == NULL )
		    {
		      samp1[houroy][varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));
		      for ( i = 0; i < gridsize; i++ )
			samp1[houroy][varID][levelID].ptr[i] = nsets[houroy];
		    }
		  
		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], vars1[houroy][varID][levelID].missval) )
		      samp1[houroy][varID][levelID].ptr[i]++;
		}

	      if ( operfunc == func_std || operfunc == func_var )
		{
		  farsumq(&vars2[houroy][varID][levelID], field);
		  farsum(&vars1[houroy][varID][levelID], field);
		}
	      else
		{
		  farfun(&vars1[houroy][varID][levelID], field, operfunc);
		}
	    }
	}

      if ( nsets[houroy] == 0 && (operfunc == func_std || operfunc == func_var) )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      farmoq(&vars2[houroy][varID][levelID], vars1[houroy][varID][levelID]);
	  }

      nsets[houroy]++;
      tsID++;
    }

  for ( houroy = 0; houroy < NHOUR; houroy++ )
    if ( nsets[houroy] )
      {
	if ( operfunc == func_mean || operfunc == func_avg )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( samp1[houroy][varID][levelID].ptr == NULL )
		    farcmul(&vars1[houroy][varID][levelID], 1.0/nsets[houroy]);
		  else
		    fardiv(&vars1[houroy][varID][levelID], samp1[houroy][varID][levelID]);
		}
	    }
	else if ( operfunc == func_std || operfunc == func_var )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( samp1[houroy][varID][levelID].ptr == NULL )
		    {
		      if ( operfunc == func_std )
			farcstd(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], 1.0/nsets[houroy]);
		      else
			farcvar(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], 1.0/nsets[houroy]);
		    }
		  else
		    {
		      farinv(&samp1[houroy][varID][levelID]);
		      if ( operfunc == func_std )
			farstd(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], samp1[houroy][varID][levelID]);
		      else
			farvar(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], samp1[houroy][varID][levelID]);
		    }
		}
	    }

	taxisDefVdate(taxisID2, vdates[houroy]);
	taxisDefVtime(taxisID2, vtimes[houroy]);
	streamDefTimestep(streamID2, otsID);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID && vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;

	    streamDefRecord(streamID2, varID, levelID);
	    streamWriteRecord(streamID2, vars1[houroy][varID][levelID].ptr,
			      vars1[houroy][varID][levelID].nmiss);
	  }

	otsID++;
      }

  for ( houroy = 0; houroy < NHOUR; houroy++ )
    {
      if ( vars1[houroy] != NULL )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  free(vars1[houroy][varID][levelID].ptr);
		  if ( samp1[houroy][varID][levelID].ptr ) free(samp1[houroy][varID][levelID].ptr);
		  if ( operfunc == func_std || operfunc == func_var ) free(vars2[houroy][varID][levelID].ptr);
		}
	      
	      free(vars1[houroy][varID]);
	      free(samp1[houroy][varID]);
	      if ( operfunc == func_std || operfunc == func_var ) free(vars2[houroy][varID]);
	    }

	  free(samp1[houroy]);
	  free(vars1[houroy]);
	  if ( operfunc == func_std || operfunc == func_var ) free(vars2[houroy]);
	}
    }

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
