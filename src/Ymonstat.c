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

      Ymonstat   ymonmin         Multi-year monthly minimum
      Ymonstat   ymonmax         Multi-year monthly maximum
      Ymonstat   ymonsum         Multi-year monthly sum
      Ymonstat   ymonmean        Multi-year monthly mean
      Ymonstat   ymonavg         Multi-year monthly average
      Ymonstat   ymonvar         Multi-year monthly variance
      Ymonstat   ymonstd         Multi-year monthly standard deviation
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "field.h"


#define  NMONTH     17

void *Ymonstat(void *argument)
{
  static char func[] = "Ymonstat";
  int operatorID;
  int operfunc;
  int gridsize;
  int i;
  int varID;
  int recID;
  int gridID;
  int vdate, vtime;
  int year, month, day;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  long nsets[NMONTH];
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  int vdates[NMONTH], vtimes[NMONTH];
  int mon[NMONTH];
  int nmon = 0;
  double missval;
  field_t **vars1[NMONTH], **vars2[NMONTH], **samp1[NMONTH];
  field_t field;

  cdoInitialize(argument);

  cdoOperatorAdd("ymonmin",  func_min,  0, NULL);
  cdoOperatorAdd("ymonmax",  func_max,  0, NULL);
  cdoOperatorAdd("ymonsum",  func_sum,  0, NULL);
  cdoOperatorAdd("ymonmean", func_mean, 0, NULL);
  cdoOperatorAdd("ymonavg",  func_avg,  0, NULL);
  cdoOperatorAdd("ymonvar",  func_var,  0, NULL);
  cdoOperatorAdd("ymonstd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  for ( month = 0; month < NMONTH; month++ )
    {
      vars1[month] = NULL;
      vars2[month] = NULL;
      samp1[month] = NULL;
      nsets[month] = 0;
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

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
      if ( month < 0 || month >= NMONTH )
	cdoAbort("month %d out of range!", month);

      vdates[month] = vdate;
      vtimes[month] = vtime;

      if ( vars1[month] == NULL )
	{
	  mon[nmon++] = month;
	  vars1[month] = (field_t **) malloc(nvars*sizeof(field_t *));
	  samp1[month] = (field_t **) malloc(nvars*sizeof(field_t *));
	  if ( operfunc == func_std || operfunc == func_var )
	    vars2[month] = (field_t **) malloc(nvars*sizeof(field_t *));

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      missval  = vlistInqVarMissval(vlistID1, varID);

	      vars1[month][varID] = (field_t *)  malloc(nlevel*sizeof(field_t));
	      samp1[month][varID] = (field_t *)  malloc(nlevel*sizeof(field_t));
	      if ( operfunc == func_std || operfunc == func_var )
		vars2[month][varID] = (field_t *)  malloc(nlevel*sizeof(field_t));
	      
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  vars1[month][varID][levelID].grid    = gridID;
		  vars1[month][varID][levelID].nmiss   = 0;
		  vars1[month][varID][levelID].missval = missval;
		  vars1[month][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		  samp1[month][varID][levelID].grid    = gridID;
		  samp1[month][varID][levelID].nmiss   = 0;
		  samp1[month][varID][levelID].missval = missval;
		  samp1[month][varID][levelID].ptr     = NULL;
		  if ( operfunc == func_std || operfunc == func_var )
		    {
		      vars2[month][varID][levelID].grid    = gridID;
		      vars2[month][varID][levelID].nmiss   = 0;
		      vars2[month][varID][levelID].missval = missval;
		      vars2[month][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		    }
		}
	    }
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  recVarID[recID]   = varID;
	  recLevelID[recID] = levelID;

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  if ( nsets[month] == 0 )
	    {
	      streamReadRecord(streamID1, vars1[month][varID][levelID].ptr, &nmiss);
	      vars1[month][varID][levelID].nmiss = nmiss;

	      if ( nmiss > 0 || samp1[month][varID][levelID].ptr )
		{
		  if ( samp1[month][varID][levelID].ptr == NULL )
		    samp1[month][varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));

		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(vars1[month][varID][levelID].ptr[i],
				      vars1[month][varID][levelID].missval) )
		      samp1[month][varID][levelID].ptr[i] = 0;
		    else
		      samp1[month][varID][levelID].ptr[i] = 1;
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &field.nmiss);
	      field.grid    = vars1[month][varID][levelID].grid;
	      field.missval = vars1[month][varID][levelID].missval;

	      if ( field.nmiss > 0 || samp1[month][varID][levelID].ptr )
		{
		  if ( samp1[month][varID][levelID].ptr == NULL )
		    {
		      samp1[month][varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));
		      for ( i = 0; i < gridsize; i++ )
			samp1[month][varID][levelID].ptr[i] = nsets[month];
		    }
		  
		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], vars1[month][varID][levelID].missval) )
		      samp1[month][varID][levelID].ptr[i]++;
		}

	      if ( operfunc == func_std || operfunc == func_var )
		{
		  farsumq(&vars2[month][varID][levelID], field);
		  farsum(&vars1[month][varID][levelID], field);
		}
	      else
		{
		  farfun(&vars1[month][varID][levelID], field, operfunc);
		}
	    }
	}

      if ( nsets[month] == 0 && (operfunc == func_std || operfunc == func_var) )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      farmoq(&vars2[month][varID][levelID], vars1[month][varID][levelID]);
	  }

      nsets[month]++;
      tsID++;
    }

  for ( i = 0; i < nmon; i++ )
    {
      month = mon[i];
      if ( nsets[month] == 0 ) cdoAbort("Internal problem, nsets[%d] not set!", month);

      if ( operfunc == func_mean || operfunc == func_avg )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[month][varID][levelID].ptr == NULL )
		  farcmul(&vars1[month][varID][levelID], 1.0/nsets[month]);
		else
		  fardiv(&vars1[month][varID][levelID], samp1[month][varID][levelID]);
	      }
	  }
      else if ( operfunc == func_std || operfunc == func_var )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[month][varID][levelID].ptr == NULL )
		  {
		    if ( operfunc == func_std )
		      farcstd(&vars1[month][varID][levelID], vars2[month][varID][levelID], 1.0/nsets[month]);
		    else
		      farcvar(&vars1[month][varID][levelID], vars2[month][varID][levelID], 1.0/nsets[month]);
		  }
		else
		  {
		    farinv(&samp1[month][varID][levelID]);
		    if ( operfunc == func_std )
		      farstd(&vars1[month][varID][levelID], vars2[month][varID][levelID], samp1[month][varID][levelID]);
		    else
		      farvar(&vars1[month][varID][levelID], vars2[month][varID][levelID], samp1[month][varID][levelID]);
		  }
	      }
	  }

      taxisDefVdate(taxisID2, vdates[month]);
      taxisDefVtime(taxisID2, vtimes[month]);
      streamDefTimestep(streamID2, otsID++);

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  varID    = recVarID[recID];
	  levelID  = recLevelID[recID];
	  
	  if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
	    {
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, vars1[month][varID][levelID].ptr,
				vars1[month][varID][levelID].nmiss);
	    }
	}
    }

  for ( month = 0; month < NMONTH; month++ )
    {
      if ( vars1[month] != NULL )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  free(vars1[month][varID][levelID].ptr);
		  if ( samp1[month][varID][levelID].ptr ) free(samp1[month][varID][levelID].ptr);
		  if ( operfunc == func_std || operfunc == func_var ) free(vars2[month][varID][levelID].ptr);
		}
	      
	      free(vars1[month][varID]);
	      free(samp1[month][varID]);
	      if ( operfunc == func_std || operfunc == func_var ) free(vars2[month][varID]);
	    }

	  free(vars1[month]);
	  free(samp1[month]);
	  if ( operfunc == func_std || operfunc == func_var ) free(vars2[month]);
	}
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  cdoFinish();

  return (0);
}
