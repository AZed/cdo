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

      Seasstat   seasmin         Seasonal minimum
      Seasstat   seasmax         Seasonal maximum
      Seasstat   seassum         Seasonal sum
      Seasstat   seasmean        Seasonal mean
      Seasstat   seasavg         Seasonal average
      Seasstat   seasvar         Seasonal variance
      Seasstat   seasstd         Seasonal standard deviation
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


void *Seasstat(void *argument)
{
  static char func[] = "Seasstat";
  int operatorID;
  int operfunc;
  int gridsize;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int vdate1 = 0, vtime1 = 0;
  int nrecs, nrecords;
  int gridID, varID, levelID, recID;
  int tsID;
  int otsID;
  long nsets;
  int i;
  int year, month, seas, seas0 = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  int newseas, oldmon = 0, newmon;
  int nseason = 0;
  double missval;
  field_t **vars1 = NULL, **vars2 = NULL, **samp1 = NULL;
  field_t field;
  int season_start;
  const char *seas_name[4];

  cdoInitialize(argument);

  cdoOperatorAdd("seasmin",  func_min,  0, NULL);
  cdoOperatorAdd("seasmax",  func_max,  0, NULL);
  cdoOperatorAdd("seassum",  func_sum,  0, NULL);
  cdoOperatorAdd("seasmean", func_mean, 0, NULL);
  cdoOperatorAdd("seasavg",  func_avg,  0, NULL);
  cdoOperatorAdd("seasvar",  func_var,  0, NULL);
  cdoOperatorAdd("seasstd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  season_start = get_season_start();
  get_season_name(seas_name);

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

  vars1 = (field_t **) malloc(nvars*sizeof(field_t *));
  samp1 = (field_t **) malloc(nvars*sizeof(field_t *));
  if ( operfunc == func_std || operfunc == func_var )
    vars2 = (field_t **) malloc(nvars*sizeof(field_t *));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval  = vlistInqVarMissval(vlistID1, varID);

      vars1[varID] = (field_t *)  malloc(nlevel*sizeof(field_t));
      samp1[varID] = (field_t *)  malloc(nlevel*sizeof(field_t));
      if ( operfunc == func_std || operfunc == func_var )
	vars2[varID] = (field_t *)  malloc(nlevel*sizeof(field_t));

      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  vars1[varID][levelID].grid    = gridID;
	  vars1[varID][levelID].nmiss   = 0;
	  vars1[varID][levelID].missval = missval;
	  vars1[varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
	  samp1[varID][levelID].grid    = gridID;
	  samp1[varID][levelID].nmiss   = 0;
	  samp1[varID][levelID].missval = missval;
	  samp1[varID][levelID].ptr     = NULL;
	  if ( operfunc == func_std || operfunc == func_var )
	    {
	      vars2[varID][levelID].grid    = gridID;
	      vars2[varID][levelID].nmiss   = 0;
	      vars2[varID][levelID].missval = missval;
	      vars2[varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
	    }
	}
    }

  tsID    = 0;
  otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      newseas = FALSE;
      while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
	{
	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);
	  year  =  vdate / 10000;
	  month = (vdate - year*10000) / 100;
	  if ( month < 1 || month > 12 )
	    cdoAbort("Month %d out of range!", month);

	  newmon = month;

	  if ( season_start == START_DEC )
	    {
	      if ( newmon == 12 ) newmon = 0;

	      if ( month <= 12 )
		seas = (month % 12) / 3;
	      else
		seas = month - 13;
	    }
	  else
	    {
	      if ( month <= 12 )
		seas = (month - 1) / 3;
	      else
		seas = month - 13;
	    }

	  if ( seas < 0 || seas > 3 )
	    cdoAbort("Season %d out of range!", seas+1);

	  if ( nsets == 0 )
	    {
	      nseason++;
	      vdate0 = vdate;
	      vtime0 = vtime;
	      seas0  = seas;
	      oldmon = newmon;
	    }

	  if ( newmon < oldmon ) newseas = TRUE;

	  if ( (seas != seas0) || newseas ) break;

	  oldmon = newmon;

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);

	      if ( tsID == 0 )
		{
		  recVarID[recID]   = varID;
		  recLevelID[recID] = levelID;
		}

	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	      if ( nsets == 0 )
		{
		  streamReadRecord(streamID1, vars1[varID][levelID].ptr, &nmiss);
		  vars1[varID][levelID].nmiss = nmiss;

		  if ( nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			samp1[varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));

		      for ( i = 0; i < gridsize; i++ )
			if ( DBL_IS_EQUAL(vars1[varID][levelID].ptr[i],
					  vars1[varID][levelID].missval) )
			  samp1[varID][levelID].ptr[i] = 0;
			else
			  samp1[varID][levelID].ptr[i] = 1;
		    }
		}
	      else
		{
		  streamReadRecord(streamID1, field.ptr, &field.nmiss);
		  field.grid    = vars1[varID][levelID].grid;
		  field.missval = vars1[varID][levelID].missval;

		  if ( field.nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			{
			  samp1[varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));
			  for ( i = 0; i < gridsize; i++ )
			    samp1[varID][levelID].ptr[i] = nsets;
			}

		      for ( i = 0; i < gridsize; i++ )
			if ( !DBL_IS_EQUAL(field.ptr[i], vars1[varID][levelID].missval) )
			  samp1[varID][levelID].ptr[i]++;
		    }

		  if ( operfunc == func_std || operfunc == func_var )
		    {
		      farsumq(&vars2[varID][levelID], field);
		      farsum(&vars1[varID][levelID], field);
		    }
		  else
		    {
		      farfun(&vars1[varID][levelID], field, operfunc);
		    }
		}
	    }

	  if ( nsets == 0 && (operfunc == func_std || operfunc == func_var) )
	    for ( varID = 0; varID < nvars; varID++ )
	      {
		if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
		nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		for ( levelID = 0; levelID < nlevel; levelID++ )
		  farmoq(&vars2[varID][levelID], vars1[varID][levelID]);
	      }

	  vdate1 = vdate;
	  vtime1 = vtime;
	  nsets++;
	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      if ( operfunc == func_mean || operfunc == func_avg )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[varID][levelID].ptr == NULL )
		  farcmul(&vars1[varID][levelID], 1.0/nsets);
		else
		  fardiv(&vars1[varID][levelID], samp1[varID][levelID]);
	      }
	  }
      else if ( operfunc == func_std || operfunc == func_var )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[varID][levelID].ptr == NULL )
		  {
		    if ( operfunc == func_std )
		      farcstd(&vars1[varID][levelID], vars2[varID][levelID], 1.0/nsets);
		    else
		      farcvar(&vars1[varID][levelID], vars2[varID][levelID], 1.0/nsets);
		  }
		else
		  {
		    farinv(&samp1[varID][levelID]);
		    if ( operfunc == func_std )
		      farstd(&vars1[varID][levelID], vars2[varID][levelID], samp1[varID][levelID]);
		    else
		      farvar(&vars1[varID][levelID], vars2[varID][levelID], samp1[varID][levelID]);
		  }
	      }
	  }

      if ( cdoVerbose )
	{
	  char vdatestr0[32], vtimestr0[32];
	  char vdatestr1[32], vtimestr1[32];
	  date2str(vdate0, vdatestr0, sizeof(vdatestr0));
	  time2str(vtime0, vtimestr0, sizeof(vtimestr0));
	  date2str(vdate1, vdatestr1, sizeof(vdatestr1));
	  time2str(vtime1, vtimestr1, sizeof(vtimestr1));
	  cdoPrint("season %3d %3s start %s %s end %s %s ntimesteps %d", 
		   nseason, seas_name[seas0],
		   vdatestr0, vtimestr0, vdatestr1, vtimestr1, nsets);
	}

      taxisDefVdate(taxisID2, vdate1);
      taxisDefVtime(taxisID2, vtime1);
      streamDefTimestep(streamID2, otsID++);

      if ( nsets < 3 )
	{
	  char vdatestr[32];
	  date2str(vdate0, vdatestr, sizeof(vdatestr));
	  cdoWarning("Season %3d (%s) has only %d input time step%s!", 
		     otsID, vdatestr, nsets, nsets == 1 ? "" : "s");
	}

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  varID   = recVarID[recID];
	  levelID = recLevelID[recID];

	  if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
	    {
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, vars1[varID][levelID].ptr,  vars1[varID][levelID].nmiss);
	    }
	}

      if ( nrecs == 0 ) break;
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  free(vars1[varID][levelID].ptr);
	  if ( samp1[varID][levelID].ptr ) free(samp1[varID][levelID].ptr);
	  if ( operfunc == func_std || operfunc == func_var ) free(vars2[varID][levelID].ptr);
	}

      free(vars1[varID]);
      free(samp1[varID]);
      if ( operfunc == func_std || operfunc == func_var ) free(vars2[varID]);
    }

  free(vars1);
  free(samp1);
  if ( operfunc == func_std || operfunc == func_var ) free(vars2);

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
