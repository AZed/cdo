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

      Runstat    runmin          Running minimum
      Runstat    runmax          Running maximum
      Runstat    runsum          Running sum
      Runstat    runmean         Running mean
      Runstat    runavg          Running average
      Runstat    runvar          Running variance
      Runstat    runstd          Running standard deviation
*/

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "field.h"


void datetime_avg(int calendar, int ndates, datetime_t *datetime)
{
  int vdate, vtime;
  juldate_t juldate1, juldate2, juldatem;
  double seconds;
  /*
  for ( i = 0; i < ndates; i++ )
    fprintf(stdout, "%4d %d %d\n", i+1, datetime[i].date, datetime[i].time);
  */
  if ( ndates%2 == 0 )
    {
      /*
      vdate = datetime[ndates-1].date;
      vtime = datetime[ndates-1].time;
      */
      vdate = datetime[ndates/2-1].date;
      vtime = datetime[ndates/2-1].time;
      juldate1 = juldate_encode(calendar, vdate, vtime);

      vdate = datetime[ndates/2].date;
      vtime = datetime[ndates/2].time;
      juldate2 = juldate_encode(calendar, vdate, vtime);

      seconds = juldate_to_seconds(juldate_sub(juldate2, juldate1)) / 2;
      juldatem = juldate_add_seconds(NINT(seconds), juldate1);
      juldate_decode(calendar, juldatem, &vdate, &vtime);
    }
  else
    {
      vdate = datetime[ndates/2].date;
      vtime = datetime[ndates/2].time;
    }

  datetime[ndates].date = vdate;
  datetime[ndates].time = vtime;
  /*
  fprintf(stdout, "res: %d %d\n\n", datetime[ndates].date, datetime[ndates].time);
  */
}


void *Runstat(void *argument)
{
  static char func[] = "Runstat";
  int operatorID;
  int operfunc;
  int gridsize;
  int i;
  int varID;
  int recID;
  int gridID;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  int inp, its, ndates = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  double missval;
  field_t ***vars1 = NULL, ***vars2 = NULL, ***samp1 = NULL;
  datetime_t *datetime;
  int taxisID1, taxisID2;
  int calendar;
  int runstat_date = DATE_MIDDLE;
  char *envstr;

  cdoInitialize(argument);

  envstr = getenv("RUNSTAT_DATE");
  if ( envstr )
    {
      int env_date = -1;

      if      ( memcmp(envstr, "first", 5) == 0 ||
		memcmp(envstr, "FIRST", 5) == 0 ||
		memcmp(envstr, "First", 5) == 0 )  env_date = DATE_FIRST;
      else if ( memcmp(envstr, "last", 4) == 0 ||
		memcmp(envstr, "LAST", 4) == 0 ||
		memcmp(envstr, "Last", 4) == 0 )   env_date = DATE_LAST;
      else if ( memcmp(envstr, "middle", 6) == 0 ||
		memcmp(envstr, "MIDDLE", 6) == 0 ||
		memcmp(envstr, "Middle", 6) == 0 ) env_date = DATE_MIDDLE;

      if ( env_date >= 0 )
	{
	  runstat_date = env_date;

	  if ( cdoVerbose )
	    cdoPrint("Set RUNSTAT_DATE to %s", envstr);
	}
    }

  cdoOperatorAdd("runmin",  func_min,  0, NULL);
  cdoOperatorAdd("runmax",  func_max,  0, NULL);
  cdoOperatorAdd("runsum",  func_sum,  0, NULL);
  cdoOperatorAdd("runmean", func_mean, 0, NULL);
  cdoOperatorAdd("runavg",  func_avg,  0, NULL);
  cdoOperatorAdd("runvar",  func_var,  0, NULL);
  cdoOperatorAdd("runstd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  operatorInputArg("number of timesteps");
  ndates = atoi(operatorArgv()[0]);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  calendar = taxisInqCalendar(taxisID1);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  datetime = (datetime_t *) malloc((ndates+1)*sizeof(datetime_t));
  vars1 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));
  samp1 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));
  if ( operfunc == func_std || operfunc == func_var )
    vars2 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));

  for ( its = 0; its < ndates; its++ )
    {
      vars1[its] = (field_t **) malloc(nvars*sizeof(field_t *));
      samp1[its] = (field_t **) malloc(nvars*sizeof(field_t *));
      if ( operfunc == func_std || operfunc == func_var )
	vars2[its] = (field_t **) malloc(nvars*sizeof(field_t *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  missval  = vlistInqVarMissval(vlistID1, varID);

	  vars1[its][varID] = (field_t *)  malloc(nlevel*sizeof(field_t));
	  samp1[its][varID] = (field_t *)  malloc(nlevel*sizeof(field_t));
	  if ( operfunc == func_std || operfunc == func_var )
	    vars2[its][varID] = (field_t *)  malloc(nlevel*sizeof(field_t));

	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      vars1[its][varID][levelID].grid    = gridID;
	      vars1[its][varID][levelID].nmiss   = 0;
	      vars1[its][varID][levelID].missval = missval;
	      vars1[its][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
	      samp1[its][varID][levelID].grid    = gridID;
	      samp1[its][varID][levelID].nmiss   = 0;
	      samp1[its][varID][levelID].missval = missval;
	      samp1[its][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
	      if ( operfunc == func_std || operfunc == func_var )
		{
		  vars2[its][varID][levelID].grid    = gridID;
		  vars2[its][varID][levelID].nmiss   = 0;
		  vars2[its][varID][levelID].missval = missval;
		  vars2[its][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		}
	    }
	}
    }

  for ( tsID = 0; tsID < ndates; tsID++ )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 )
	cdoAbort("File has less then %d timesteps!", ndates);

      datetime[tsID].date = taxisInqVdate(taxisID1);
      datetime[tsID].time = taxisInqVtime(taxisID1);
	
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;
	    }
	  
	  streamReadRecord(streamID1, vars1[tsID][varID][levelID].ptr, &nmiss);
	  vars1[tsID][varID][levelID].nmiss = nmiss;

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  missval = vars1[ndates-1][varID][levelID].missval;

	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(vars1[tsID][varID][levelID].ptr[i], missval) )
	      samp1[tsID][varID][levelID].ptr[i] = 0;
	    else
	      {
		samp1[tsID][varID][levelID].ptr[i] = 1;
		for ( inp = 0; inp < tsID; inp++ )
		  samp1[inp][varID][levelID].ptr[i]++;
	      }

	  if ( operfunc == func_std || operfunc == func_var )
	    {
	      farmoq(&vars2[tsID][varID][levelID], vars1[tsID][varID][levelID]);
	      for ( inp = 0; inp < tsID; inp++ )
		{
		  farsumq(&vars2[inp][varID][levelID], vars1[tsID][varID][levelID]);
		  farsum(&vars1[inp][varID][levelID], vars1[tsID][varID][levelID]);
		}
	    }
	  else
	    {
	      for ( inp = 0; inp < tsID; inp++ )
		{
		  farfun(&vars1[inp][varID][levelID], vars1[tsID][varID][levelID], operfunc);
		}
	    }
	}
    }

  otsID = 0;
  while ( TRUE )
    {
      if ( operfunc == func_mean || operfunc == func_avg )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[0][varID][levelID].ptr == NULL )
		  farcmul(&vars1[0][varID][levelID], 1.0/ndates);
		else
		  fardiv(&vars1[0][varID][levelID], samp1[0][varID][levelID]);
	      }
	  }
      else if ( operfunc == func_std || operfunc == func_var )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[0][varID][levelID].ptr == NULL )
		  {
		    if ( operfunc == func_std )
		      farcstd(&vars1[0][varID][levelID], vars2[0][varID][levelID], 1.0/ndates);
		    else
		      farcvar(&vars1[0][varID][levelID], vars2[0][varID][levelID], 1.0/ndates);
		  }
		else
		  {
		    farinv(&samp1[0][varID][levelID]);
		    if ( operfunc == func_std )
		      farstd(&vars1[0][varID][levelID], vars2[0][varID][levelID], samp1[0][varID][levelID]);
		    else
		      farvar(&vars1[0][varID][levelID], vars2[0][varID][levelID], samp1[0][varID][levelID]);
		  }
	      }
	  }

      if ( runstat_date == DATE_MIDDLE )
	{
	  datetime_avg(calendar, ndates, datetime);
	}
      else if ( runstat_date == DATE_FIRST )
	{
	  datetime[ndates].date = datetime[0].date;
	  datetime[ndates].time = datetime[0].time;
	}
      else if ( runstat_date == DATE_LAST )
	{
	  datetime[ndates].date = datetime[ndates-1].date;
	  datetime[ndates].time = datetime[ndates-1].time;
	}

      taxisDefVdate(taxisID2, datetime[ndates].date);
      taxisDefVtime(taxisID2, datetime[ndates].time);
      streamDefTimestep(streamID2, otsID++);

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  varID    = recVarID[recID];
	  levelID  = recLevelID[recID];

	  if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
	    {
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, vars1[0][varID][levelID].ptr, vars1[0][varID][levelID].nmiss);
	    }
	}

      datetime[ndates] = datetime[0];
      vars1[ndates] = vars1[0];
      samp1[ndates] = samp1[0];
      if ( operfunc == func_std || operfunc == func_var )
        vars2[ndates] = vars2[0];

      for ( inp = 0; inp < ndates; inp++ )
	{
	  datetime[inp] = datetime[inp+1];
	  vars1[inp] = vars1[inp+1];
	  samp1[inp] = samp1[inp+1];
	  if ( operfunc == func_std || operfunc == func_var )
	    vars2[inp] = vars2[inp+1];
	}

      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      datetime[ndates-1].date = taxisInqVdate(taxisID1);
      datetime[ndates-1].time = taxisInqVtime(taxisID1);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  
	  streamReadRecord(streamID1, vars1[ndates-1][varID][levelID].ptr, &nmiss);
	  vars1[ndates-1][varID][levelID].nmiss = nmiss;

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  missval = vars1[ndates-1][varID][levelID].missval;

	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(vars1[ndates-1][varID][levelID].ptr[i], missval) )
	      samp1[ndates-1][varID][levelID].ptr[i] = 0;
	    else
	      {
		samp1[ndates-1][varID][levelID].ptr[i] = 1;
		for ( inp = 0; inp < ndates-1; inp++ )
		  samp1[inp][varID][levelID].ptr[i]++;
	      }

	  if ( operfunc == func_std || operfunc == func_var )
	    {
	      for ( inp = 0; inp < ndates-1; inp++ )
		{
		  farsumq(&vars2[inp][varID][levelID], vars1[ndates-1][varID][levelID]);
		  farsum(&vars1[inp][varID][levelID], vars1[ndates-1][varID][levelID]);
		}
	      farmoq(&vars2[ndates-1][varID][levelID], vars1[ndates-1][varID][levelID]);
	    }
	  else
	    {
	      for ( inp = 0; inp < ndates-1; inp++ )
		{
		  farfun(&vars1[inp][varID][levelID], vars1[ndates-1][varID][levelID], operfunc);
		}
	    }
	}

      tsID++;
    }

  for ( its = 0; its < ndates; its++ )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      free(vars1[its][varID][levelID].ptr);
	      if ( samp1[its][varID][levelID].ptr ) free(samp1[its][varID][levelID].ptr);
	      if ( operfunc == func_std || operfunc == func_var ) free(vars2[its][varID][levelID].ptr);
	    }

	  free(vars1[its][varID]);
	  free(samp1[its][varID]);
	  if ( operfunc == func_std || operfunc == func_var ) free(vars2[its][varID]);
	}

      free(vars1[its]);
      free(samp1[its]);
      if ( operfunc == func_std || operfunc == func_var ) free(vars2[its]);
    }

  free(datetime);
  free(vars1);
  free(samp1);
  if ( operfunc == func_std || operfunc == func_var ) free(vars2);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
