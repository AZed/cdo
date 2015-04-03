/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2013 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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
      Runstat    runvar1         Running variance [Divisor is (n-1)]
      Runstat    runstd          Running standard deviation
      Runstat    runstd1         Running standard deviation [Divisor is (n-1)]
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void datetime_avg_dtinfo(int calendar, int ndates, dtinfo_t *dtinfo)
{
  int vdate, vtime;
  juldate_t juldate1, juldate2, juldatem;
  double seconds;
  /*
  for ( i = 0; i < ndates; i++ )
    fprintf(stdout, "%4d %d %d\n", i+1, dtinfo[i].v.date, dtinfo[i].v.time);
  */
  if ( ndates%2 == 0 )
    {
      /*
      vdate = dtinfo[ndates-1].v.date;
      vtime = dtinfo[ndates-1].v.time;
      */
      vdate = dtinfo[ndates/2-1].v.date;
      vtime = dtinfo[ndates/2-1].v.time;
      juldate1 = juldate_encode(calendar, vdate, vtime);

      vdate = dtinfo[ndates/2].v.date;
      vtime = dtinfo[ndates/2].v.time;
      juldate2 = juldate_encode(calendar, vdate, vtime);

      seconds = juldate_to_seconds(juldate_sub(juldate2, juldate1)) / 2;
      juldatem = juldate_add_seconds(NINT(seconds), juldate1);
      juldate_decode(calendar, juldatem, &vdate, &vtime);
    }
  else
    {
      vdate = dtinfo[ndates/2].v.date;
      vtime = dtinfo[ndates/2].v.time;
    }

  dtinfo[ndates].v.date = vdate;
  dtinfo[ndates].v.time = vtime;
  /*
  fprintf(stdout, "res: %d %d\n\n", dtinfo[ndates].v.date, dtinfo[ndates].v.time);
  */
}


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


void get_timestat_date(int *tstat_date)
{
  char *envstr;

  envstr = getenv("TIMESTAT_DATE");
  if ( envstr == NULL ) envstr = getenv("RUNSTAT_DATE");
  if ( envstr )
    {
      int env_date = -1;
      char envstrl[8];

      memcpy(envstrl, envstr, 8);
      envstrl[7] = 0;
      strtolower(envstrl);

      if      ( memcmp(envstrl, "first", 5)  == 0 )  env_date = DATE_FIRST;
      else if ( memcmp(envstrl, "last", 4)   == 0 )  env_date = DATE_LAST;
      else if ( memcmp(envstrl, "middle", 6) == 0 )  env_date = DATE_MIDDLE;

      if ( env_date >= 0 )
	{
	  *tstat_date = env_date;

	  if ( cdoVerbose ) cdoPrint("Set TIMESTAT_DATE to %s", envstr);
	}
    }
}


void *Runstat(void *argument)
{
  int operatorID;
  int operfunc;
  int gridsize, gridsizemax;
  int i;
  int varID;
  int recID;
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
  int lmean = FALSE, lvarstd = FALSE, lstd = FALSE;
  int *imask;
  double missval;
  double divisor;
  field_t ***vars1 = NULL, ***vars2 = NULL, ***samp1 = NULL;
  dtinfo_t *dtinfo;
  int taxisID1, taxisID2;
  int calendar;
  int runstat_nomiss = 0;
  int timestat_date = DATE_MIDDLE;
  char *envstr;

  cdoInitialize(argument);

  envstr = getenv("RUNSTAT_NOMISS");
  if ( envstr )
    {
      char *endptr;
      int envval = (int) strtol(envstr, &endptr, 10);
      if ( envval == 1 ) runstat_nomiss = 1;
    }

  get_timestat_date(&timestat_date);

  cdoOperatorAdd("runmin",  func_min,  0, NULL);
  cdoOperatorAdd("runmax",  func_max,  0, NULL);
  cdoOperatorAdd("runsum",  func_sum,  0, NULL);
  cdoOperatorAdd("runmean", func_mean, 0, NULL);
  cdoOperatorAdd("runavg",  func_avg,  0, NULL);
  cdoOperatorAdd("runvar",  func_var,  0, NULL);
  cdoOperatorAdd("runvar1", func_var1, 0, NULL);
  cdoOperatorAdd("runstd",  func_std,  0, NULL);
  cdoOperatorAdd("runstd1", func_std1, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);

  operatorInputArg("number of timesteps");
  ndates = atoi(operatorArgv()[0]);

  lmean   = operfunc == func_mean || operfunc == func_avg;
  lstd    = operfunc == func_std || operfunc == func_std1;
  lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  divisor = operfunc == func_std1 || operfunc == func_var1;

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  calendar = taxisInqCalendar(taxisID1);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  dtinfo = (dtinfo_t *) malloc((ndates+1)*sizeof(dtinfo_t));
  vars1 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));
  if ( !runstat_nomiss )
    samp1 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));
  if ( lvarstd )
    vars2 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));

  for ( its = 0; its < ndates; its++ )
    {
      vars1[its] = field_malloc(vlistID1, FIELD_PTR);
      if ( !runstat_nomiss )
	samp1[its] = field_malloc(vlistID1, FIELD_PTR);
      if ( lvarstd )
	vars2[its] = field_malloc(vlistID1, FIELD_PTR);
    }

  gridsizemax = vlistGridsizeMax(vlistID1);
  imask = (int *) malloc(gridsizemax*sizeof(int));

  for ( tsID = 0; tsID < ndates; tsID++ )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) cdoAbort("File has less then %d timesteps!", ndates);

      taxisInqDTinfo(taxisID1, &dtinfo[tsID]);
	
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

	  if ( runstat_nomiss && nmiss > 0 ) cdoAbort("Missing values supported swichted off!");

	  if ( !runstat_nomiss )
	    {
	      gridsize = gridInqSize(vars1[0][varID][levelID].grid);
	      missval  = vars1[0][varID][levelID].missval;

	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(vars1[tsID][varID][levelID].ptr[i], missval) )
		  imask[i] = 0;
		else
		  imask[i] = 1;

	      for ( i = 0; i < gridsize; i++ )
		samp1[tsID][varID][levelID].ptr[i] = (double) imask[i];

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, inp)
#endif
	      for ( inp = 0; inp < tsID; inp++ )
		{
		  double *ptr = samp1[inp][varID][levelID].ptr;
		  for ( i = 0; i < gridsize; i++ )
		    if ( imask[i] > 0 ) ptr[i]++;
		}
	    }

	  if ( lvarstd )
	    {
	      farmoq(&vars2[tsID][varID][levelID], vars1[tsID][varID][levelID]);
#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
	      for ( inp = 0; inp < tsID; inp++ )
		{
		  farsumq(&vars2[inp][varID][levelID], vars1[tsID][varID][levelID]);
		  farsum(&vars1[inp][varID][levelID], vars1[tsID][varID][levelID]);
		}
	    }
	  else
	    {
#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
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
      if ( lmean )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( runstat_nomiss )
		  farcmul(&vars1[0][varID][levelID], 1.0/ndates);
		else
		  fardiv(&vars1[0][varID][levelID], samp1[0][varID][levelID]);
	      }
	  }
      else if ( lvarstd )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( runstat_nomiss )
		  {
		    if ( lstd )
		      farcstdx(&vars1[0][varID][levelID], vars2[0][varID][levelID], ndates, divisor);
		    else
		      farcvarx(&vars1[0][varID][levelID], vars2[0][varID][levelID], ndates, divisor);
		  }
		else
		  {
		    if ( lstd )
		      farstdx(&vars1[0][varID][levelID], vars2[0][varID][levelID], samp1[0][varID][levelID], divisor);
		    else
		      farvarx(&vars1[0][varID][levelID], vars2[0][varID][levelID], samp1[0][varID][levelID], divisor);
		  }
	      }
	  }

      if      ( timestat_date == DATE_MIDDLE ) datetime_avg_dtinfo(calendar, ndates, dtinfo);
      else if ( timestat_date == DATE_FIRST  ) dtinfo[ndates].v = dtinfo[0].v;
      else if ( timestat_date == DATE_LAST   ) dtinfo[ndates].v = dtinfo[ndates-1].v;

      if ( taxisHasBounds(taxisID2) )
	{
	  dtinfo[ndates].b[0] = dtinfo[0].b[0];
	  dtinfo[ndates].b[1] = dtinfo[ndates-1].b[1];
	}

      taxisDefDTinfo(taxisID2, dtinfo[ndates]);
      streamDefTimestep(streamID2, otsID);

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  varID    = recVarID[recID];
	  levelID  = recLevelID[recID];

	  if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, vars1[0][varID][levelID].ptr, vars1[0][varID][levelID].nmiss);
	}

      otsID++;

      dtinfo[ndates] = dtinfo[0];
      vars1[ndates] = vars1[0];
      if ( !runstat_nomiss )
	samp1[ndates] = samp1[0];
      if ( lvarstd )
        vars2[ndates] = vars2[0];

      for ( inp = 0; inp < ndates; inp++ )
	{
	  dtinfo[inp] = dtinfo[inp+1];
	  vars1[inp] = vars1[inp+1];
	  if ( !runstat_nomiss )
	    samp1[inp] = samp1[inp+1];
	  if ( lvarstd )
	    vars2[inp] = vars2[inp+1];
	}

      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      taxisInqDTinfo(taxisID1, &dtinfo[ndates-1]);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  
	  streamReadRecord(streamID1, vars1[ndates-1][varID][levelID].ptr, &nmiss);
	  vars1[ndates-1][varID][levelID].nmiss = nmiss;

	  if ( runstat_nomiss && nmiss > 0 ) cdoAbort("Missing values supported swichted off!");

	  if ( !runstat_nomiss )
	    {
	      gridsize = gridInqSize(vars1[0][varID][levelID].grid);
	      missval  = vars1[0][varID][levelID].missval;

	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(vars1[ndates-1][varID][levelID].ptr[i], missval) )
		  imask[i] = 0;
		else
		  imask[i] = 1;

	      for ( i = 0; i < gridsize; i++ )
		samp1[ndates-1][varID][levelID].ptr[i] = (double) imask[i];

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, inp)
#endif
	      for ( inp = 0; inp < ndates-1; inp++ )
		{
		  double *ptr = samp1[inp][varID][levelID].ptr;
		  for ( i = 0; i < gridsize; i++ )
		    if ( imask[i] > 0 ) ptr[i]++;
		}
	    }

	  if ( lvarstd )
	    {
	      farmoq(&vars2[ndates-1][varID][levelID], vars1[ndates-1][varID][levelID]);
#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
	      for ( inp = 0; inp < ndates-1; inp++ )
		{
		  farsumq(&vars2[inp][varID][levelID], vars1[ndates-1][varID][levelID]);
		  farsum(&vars1[inp][varID][levelID], vars1[ndates-1][varID][levelID]);
		}
	    }
	  else
	    {
#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
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
      field_free(vars1[its], vlistID1);
      if ( !runstat_nomiss ) field_free(samp1[its], vlistID1);
      if ( lvarstd ) field_free(vars2[its], vlistID1);
    }

  free(dtinfo);
  free(vars1);
  if ( !runstat_nomiss ) free(samp1);
  if ( lvarstd ) free(vars2);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);
  if ( imask )      free(imask);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
