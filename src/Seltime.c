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

      Seltime    seltimestep     Select timesteps
      Seltime    seltime         Select times
      Seltime    selhour         Select hours
      Seltime    selday          Select days
      Seltime    selmon          Select months
      Seltime    selyear         Select years
      Seltime    selseas         Select seasons
      Seltime    seldate         Select dates
      Seltime    selsmon         Select single month
*/

#include <ctype.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "error.h"
#include "util.h"
#include "functs.h"
#include "list.h"


#define  NOPERATORS  32

void *Seltime(void *argument)
{
  const char func[] = "Seltime";
  int SELTIMESTEP, SELDATE, SELTIME, SELHOUR, SELDAY, SELMON, SELYEAR, SELSEAS, SELSMON;
  int operatorID;
  int operfunc, intval;
  int moddat[NOPERATORS];
  int streamID1, streamID2;
  int tsID, tsID2, nrecs;
  int recID, varID, levelID;
  int *intarr, nsel = 0, selival;
  int vlistID1 = -1, vlistID2 = -1;
  int taxisID1, taxisID2;
  int vdate, vtime;
  int copytimestep;
  int copy_nts2 = FALSE;
  int i, isel;
  int lcopy = FALSE;
  int gridsize;
  int status;
  int nmiss;
  int lnts1;
  int ncts = 0, nts, it;
  int *selfound = NULL;
  int year = 1, month = 1, day = 1, hour = 0, minute = 0, second = 0;
  int nts1 = 0, nts2 = 0;
  int its1 = 0, its2 = 0;
  double selfval = 0, *fltarr, fval = 0;
  double *array = NULL;
  LIST *ilist = listNew(INT_LIST);
  LIST *flist = listNew(FLT_LIST);
  int gridID;
  int nvars, nlevel;
  int nconst, lconstout = FALSE;
  int process_nts1 = FALSE, process_nts2 = FALSE;
  int *vdate_list = NULL, *vtime_list = NULL;
  double missval;
  double *single;
  field_t ***vars = NULL;

  cdoInitialize(argument);

  SELTIMESTEP = cdoOperatorAdd("seltimestep", func_step,     1, "timesteps");
  SELDATE     = cdoOperatorAdd("seldate",     func_datetime, 1, "start date and end date (format YYYY-MM-DDThh:mm:ss)");
  SELTIME     = cdoOperatorAdd("seltime",     func_time,     1, "times (format hh:mm:ss)");
  SELHOUR     = cdoOperatorAdd("selhour",     func_time, 10000, "hours");
  SELDAY      = cdoOperatorAdd("selday",      func_date,     1, "days");
  SELMON      = cdoOperatorAdd("selmon",      func_date,   100, "months");
  SELYEAR     = cdoOperatorAdd("selyear",     func_date, 10000, "years");
  SELSEAS     = cdoOperatorAdd("selseas",     func_date,   100, "seasons");
  SELSMON     = cdoOperatorAdd("selsmon",     func_date,   100, "month[,nts1[,nts2]]");

  moddat[SELTIMESTEP] =          1;
  /*  moddat[SELDATE]     = 1000000000; */
  moddat[SELDATE]     =          0;
  moddat[SELTIME]     =    1000000;
  moddat[SELHOUR]     =      10000;
  moddat[SELDAY]      =        100;
  moddat[SELMON]      =        100;
  moddat[SELYEAR]     = 1000000000;
  moddat[SELSEAS]     =        100;
  moddat[SELSMON]     =        100;

  operatorID = cdoOperatorID();

  operfunc = cdoOperatorFunc(operatorID);
  intval   = cdoOperatorIntval(operatorID);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorID == SELSEAS )
    {
      char Seas[3];
      int seas[4] = {FALSE, FALSE, FALSE, FALSE};
      int imon[17]; /* 1-16 ! */
      int ival;
      size_t len;
      int season_start;

      season_start = get_season_start();
      nsel = operatorArgc();
      if ( isdigit(*operatorArgv()[0]))
	for ( i = 0; i < nsel; i++ )
	  {
	    ival = atoi(operatorArgv()[i]);
	    if      ( ival == 1 || ival == 13 ) seas[0] = TRUE;
	    else if ( ival == 2 || ival == 14 ) seas[1] = TRUE;
	    else if ( ival == 3 || ival == 15 ) seas[2] = TRUE;
	    else if ( ival == 4 || ival == 16 ) seas[3] = TRUE;
	    else cdoAbort("Season %d not available!", ival);
	  }
      else
	for ( i = 0; i < nsel; i++ )
	  {
	    len = strlen(operatorArgv()[i]);
	    if ( len > 3 ) len = 3;
	    while ( len-- > 0 ) Seas[len] = toupper(operatorArgv()[i][len]);
	    if ( season_start == START_DEC )
	      {
		if      ( memcmp(Seas, "DJF", 3) == 0 ) seas[0] = TRUE;
		else if ( memcmp(Seas, "MAM", 3) == 0 ) seas[1] = TRUE;
		else if ( memcmp(Seas, "JJA", 3) == 0 ) seas[2] = TRUE;
		else if ( memcmp(Seas, "SON", 3) == 0 ) seas[3] = TRUE;
		else cdoAbort("Season %s not available!", operatorArgv()[i]);
	      }
	    else
	      {
		if      ( memcmp(Seas, "JFM", 3) == 0 ) seas[0] = TRUE;
		else if ( memcmp(Seas, "AMJ", 3) == 0 ) seas[1] = TRUE;
		else if ( memcmp(Seas, "JAS", 3) == 0 ) seas[2] = TRUE;
		else if ( memcmp(Seas, "OND", 3) == 0 ) seas[3] = TRUE;
		else cdoAbort("Season %s not available!", operatorArgv()[i]);
	      }
	  }

      for ( i = 0; i < 17; ++i ) imon[i] = 0;

      if ( season_start == START_DEC )
	{
	  if ( seas[0] ) { imon[12]++; imon[ 1]++; imon[ 2]++; imon[13]++; }
	  if ( seas[1] ) { imon[ 3]++; imon[ 4]++; imon[ 5]++; imon[14]++; }
	  if ( seas[2] ) { imon[ 6]++; imon[ 7]++; imon[ 8]++; imon[15]++; }
	  if ( seas[3] ) { imon[ 9]++; imon[10]++; imon[11]++; imon[16]++; }
	}
      else
	{
	  if ( seas[0] ) { imon[ 1]++; imon[ 2]++; imon[ 3]++; imon[13]++; }
	  if ( seas[1] ) { imon[ 4]++; imon[ 5]++; imon[ 6]++; imon[14]++; }
	  if ( seas[2] ) { imon[ 7]++; imon[ 8]++; imon[ 9]++; imon[15]++; }
	  if ( seas[3] ) { imon[10]++; imon[11]++; imon[12]++; imon[16]++; }
	}

      nsel = 0;
      for ( i = 1; i < 17; ++i )
	{
	  if ( imon[i] )
	    listSetInt(ilist, nsel++, i);
	}
    }
  else if ( operatorID == SELDATE )
    {
      int set2 = TRUE;

      nsel = operatorArgc();
      if ( nsel < 1 ) cdoAbort("Not enough arguments!");
      for ( i = 0; i < nsel; i++)
	{
	  if      ( operatorArgv()[i][0] == '-' && operatorArgv()[i][1] == 0 )
	    {
	      if ( i == 0 )
		fval = -99999999999.;
	      else
		fval =  99999999999.;

	      listSetFlt(flist, i,  fval);
	    }
	  else if ( strchr(operatorArgv()[i], '-') == NULL )
	    {
	      fval = atof(operatorArgv()[i]);
	      listSetFlt(flist, i, fval);
	    }
	  else
	    {
	      year = 1; month = 1; day = 1; hour = 0; minute = 0, second = 0;
	      if ( strchr(operatorArgv()[i], 'T') )
		{
		  status = sscanf(operatorArgv()[i], "%d-%d-%dT%d:%d:%d",
				  &year, &month, &day, &hour, &minute, &second);
		  fval = cdiEncodeTime(hour, minute, second);
		  if ( fabs(fval) > 0 ) fval /= 1000000;
		  fval += cdiEncodeDate(year, month, day);
		  listSetFlt(flist, i, fval);
		  set2 = FALSE;
		}
	      else
		{
		  status = sscanf(operatorArgv()[i], "%d-%d-%d", &year, &month, &day);
		  fval = cdiEncodeDate(year, month, day);

		  if ( nsel > 1 && i > 0 ) fval += 0.999;

		  listSetFlt(flist, i, fval);
		}
	    }
	}

      if ( nsel == 1 && set2 == TRUE )
	{
	  fval += 0.999;
	  listSetFlt(flist, nsel, fval);
	  nsel = 2;
	}
    }
  else if ( operatorID == SELTIME )
    {
      nsel = operatorArgc();
      if ( nsel < 1 ) cdoAbort("Not enough arguments!");
      for ( i = 0; i < nsel; i++ )
	{
	  if ( strchr(operatorArgv()[i], ':') )
	    {
	      sscanf(operatorArgv()[i], "%d:%d:%d", &hour, &minute, &second);
	      listSetInt(ilist, i, cdiEncodeTime(hour, minute, second));
	    }
	  else
	    {
	      listSetInt(ilist, i, atoi(operatorArgv()[i]));
	    }
	}
    }
  else
    nsel = args2intlist(operatorArgc(), operatorArgv(), ilist);

  if ( nsel < 1 ) cdoAbort("No timestep selected!");

  intarr = (int *) listArrayPtr(ilist);
  fltarr = (double *) listArrayPtr(flist);

  if ( operatorID == SELSMON )
    {
      if ( nsel > 1 ) nts1 = intarr[1];
      if ( nsel > 2 ) nts2 = intarr[2];
      else            nts2 = nts1;

      if ( nsel > 3 ) cdoAbort("Too many parameters");

      if ( cdoVerbose )
	cdoPrint("mon=%d  nts1=%d  nts2=%d", intarr[0], nts1, nts2);

      nsel = 1;
    }

  if ( nsel )
    {
      selfound = (int *) malloc(nsel*sizeof(int));
      for ( i = 0; i < nsel; i++ ) selfound[i] = FALSE;
    }

  if ( cdoVerbose )
    {
      for ( i = 0; i < nsel; i++ )
	if ( operatorID == SELDATE )
	  cdoPrint("fltarr entry: %d %14.4f", i+1, fltarr[i]);
	else
	  cdoPrint("intarr entry: %d %d", i+1, intarr[i]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);
  if ( nsel == 1 && operfunc == func_step )  vlistDefNtsteps(vlistID2, 1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double *) malloc(gridsize*sizeof(double));
    }

  nvars = vlistNvars(vlistID1);
  nconst = 0;
  for ( varID = 0; varID < nvars; varID++ )
    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) nconst++;
      
  lnts1 = operatorID == SELSMON && nts1 > 0;

  if ( lnts1 || nconst )
    {
      if ( lnts1 )
	{
	  vdate_list = (int *) malloc(nts1*sizeof(int));
	  vtime_list = (int *) malloc(nts1*sizeof(int));
	}
      else
	{
	  nts1 = 1;
	}

      vars  = (field_t ***) malloc(nts1*sizeof(field_t **));

      for ( tsID = 0; tsID < nts1; tsID++ )
	{
	  vars[tsID] = (field_t **) malloc(nvars*sizeof(field_t *));

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( lnts1 || (vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT) )
		{
		  gridID  = vlistInqVarGrid(vlistID1, varID);
		  missval = vlistInqVarMissval(vlistID1, varID);
		  nlevel  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		  gridsize = gridInqSize(gridID);
		  
		  vars[tsID][varID] = (field_t *) malloc(nlevel*sizeof(field_t));

		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    {
		      vars[tsID][varID][levelID].grid    = gridID;
		      vars[tsID][varID][levelID].missval = missval;
		      vars[tsID][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		    }
		}
	    }
	}
    }

  tsID  = 0;
  tsID2 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      copytimestep = FALSE;
      selival = -1;

      if ( operfunc == func_step )
	{
	  selival = tsID + 1;
	  if ( selival > intarr[nsel-1] ) break;
	}
      else if ( operfunc == func_date )
	{
	  selival = (vdate/intval)%moddat[operatorID];
	}
      else if ( operfunc == func_time )
	{
	  selival = (vtime/intval)%moddat[operatorID];
	}
      else if ( operfunc == func_datetime )
	{
	  selfval = vdate + vtime/1000000.;
	}

      if ( operatorID == SELDATE )
	{
	  if ( selfval >= fltarr[0] && selfval <= fltarr[nsel-1] )
	    {
	      copytimestep = TRUE;
	      selfound[0]      = TRUE;
	      selfound[nsel-1] = TRUE;
	    }
	}
      else
	{
	  for ( i = 0; i < nsel; i++ )
	    if ( selival == intarr[i] )
	      {
		copytimestep = TRUE;
		selfound[i] = TRUE;
		break;
	      }
	}

      if ( operatorID == SELSMON && copytimestep == FALSE )
	{
	  copy_nts2 = FALSE;

	  if ( process_nts1 == TRUE )
	    {
	      process_nts2 = TRUE;
	      its2 = 0;
	      process_nts1 = FALSE;
	    }

	  if ( process_nts2 == TRUE )
	    {
	      if ( its2++ < nts2 )
		{
		  copy_nts2 = TRUE;
		}
	      else
		process_nts2 = FALSE;
	    }
	}

      if ( copytimestep || copy_nts2 )
	{
	  if ( lnts1 && ncts == 0 )
	    {
	      nts = nts1;
	      if ( its1 < nts1 )
		{
		  nts = its1;
		  cdoWarning("%d timesteps missing before month %d!", nts1-its1, intarr[0]);
		}

	      for ( it = 0; it < nts; it++ )
		{
		  taxisDefVdate(taxisID2, vdate_list[it]);
		  taxisDefVtime(taxisID2, vtime_list[it]);
		  streamDefTimestep(streamID2, tsID2++);
		  
		  for ( varID = 0; varID < nvars; varID++ )
		    {
		      if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT && tsID2 > 1 ) continue;
		      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		      for ( levelID = 0; levelID < nlevel; levelID++ )
			{
			  streamDefRecord(streamID2, varID, levelID);
			  single = vars[it][varID][levelID].ptr;
			  nmiss  = vars[it][varID][levelID].nmiss;
			  streamWriteRecord(streamID2, single, nmiss);
			}
		    }
		}

	      its1 = 0;
	    }

	  ncts++;
	  if ( process_nts2 == FALSE )
	    {
	      its2 = 0;
	      process_nts1 = TRUE;
	    }

	  taxisCopyTimestep(taxisID2, taxisID1);

	  streamDefTimestep(streamID2, tsID2++);

	  if ( tsID > 0 && lconstout )
	    {
	      lconstout = FALSE;
	      nts = nts1 - 1;
	      for ( varID = 0; varID < nvars; varID++ )
		{
		  if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT )
		    {
		      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		      for ( levelID = 0; levelID < nlevel; levelID++ )
			{
			  streamDefRecord(streamID2, varID, levelID);
			  single = vars[nts][varID][levelID].ptr;
			  nmiss  = vars[nts][varID][levelID].nmiss;
			  streamWriteRecord(streamID2, single, nmiss);
			}
		    }
		}
	    }

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);
	      streamDefRecord(streamID2, varID, levelID);
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
      else
	{
	  ncts = 0;

	  if ( lnts1 || tsID == 0 )
	    {
	      if ( tsID == 0 && nconst && (!lnts1) ) lconstout = TRUE;

	      nts = nts1-1;
	      if ( lnts1 )
		{
		  if ( its1 <= nts )
		    nts = its1;
		  else
		    for ( it = 0; it < nts; it++ )
		      {
			vdate_list[it] = vdate_list[it+1];
			vtime_list[it] = vtime_list[it+1];
			for ( varID = 0; varID < nvars; varID++ )
			  {
			    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
			    gridID   = vlistInqVarGrid(vlistID1, varID);
			    gridsize = gridInqSize(gridID);
			    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
			    for ( levelID = 0; levelID < nlevel; levelID++ )
			      {
				memcpy(vars[it][varID][levelID].ptr,
				       vars[it+1][varID][levelID].ptr,
				       gridsize*sizeof(double));
				vars[it][varID][levelID].nmiss = vars[it+1][varID][levelID].nmiss;
			      }
			  }
		      }

		  vdate_list[nts] = taxisInqVdate(taxisID1);
		  vtime_list[nts] = taxisInqVtime(taxisID1);

		  its1++;
		}

	      for ( recID = 0; recID < nrecs; recID++ )
		{
		  streamInqRecord(streamID1, &varID, &levelID);
		  if ( lnts1 || (vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT) )
		    {
		      single = vars[nts][varID][levelID].ptr;
		      streamReadRecord(streamID1, single, &nmiss);
		      vars[nts][varID][levelID].nmiss = nmiss;
		    }
		}
	    }
	}
       
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);
 
  if ( operatorID == SELSMON )
    if ( its2 < nts2 )
      cdoWarning("%d timesteps missing after the last month!", nts2-its2);

  if ( ! lcopy )
    if ( array ) free(array);

  for ( isel = 0; isel < nsel; isel++ )
    {
      if ( selfound[isel] == FALSE )
	{
	  if ( operatorID == SELTIMESTEP )
	    {
	      cdoWarning("Time step %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELDATE )
	    {
	      if ( isel == 0 )
		cdoWarning("Date between %14.4f and %14.4f not found!", fltarr[0], fltarr[nsel-1]);
	    }
	  else if ( operatorID == SELTIME )
	    {
	      cdoWarning("Time %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELHOUR )
	    {
	      cdoWarning("Hour %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELDAY )
	    {
	      cdoWarning("Day %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELMON )
	    {
	      cdoWarning("Month %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELYEAR )
	    {
	      cdoWarning("Year %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELSEAS )
	    {
	      if ( isel < 3 )
		cdoWarning("Month %d not found!", intarr[isel]);
	    }
	}
    }

  if ( selfound ) free(selfound);

  listDelete(ilist);

  if ( lnts1 || nconst )
    {
      for ( tsID = 0; tsID < nts1; tsID++ )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( lnts1 || (vlistInqVarTime(vlistID2, varID) == TIME_CONSTANT) )
		{
		  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    if ( vars[tsID][varID][levelID].ptr )
		      free(vars[tsID][varID][levelID].ptr);

		  free(vars[tsID][varID]);
		}
	    }
	  free(vars[tsID]);
	}

      if ( vars  ) free(vars);
      if ( vdate_list ) free(vdate_list);
      if ( vtime_list ) free(vtime_list);
    }

  vlistDestroy(vlistID2);

  cdoFinish();

  return (NULL);
}
