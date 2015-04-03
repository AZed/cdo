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

      Yearmonstat   yearmonmean        Yearly mean from monthly data
      Yearmonstat   yearmonavg         Yearly average from monthly data
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Yearmonstat(void *argument)
{
  int operatorID;
  int operfunc;
  int gridsize;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int nrecs, nrecords;
  int varID, levelID, recID;
  int tsID;
  int otsID;
  long nsets;
  double dsets;
  int i;
  int dpm;
  int year0 = 0, month0 = 0;
  int year, month, day;
  int calendar;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  char vdatestr[32], vtimestr[32];
  field_t **vars1 = NULL, **samp1 = NULL;
  field_t field;
  int timestat_date = DATE_MIDDLE;
  dtinfo_t dtinfo[13];

  cdoInitialize(argument);

  get_timestat_date(&timestat_date);

  cdoOperatorAdd("yearmonmean",  func_mean, 0, NULL);
  cdoOperatorAdd("yearmonavg",   func_avg,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);

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

  recVarID   = (int*) malloc(nrecords*sizeof(int));
  recLevelID = (int*) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);

  field_init(&field);
  field.ptr = (double*) malloc(gridsize*sizeof(double));

  vars1 = field_malloc(vlistID1, FIELD_PTR);
  samp1 = field_malloc(vlistID1, FIELD_NONE);

  tsID    = 0;
  otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      dsets = 0;
      while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
	{
	  taxisInqDTinfo(taxisID1, &dtinfo[nsets]);
	  vdate = dtinfo[nsets].v.date;
	  vtime = dtinfo[nsets].v.time;
	  cdiDecodeDate(vdate, &year, &month, &day);

	  if ( nsets == 0 ) year0 = year;

	  if ( year != year0 ) break;

	  if ( nsets > 0 && month == month0 )
	    {
	      date2str(vdate0, vdatestr, sizeof(vdatestr));
	      time2str(vtime0, vtimestr, sizeof(vtimestr));
	      cdoWarning("   last timestep: %s %s", vdatestr, vtimestr);
	      date2str(vdate, vdatestr, sizeof(vdatestr));
	      time2str(vtime, vtimestr, sizeof(vtimestr));
	      cdoWarning("current timestep: %s %s", vdatestr, vtimestr);
	      cdoAbort("Month does not change!");
	    }

	  dpm = days_per_month(calendar, year, month);

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

		  farcmul(&vars1[varID][levelID], dpm);

		  if ( nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			samp1[varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));

		      for ( i = 0; i < gridsize; i++ )
			if ( DBL_IS_EQUAL(vars1[varID][levelID].ptr[i], vars1[varID][levelID].missval) )
			  samp1[varID][levelID].ptr[i] = 0;
			else
			  samp1[varID][levelID].ptr[i] = dpm;
		    }
		}
	      else
		{
		  streamReadRecord(streamID1, field.ptr, &field.nmiss);
		  field.grid    = vars1[varID][levelID].grid;
		  field.missval = vars1[varID][levelID].missval;

		  farcmul(&field, dpm);

		  if ( field.nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			{
			  samp1[varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));
			  for ( i = 0; i < gridsize; i++ )
			    samp1[varID][levelID].ptr[i] = dsets;
			}

		      for ( i = 0; i < gridsize; i++ )
			if ( !DBL_IS_EQUAL(field.ptr[i], vars1[varID][levelID].missval) )
			  samp1[varID][levelID].ptr[i] += dpm;
		    }

		  farfun(&vars1[varID][levelID], field, operfunc);
		}
	    }

	  month0 = month;
	  vdate0 = vdate;
	  vtime0 = vtime;
	  nsets++;
	  dsets += dpm;
	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      if ( samp1[varID][levelID].ptr == NULL )
		farcmul(&vars1[varID][levelID], 1.0/dsets);
	      else
		fardiv(&vars1[varID][levelID], samp1[varID][levelID]);
	    }
	}

      if ( cdoVerbose )
	{
	  date2str(vdate0, vdatestr, sizeof(vdatestr));
	  time2str(vtime0, vtimestr, sizeof(vtimestr));
	  cdoPrint("%s %s  nsets = %d", vdatestr, vtimestr, nsets);
	}

      if      ( timestat_date == DATE_MIDDLE ) datetime_avg_dtinfo(calendar, nsets, dtinfo);
      else if ( timestat_date == DATE_FIRST  ) dtinfo[nsets].v = dtinfo[0].v;
      else if ( timestat_date == DATE_LAST   ) dtinfo[nsets].v = dtinfo[nsets-1].v;

      if ( taxisHasBounds(taxisID2) )
	{
	  dtinfo[nsets].b[0] = dtinfo[0].b[0];
	  dtinfo[nsets].b[1] = dtinfo[nsets-1].b[1];
	}

      taxisDefDTinfo(taxisID2, dtinfo[nsets]);
      streamDefTimestep(streamID2, otsID);

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  varID   = recVarID[recID];
	  levelID = recLevelID[recID];

	  if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, vars1[varID][levelID].ptr,  vars1[varID][levelID].nmiss);
	}

      if ( nrecs == 0 ) break;
      otsID++;
    }


  field_free(vars1, vlistID1);
  field_free(samp1, vlistID1);

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  cdoFinish();

  return (0);
}
