/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

      Yseasstat  yseasmin        Multi-year seasonally minimum
      Yseasstat  yseasmax        Multi-year seasonally maximum
      Yseasstat  yseassum        Multi-year seasonally sum
      Yseasstat  yseasmean       Multi-year seasonally mean
      Yseasstat  yseasavg        Multi-year seasonally average
      Yseasstat  yseasvar        Multi-year seasonally variance
      Yseasstat  yseasstd        Multi-year seasonally standard deviation
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


#define  NSEAS       4

typedef struct {
  int vdate;
  int vtime;
}
date_time_t;


static 
void set_date(int vdate_new, int vtime_new, date_time_t *datetime)
{
  int year, month, day;

  cdiDecodeDate(vdate_new, &year, &month, &day);
  if ( month == 12 ) vdate_new = cdiEncodeDate(year-1, month, day);

  if ( vdate_new > datetime->vdate )
    {
      datetime->vdate = vdate_new;
      datetime->vtime = vtime_new;
    }
}


void *Yseasstat(void *argument)
{
  int operatorID;
  int operfunc;
  int gridsize;
  int i;
  int varID;
  int recID;
  int vdate, vtime;
  int year, month, day, seas;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  long nsets[NSEAS];
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  date_time_t datetime[NSEAS];
  field_t **vars1[NSEAS], **vars2[NSEAS], **samp1[NSEAS];
  field_t field;
  int season_start;

  cdoInitialize(argument);

  cdoOperatorAdd("yseasmin",  func_min,  0, NULL);
  cdoOperatorAdd("yseasmax",  func_max,  0, NULL);
  cdoOperatorAdd("yseassum",  func_sum,  0, NULL);
  cdoOperatorAdd("yseasmean", func_mean, 0, NULL);
  cdoOperatorAdd("yseasavg",  func_avg,  0, NULL);
  cdoOperatorAdd("yseasvar",  func_var,  0, NULL);
  cdoOperatorAdd("yseasstd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);

  season_start = get_season_start();
  for ( seas = 0; seas < NSEAS; seas++ )
    {
      vars1[seas]  = NULL;
      vars2[seas]  = NULL;
      samp1[seas]  = NULL;
      nsets[seas]  = 0;
      datetime[seas].vdate = 0;
      datetime[seas].vtime = 0;
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID2) ) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int*) malloc(nrecords*sizeof(int));
  recLevelID = (int*) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);
  field_init(&field);
  field.ptr = (double*) malloc(gridsize*sizeof(double));

  tsID = 0;
  otsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);
      cdiDecodeDate(vdate, &year, &month, &day);
      if ( month < 0 || month > 16 )
	cdoAbort("Month %d out of range!", month);

      if ( season_start == START_DEC )
	{
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

      set_date(vdate, vtime, &datetime[seas]);

      if ( vars1[seas] == NULL )
	{
	  vars1[seas] = field_malloc(vlistID1, FIELD_PTR);
	  samp1[seas] = field_malloc(vlistID1, FIELD_NONE);
	  if ( operfunc == func_std || operfunc == func_var )
	    vars2[seas] = field_malloc(vlistID1, FIELD_PTR);
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

	  if ( nsets[seas] == 0 )
	    {
	      streamReadRecord(streamID1, vars1[seas][varID][levelID].ptr, &nmiss);
	      vars1[seas][varID][levelID].nmiss = nmiss;

	      if ( nmiss > 0 || samp1[seas][varID][levelID].ptr )
		{
		  if ( samp1[seas][varID][levelID].ptr == NULL )
		    samp1[seas][varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));

		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(vars1[seas][varID][levelID].ptr[i],
				      vars1[seas][varID][levelID].missval) )
		      samp1[seas][varID][levelID].ptr[i] = 0;
		    else
		      samp1[seas][varID][levelID].ptr[i] = 1;
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &field.nmiss);
	      field.grid    = vars1[seas][varID][levelID].grid;
	      field.missval = vars1[seas][varID][levelID].missval;

	      if ( field.nmiss > 0 || samp1[seas][varID][levelID].ptr )
		{
		  if ( samp1[seas][varID][levelID].ptr == NULL )
		    {
		      samp1[seas][varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));
		      for ( i = 0; i < gridsize; i++ )
			samp1[seas][varID][levelID].ptr[i] = nsets[seas];
		    }
		  
		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], vars1[seas][varID][levelID].missval) )
		      samp1[seas][varID][levelID].ptr[i]++;
		}

	      if ( operfunc == func_std || operfunc == func_var )
		{
		  farsumq(&vars2[seas][varID][levelID], field);
		  farsum(&vars1[seas][varID][levelID], field);
		}
	      else
		{
		  farfun(&vars1[seas][varID][levelID], field, operfunc);
		}
	    }
	}

      if ( nsets[seas] == 0 && (operfunc == func_std || operfunc == func_var) )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      farmoq(&vars2[seas][varID][levelID], vars1[seas][varID][levelID]);
	  }

      nsets[seas]++;
      tsID++;
    }

  for ( seas = 0; seas < NSEAS; seas++ )
    if ( nsets[seas] )
      {
	if ( operfunc == func_mean || operfunc == func_avg )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( samp1[seas][varID][levelID].ptr == NULL )
		    farcmul(&vars1[seas][varID][levelID], 1.0/nsets[seas]);
		  else
		    fardiv(&vars1[seas][varID][levelID], samp1[seas][varID][levelID]);
		}
	    }
	else if ( operfunc == func_std || operfunc == func_var )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( samp1[seas][varID][levelID].ptr == NULL )
		    {
		      if ( operfunc == func_std )
			farcstd(&vars1[seas][varID][levelID], vars2[seas][varID][levelID], 1.0/nsets[seas]);
		      else
			farcvar(&vars1[seas][varID][levelID], vars2[seas][varID][levelID], 1.0/nsets[seas]);
		    }
		  else
		    {
		      farinv(&samp1[seas][varID][levelID]);
		      if ( operfunc == func_std )
			farstd(&vars1[seas][varID][levelID], vars2[seas][varID][levelID], samp1[seas][varID][levelID]);
		      else
			farvar(&vars1[seas][varID][levelID], vars2[seas][varID][levelID], samp1[seas][varID][levelID]);
		    }
		}
	    }

	taxisDefVdate(taxisID2, datetime[seas].vdate);
	taxisDefVtime(taxisID2, datetime[seas].vtime);
	streamDefTimestep(streamID2, otsID);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	    streamDefRecord(streamID2, varID, levelID);
	    streamWriteRecord(streamID2, vars1[seas][varID][levelID].ptr,
			      vars1[seas][varID][levelID].nmiss);
	  }

	otsID++;
      }

  for ( seas = 0; seas < NSEAS; seas++ )
    {
      if ( vars1[seas] != NULL )
	{
	  field_free(vars1[seas], vlistID1);
	  field_free(samp1[seas], vlistID1);
	  if ( operfunc == func_std || operfunc == func_var ) free(vars2[seas]);
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
