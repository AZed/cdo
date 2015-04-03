/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2009 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Detrend    detrend         Detrend
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  NALLOC_INC  1000


static
void detrend(int nts, double missval1, double *array1, double *array2)
{
  int j;
  int n;
  double sumj, sumjj;
  double sumx, sumjx;
  double work1, work2;
  double missval2 = missval1;

  sumx = sumjx = 0;
  sumj = sumjj = n = 0;
  for ( j = 0; j < nts; j++ )
    if ( !DBL_IS_EQUAL(array1[j], missval1) )
      {
	sumx  += array1[j];
	sumjx += j * array1[j];
	sumj  += j;
	sumjj += j * j;
	n++;
      }

  work1 = DIV(SUB(sumjx, DIV(MUL(sumx, sumj), n) ),
	      SUB(sumjj, DIV(MUL(sumj, sumj), n)) );
  work2 = SUB(DIV(sumx, n), MUL(work1, DIV(sumj, n)));

  for ( j = 0; j < nts; j++ )
    array2[j] = SUB(array1[j], ADD(work2, MUL(j, work1)));

}


void *Detrend(void *argument)
{
  static char func[] = "Detrend";
  int gridsize;
  int nrecs;
  int gridID, varID, levelID, recID;
  int tsID;
  int i;
  int nts;
  int nalloc = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *vdate = NULL, *vtime = NULL;
  double missval;
  double *array1, *array2;
  FIELD ***vars = NULL;

  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  nvars = vlistNvars(vlistID1);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID >= nalloc )
	{
	  nalloc += NALLOC_INC;
	  vdate = (int *) realloc(vdate, nalloc*sizeof(int));
	  vtime = (int *) realloc(vtime, nalloc*sizeof(int));
	  vars  = (FIELD ***) realloc(vars, nalloc*sizeof(FIELD **));
	}

      vdate[tsID] = taxisInqVdate(taxisID1);
      vtime[tsID] = taxisInqVtime(taxisID1);

      vars[tsID] = (FIELD **) malloc(nvars*sizeof(FIELD *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  missval  = vlistInqVarMissval(vlistID1, varID);
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

	  vars[tsID][varID] = (FIELD *) malloc(nlevel*sizeof(FIELD));

	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      vars[tsID][varID][levelID].grid    = gridID;
	      vars[tsID][varID][levelID].missval = missval;
	      vars[tsID][varID][levelID].ptr     = NULL;
	    }
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  vars[tsID][varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));
	  streamReadRecord(streamID1, vars[tsID][varID][levelID].ptr, &nmiss);
	  vars[tsID][varID][levelID].nmiss = nmiss;
	}

      tsID++;
    }

  nts = tsID;

  array1 = (double *) malloc(nts*sizeof(double));
  array2 = (double *) malloc(nts*sizeof(double));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      missval  = vlistInqVarMissval(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  for ( i = 0; i < gridsize; i++ )
	    {
	      for ( tsID = 0; tsID < nts; tsID++ )
		array1[tsID] = vars[tsID][varID][levelID].ptr[i];

	      detrend(nts, missval, array1, array2);

	      for ( tsID = 0; tsID < nts; tsID++ )
		vars[tsID][varID][levelID].ptr[i] = array2[tsID];
	    }
	}
    }

  if ( array1 ) free(array1);
  if ( array2 ) free(array2);

  for ( tsID = 0; tsID < nts; tsID++ )
    {
      taxisDefVdate(taxisID2, vdate[tsID]);
      taxisDefVtime(taxisID2, vtime[tsID]);
      streamDefTimestep(streamID2, tsID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      if ( vars[tsID][varID][levelID].ptr )
		{
		  nmiss = vars[tsID][varID][levelID].nmiss;
		  streamDefRecord(streamID2, varID, levelID);
		  streamWriteRecord(streamID2, vars[tsID][varID][levelID].ptr, nmiss);
		  free(vars[tsID][varID][levelID].ptr);
		}
	    }
	  free(vars[tsID][varID]);
	}
      free(vars[tsID]);
    }

  if ( vars  ) free(vars);
  if ( vdate ) free(vdate);
  if ( vtime ) free(vtime);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
