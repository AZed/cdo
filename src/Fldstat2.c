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

      Fldstat2    fldcor         Correlation of two fields
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "functs.h"


/* routine corr copied from PINGO */
static
double corr(double * restrict in0, double * restrict in1,
	    const double * restrict weight, double missval, long gridsize)
{
  long i;
  double sum0, sum1, sum00, sum01, sum11, wsum0;
  double out;
  double missval1 = missval, missval2 = missval;

  sum0 = sum1 = sum00 = sum01 = sum11 = 0;
  wsum0 = 0;
	
  for ( i = 0; i < gridsize; ++i )
    {
      if ( weight[i] != missval &&
	   in0[i] != missval && in1[i] != missval)
	    {
	      sum0  += weight[i] * in0[i];
	      sum1  += weight[i] * in1[i];
	      sum00 += weight[i] * in0[i] * in0[i];
	      sum01 += weight[i] * in0[i] * in1[i];
	      sum11 += weight[i] * in1[i] * in1[i];
	      wsum0 += weight[i];
	    }
    }

  out = wsum0 ?
        FDIV((sum01 * wsum0 - sum0 * sum1),
	     FROOT((sum00 * wsum0 - sum0 * sum0) *
		   (sum11 * wsum0 - sum1 * sum1))) : missval;

  return (out);
}


void *Fldstat2(void *argument)
{
  static char func[] = "Fldstat2";
  int operatorID;
  int operfunc;
  int streamID1, streamID2, streamID3;
  int vlistID1, vlistID2, vlistID3;
  int gridID, lastgridID = -1;
  int gridID3;
  int wstatus = FALSE;
  int code = 0, oldcode = 0;
  int index, ngrids;
  int recID, nrecs, nrecs2;
  int tsID, varID, levelID;
  long gridsize;
  int needWeights = TRUE;
  int nmiss1, nmiss2, nmiss3;
  double missval;
  double slon, slat;
  double sglval;
  double *array1, *array2, *weight;
  int taxisID1, taxisID2, taxisID3;

  cdoInitialize(argument);

  cdoOperatorAdd("fldcor", 0, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  if ( operfunc == func_mean || operfunc == func_avg ||
       operfunc == func_var  || operfunc == func_std )
    needWeights = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, func_sft);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  slon = 0;
  slat = 0;
  gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  gridDefXvals(gridID3, &slon);
  gridDefYvals(gridID3, &slat);

  ngrids = vlistNgrids(vlistID1);

  for ( index = 0; index < ngrids; index++ )
    vlistChangeGridIndex(vlistID3, index, gridID3);

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));

  streamDefVlist(streamID3, vlistID3);

  gridsize = vlistGridsizeMax(vlistID1);

  array1 = (double *) malloc(gridsize*sizeof(double));
  array2 = (double *) malloc(gridsize*sizeof(double));
  weight = NULL;
  if ( needWeights )
    weight = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      nrecs2 = streamInqTimestep(streamID2, tsID);

      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamInqRecord(streamID2, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss1);
	  streamReadRecord(streamID2, array2, &nmiss2);

	  gridID = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  if ( needWeights && gridID != lastgridID )
	    {
	      lastgridID = gridID;
	      wstatus = gridWeights(gridID, weight);
	    }
	  code = vlistInqVarCode(vlistID1, varID);
	  if ( wstatus != 0 && tsID == 0 && code != oldcode )
	    cdoWarning("Using constant grid cell area weights for code %d!", oldcode=code);

	  missval = vlistInqVarMissval(vlistID1, varID);

	  sglval = corr(array1, array2, weight, missval, gridsize);

	  if ( DBL_IS_EQUAL(sglval, missval) )
	    nmiss3 = 1;
	  else
	    nmiss3 = 0;

	  streamDefRecord(streamID3, varID,  levelID);
	  streamWriteRecord(streamID3, &sglval, nmiss3);
	}
      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  if ( array1 ) free(array1);
  if ( array2 ) free(array2);
  if ( weight ) free(weight);

  cdoFinish();

  return (0);
}
