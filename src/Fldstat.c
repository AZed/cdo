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

      Fldstat    fldmin          Field minimum
      Fldstat    fldmax          Field maximum
      Fldstat    fldsum          Field sum
      Fldstat    fldmean         Field mean
      Fldstat    fldavg          Field average
      Fldstat    fldstd          Field standard deviation
      Fldstat    fldvar          Field variance
      Fldstat    fldpctl         Field percentiles
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "functs.h"


void *Fldstat(void *argument)
{
  static char func[] = "Fldstat";
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridID2, lastgrid = -1;
  int wstatus = FALSE;
  int code = 0, oldcode = 0;
  int index, ngrids;
  int recID, nrecs;
  int tsID, varID, levelID;
  int lim;
  int needWeights = FALSE;
  int nmiss;
  double slon, slat;
  double sglval;
  FIELD field;
  int taxisID1, taxisID2;
  /* RQ */
  int pn = 0;
  /* QR */

  cdoInitialize(argument);

  cdoOperatorAdd("fldmin",  func_min,  0, NULL);
  cdoOperatorAdd("fldmax",  func_max,  0, NULL);
  cdoOperatorAdd("fldsum",  func_sum,  0, NULL);
  cdoOperatorAdd("fldmean", func_mean, 0, NULL);
  cdoOperatorAdd("fldavg",  func_avg,  0, NULL);
  cdoOperatorAdd("fldvar",  func_var,  0, NULL);
  cdoOperatorAdd("fldstd",  func_std,  0, NULL);
  /* RQ */
  cdoOperatorAdd("fldpctl", func_pctl,  0, NULL);
  /* QR */

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  /* RQ */
  if ( operfunc == func_pctl )
    {
      operatorInputArg("percentile number");
      pn = atoi(operatorArgv()[0]);
      
      if ( pn < 1 || pn > 99 )
        cdoAbort("Illegal argument: percentile number %d is not in the range 1..99!", pn);
    }
  /* QR */

  if ( operfunc == func_mean || operfunc == func_avg ||
       operfunc == func_var  || operfunc == func_std )
    needWeights = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  slon = 0;
  slat = 0;
  gridID2 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  gridDefXvals(gridID2, &slon);
  gridDefYvals(gridID2, &slat);

  ngrids = vlistNgrids(vlistID1);

  for ( index = 0; index < ngrids; index++ )
    vlistChangeGridIndex(vlistID2, index, gridID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  lim = vlistGridsizeMax(vlistID1);
  field.ptr    = (double *) malloc(lim*sizeof(double));
  field.weight = NULL;
  if ( needWeights )
    field.weight = (double *) malloc(lim*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field.ptr, &field.nmiss);

	  field.grid = vlistInqVarGrid(vlistID1, varID);
	  field.size = gridInqSize(field.grid);
	  if ( needWeights && field.grid != lastgrid )
	    {
	      lastgrid = field.grid;
	      wstatus = gridWeights(field.grid, field.weight);
	    }
	  code = vlistInqVarCode(vlistID1, varID);
	  if ( wstatus != 0 && tsID == 0 && code != oldcode )
	    cdoWarning("Using constant grid cell area weights for code %d!", oldcode=code);

	  field.missval = vlistInqVarMissval(vlistID1, varID);

	  /* RQ */
	  if ( operfunc == func_pctl )
	    sglval = fldpctl(field, pn);
	  else  
	    sglval = fldfun(field, operfunc);
	  /* QR */

	  if ( cdoVerbose )
	    if ( operfunc == func_min || operfunc == func_max )
	      {
		if ( gridInqType(field.grid) == GRID_GAUSSIAN ||
		     gridInqType(field.grid) == GRID_LONLAT )
		  {
		    int i = 0, j, nlon, nlat;
		    nlon = gridInqXsize(field.grid);
		    nlat = gridInqYsize(field.grid);
		    for ( j = 0; j < nlat; ++j )
		      {
			for ( i = 0; i < nlon; ++i )
			  {
			    if ( DBL_IS_EQUAL(field.ptr[j*nlon+i], sglval) ) break;
			  }
			if ( i < nlon ) break;
		      }

		    if ( j < nlat )
		      {
			int vdate, vtime, code;
			int year, month, day, hour, minute, second;
			double level;
			double xval, yval;
			xval = gridInqXval(field.grid, i);
			yval = gridInqYval(field.grid, j);
			vdate = taxisInqVdate(taxisID1);
			vtime = taxisInqVtime(taxisID1);
			decode_date(vdate, &year, &month, &day);
			decode_time(vtime, &hour, &minute, &second);
			code = vlistInqVarCode(vlistID1, varID);
			level = zaxisInqLevel(vlistInqVarZaxis(vlistID1, varID), levelID);
			if ( tsID == 0 && recID == 0 )
			  {
			    if ( operfunc == func_min )
			      fprintf(stdout, "  Date     Time     Code  Level   Lon      Lat          Minval\n");
			    else
			      fprintf(stdout, "  Date     Time     Code  Level   Lon      Lat          Maxval\n");
			  }

			fprintf(stdout, "%4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d %3d %7g %9.7g %9.7g %12.5g\n",
				year, month, day, hour, minute, second,
				code, level, xval, yval, sglval);
		      }
		  }
	      }

	  if ( DBL_IS_EQUAL(sglval, field.missval) )
	    nmiss = 1;
	  else
	    nmiss = 0;

	  streamDefRecord(streamID2, varID,  levelID);
	  streamWriteRecord(streamID2, &sglval, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field.ptr )    free(field.ptr);
  if ( field.weight ) free(field.weight);

  cdoFinish();

  return (0);
}
