/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.1

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

      Zonstat    zonmin          Zonal minimum
      Zonstat    zonmax          Zonal maximum
      Zonstat    zonrange        Zonal range
      Zonstat    zonsum          Zonal sum
      Zonstat    zonmean         Zonal mean
      Zonstat    zonavg          Zonal average
      Zonstat    zonstd          Zonal standard deviation
      Zonstat    zonvar          Zonal variance
      Zonstat    zonpctl         Zonal percentiles
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"


void *Zonstat(void *argument)
{
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridID1 = -1, gridID2 = -1;
  int zongridID = -1;
  int nlatmax;
  int index, ngrids;
  int recID, nrecs;
  int tsID, varID, levelID;
  int lim;
  int ndiffgrids;
  int taxisID1, taxisID2;
  field_t field1, field2;
  /* RQ */
  int pn = 0;
  /* QR */

  cdoInitialize(argument);

  cdoOperatorAdd("zonmin",   func_min,   0, NULL);
  cdoOperatorAdd("zonmax",   func_max,   0, NULL);
  cdoOperatorAdd("zonrange", func_range, 0, NULL);
  cdoOperatorAdd("zonsum",   func_sum,   0, NULL);
  cdoOperatorAdd("zonmean",  func_mean,  0, NULL);
  cdoOperatorAdd("zonavg",   func_avg,   0, NULL);
  cdoOperatorAdd("zonvar",   func_var,   0, NULL);
  cdoOperatorAdd("zonstd",   func_std,   0, NULL);
  /* RQ */
  cdoOperatorAdd("zonpctl",  func_pctl,  0, NULL);
  /* QR */

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);

  /* RQ */
  if ( operfunc == func_pctl )
    {
      operatorInputArg("percentile number");
      pn = atoi(operatorArgv()[0]);
      
      if ( pn < 1 || pn > 99 )
        cdoAbort("Illegal argument: percentile number %d is not in the range 1..99!", pn);
    }
  /* QR */

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      if ( gridInqXsize(vlistGrid(vlistID1, index)) > 1 ) 
	{
	  if ( gridID1 == -1 ) gridID1 = vlistGrid(vlistID1, index);
	}
      else
	{
	  if ( zongridID == -1 ) zongridID = vlistGrid(vlistID1, index);
	}
    }

  ndiffgrids = 0;
  for ( index = 0; index < ngrids; index++ )
    {
      if ( zongridID != -1 && zongridID == vlistGrid(vlistID1, index) ) continue;
      if ( gridID1 != vlistGrid(vlistID1, index) ) ndiffgrids++;
    }
  if ( ndiffgrids > 0 ) cdoAbort("Too many different grids!");

  if ( gridInqType(gridID1) == GRID_LONLAT   ||
       gridInqType(gridID1) == GRID_GAUSSIAN ||
       gridInqType(gridID1) == GRID_GENERIC )
    {
      if ( zongridID != -1 && gridInqYsize(zongridID) == gridInqYsize(gridID1) )
	gridID2 = zongridID;
      else
	gridID2 = gridToZonal(gridID1);
    }
  else
    {
      cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));
    }

  for ( index = 0; index < ngrids; index++ )
    vlistChangeGridIndex(vlistID2, index, gridID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridID1 = vlistInqVarGrid(vlistID1, 0);
  nlatmax = gridInqYsize(gridID1); /* max nlat ? */

  lim = vlistGridsizeMax(vlistID1);
  field_init(&field2);
  field_init(&field2);
  field1.ptr  = malloc(lim*sizeof(double));
  field2.ptr  = malloc(nlatmax*sizeof(double));
  field2.grid = gridID2;

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field1.ptr, &field1.nmiss);

	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  field1.missval = vlistInqVarMissval(vlistID1, varID);
	  field2.missval = vlistInqVarMissval(vlistID1, varID);

	  if ( zongridID != -1 && zongridID == field1.grid )
	    {
	      memcpy(field2.ptr, field1.ptr, nlatmax*sizeof(double));
	      field2.nmiss = field1.nmiss;
	    }
	  else
	    {
	      if ( operfunc == func_pctl )
		zonpctl(field1, & field2, pn);
	      else  
		zonfun(field1, &field2, operfunc);
	    }

	  streamDefRecord(streamID2, varID,  levelID);
	  streamWriteRecord(streamID2, field2.ptr, field2.nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);

  cdoFinish();

  return (0);
}
