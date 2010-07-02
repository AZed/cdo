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

      Setgrid    setgrid         Set grid
      Setgrid    setgridtype     Set grid type
*/

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void *Setgrid(void *argument)
{
  static char func[] = "Setgrid";
  int SETGRID, SETGRIDTYPE, SETGRIDAREA;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int tsID, recID, varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int gridID1, gridID2 = -1;
  int ngrids, index;
  int gridsize, gridtype = -1;
  int nmiss;
  int found;
  int areasize = 0;
  int lregular = 0;
  char *gridname = NULL;
  double *areaweight = NULL;
  double *array = NULL;

  cdoInitialize(argument);

  SETGRID     = cdoOperatorAdd("setgrid",     0, 0, "grid description file or name");
  SETGRIDTYPE = cdoOperatorAdd("setgridtype", 0, 0, "grid type");
  SETGRIDAREA = cdoOperatorAdd("setgridarea", 0, 0, "filename with area weights");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));  

  if ( operatorID == SETGRID )
    {
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
    }
  else if ( operatorID == SETGRIDTYPE )
    {
      gridname = operatorArgv()[0];

      if      ( strcmp(gridname, "curvilinear") == 0 ) gridtype = GRID_CURVILINEAR;
      else if ( strcmp(gridname, "cell") == 0 )        gridtype = GRID_CELL;
      else if ( strcmp(gridname, "lonlat") == 0 )      gridtype = GRID_LONLAT;
      else if ( strcmp(gridname, "gaussian") == 0 )    gridtype = GRID_GAUSSIAN;
      else if ( strcmp(gridname, "regular") == 0 )    {gridtype = GRID_GAUSSIAN; lregular = 1;}
      else cdoAbort("Unsupported grid name: %s", gridname);
    }
  else if ( operatorID == SETGRIDAREA )
    {
      int streamID, vlistID, gridID;
      char *areafile;

      areafile = operatorArgv()[0];
      streamID = streamOpenRead(areafile);
      if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", areafile);

      vlistID = streamInqVlist(streamID);

      nrecs = streamInqTimestep(streamID, 0);
      streamInqRecord(streamID, &varID, &levelID);

      gridID = vlistInqVarGrid(vlistID, varID);
      areasize = gridInqSize(gridID);
      areaweight = (double *) malloc(areasize*sizeof(double));
  
      streamReadRecord(streamID, areaweight, &nmiss);

      streamClose(streamID);

      if ( cdoVerbose )
	{
	  int i;
	  double arrmean, arrmin, arrmax;

	  arrmean = areaweight[0];
	  arrmin  = areaweight[0];
	  arrmax  = areaweight[0];
	  for ( i = 1; i < areasize; i++ )
	    {
	      if ( areaweight[i] < arrmin ) arrmin = areaweight[i];
	      if ( areaweight[i] > arrmax ) arrmax = areaweight[i];
	      arrmean += areaweight[i];
	    }
	  arrmean = arrmean/areasize;

	  cdoPrint("areaweights: %d %#12.5g%#12.5g%#12.5g", areasize, arrmin, arrmean, arrmax);
	}
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  if ( operatorID == SETGRID )
    {
      found = 0;
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1 = vlistGrid(vlistID1, index);

	  if ( gridInqSize(gridID1) == gridInqSize(gridID2) )
	    {
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	      found++;
	    }
	}
      if ( ! found ) cdoWarning("No grid with %d points found!", gridInqSize(gridID2));
    }
  else if ( operatorID == SETGRIDTYPE )
    {
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1 = vlistGrid(vlistID1, index);

	  if ( lregular )
	    {
	      if ( gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED )
		{
		  gridID2 = gridToRegular(gridID1);
		}
	    }
	  else
	    {
	      if      ( gridtype == GRID_CURVILINEAR ) gridID2 = gridToCurvilinear(gridID1);
	      else if ( gridtype == GRID_CELL )        gridID2 = gridToCell(gridID1);
	      else cdoAbort("Unsupported grid name: %s", gridname);
	    }

	  /*	  gridCompress(gridID2); */
	  vlistChangeGridIndex(vlistID2, index, gridID2);
	}
    }
  else if ( operatorID == SETGRIDAREA )
    {
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1  = vlistGrid(vlistID1, index);
	  gridtype = gridInqType(gridID1);
	  gridsize = gridInqSize(gridID1);
	  if ( gridsize == areasize )
	    {
	      gridID2 = gridDuplicate(gridID1);
	      gridDefArea(gridID2, areaweight);
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	    }
	}
    }

  streamDefVlist(streamID2, vlistID2);
  //vlistPrint(vlistID2);

  if ( lregular )
    gridsize = vlistGridsizeMax(vlistID2);
  else
    gridsize = vlistGridsizeMax(vlistID1);

  array = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  streamReadRecord(streamID1, array, &nmiss);
	  if ( lregular )
	    {
	      gridID1 = vlistInqVarGrid(vlistID1, varID);
	      gridID2 = vlistInqVarGrid(vlistID2, varID);
	      if ( gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED )
		{
		  double missval = vlistInqVarMissval(vlistID1, varID);
		  field2regular(gridID1, gridID2, missval, array, nmiss);
		}
	    }

	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( areaweight ) free(areaweight);
  if ( array ) free(array);

  cdoFinish();

  return (0);
}
