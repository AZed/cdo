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

      Vardup     pardup          Duplicate parameters
      Vardup     parmul          Multiply parameters
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Vardup(void *argument)
{
  int PARDUP, PARMUL;
  int operatorID;
  int streamID1, streamID2;
  int nrecs, nrecords;
  int tsID, recID, varID, varID2, levelID;
  int gridsize, i;
  int vlistID1, vlistID2;
  long offset;
  int nmul = 0;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  int **varnmiss;
  double *single;
  double **vardata;
  double *array;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  PARDUP = cdoOperatorAdd("pardup", 0, 0, NULL);
  PARMUL = cdoOperatorAdd("parmul", 0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == PARDUP )
    {
      nmul = 2;
    }
  else if ( operatorID == PARMUL )
    {
      operatorInputArg("number of multiply");
      nmul = atoi(operatorArgv()[0]);
    }
  else
    cdoAbort("operator not implemented!");

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int*) malloc(nrecords*sizeof(int));
  recLevelID = (int*) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);
  array    = (double*) malloc(gridsize*sizeof(double));
  vardata  = (double **) malloc(nvars*sizeof(double *));
  varnmiss = (int **) malloc(nvars*sizeof(int *));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      vardata[varID]  = (double*) malloc(gridsize*nlevel*sizeof(double));
      varnmiss[varID] = (int*) malloc(nlevel*sizeof(int));
    }

  for ( i = 1; i < nmul; i++ )
    {
      vlistCat(vlistID2, vlistID1);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarCode(vlistID2, varID+nvars*i, -(varID+nvars*i+1));
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  recVarID[recID]   = varID;
	  recLevelID[recID] = levelID;

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single   = vardata[varID] + offset;
  
	  streamReadRecord(streamID1, single, &nmiss);
	  varnmiss[varID][levelID] = nmiss;
	}

      for ( i = 0; i < nmul; i++ )
	for ( recID = 0; recID < nrecs; recID++ )
	  {
	    varID    = recVarID[recID];
	    varID2   = varID + i*nvars;
	    levelID  = recLevelID[recID];

	    gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	    offset   = gridsize*levelID;
	    single   = vardata[varID] + offset;
	    nmiss    = varnmiss[varID][levelID];

	    memcpy(array, single, gridsize*sizeof(double));

	    streamDefRecord(streamID2,  varID2,  levelID);
	    streamWriteRecord(streamID2, array, nmiss);
	  }

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  for ( varID = 0; varID < nvars; varID++ )  free(vardata[varID]);
  for ( varID = 0; varID < nvars; varID++ )  free(varnmiss[varID]);
  free(vardata);
  free(varnmiss);
  free(array);

  cdoFinish();

  return (0);
}
