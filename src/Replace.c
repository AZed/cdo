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

      Replace    replace         Replace variables
*/

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_VARS  1024

void *Replace(void *argument)
{
  static char func[] = "Replace";
  int varID;
  int varID1, nvars1;
  int varID2, nvars2;
  int nrecs = 0;
  int tsID, recID, levelID;
  int nrecs2;
  int nchvars = 0;
  int index;
  int streamID1, streamID2, streamID3;
  int vlistID1 , vlistID2, vlistID3;
  int code1 = 0, code2;
  char varname1[128], varname2[128];
  int gridsize;
  int nmiss;
  int taxisID1, taxisID3;
  int nlev, offset;
  int nts2;
  int varlist1[MAX_VARS], varlist2[MAX_VARS];
  int **varnmiss2 = NULL;
  double **vardata2 = NULL;
  double *array = NULL, *parray;

  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID3 = taxisDuplicate(taxisID1);

  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  vlistID2 = streamInqVlist(streamID2);

  /* compare all variables in vlistID2 */

  nvars1 = vlistNvars(vlistID1);
  nvars2 = vlistNvars(vlistID2);

  for ( varID2 = 0; varID2 < nvars2; varID2++ )
    {
      code2 = vlistInqVarCode(vlistID2, varID2);
      vlistInqVarName(vlistID2, varID2, varname2);

      for ( varID1 = 0; varID1 < nvars1; varID1++ )
	{
	  code1 = vlistInqVarCode(vlistID1, varID1);
	  vlistInqVarName(vlistID1, varID1, varname1);
	  if ( strcmp(varname1, varname2) == 0 ) break;
	}

      if ( code2 > 0 && varID1 == nvars1 )
	{
	  for ( varID1 = 0; varID1 < nvars1; varID1++ )
	    {
	      code1 = vlistInqVarCode(vlistID1, varID1);
	      vlistInqVarName(vlistID1, varID1, varname1);
	      if ( code1 == code2 ) break;
	    }
	}

      if ( varID1 < nvars1 )
	{
	  int gridsize1, gridsize2, nlevel1, nlevel2;

	  gridsize1 = gridInqSize(vlistInqVarGrid(vlistID1, varID1));
	  nlevel1 = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID1));

	  gridsize2 = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
	  nlevel2 = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2));

	  if ( gridsize1 != gridsize2 )
	    cdoAbort("Variables have different gridsize!");

	  if ( nlevel1 != nlevel2 )
	    cdoAbort("Variables have different number of levels!");

	  if ( cdoVerbose )
	    cdoPrint("Variable %s (code %d) replaced by  %s (code %d)",
		     varname1, code1, varname2, code2);

	  varlist1[nchvars] = varID1;
	  varlist2[nchvars] = varID2;
	  nchvars++;
	  if ( nchvars > MAX_VARS )
	    cdoAbort("Internal problem - too many variables!");
	}
      else
	{
	  cdoPrint("Variable %s (code %d) not found!", varname2, code2);
	}
    }

  if ( nchvars )
    {
      vardata2  = (double **) malloc(nchvars*sizeof(double *));
      varnmiss2 = (int **) malloc(nchvars*sizeof(int *));
      for ( varID = 0; varID < nchvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  vardata2[varID]  = (double *) malloc(nlev*gridsize*sizeof(double));
	  varnmiss2[varID] = (int *) malloc(nlev*sizeof(int));
	}
    }

  vlistID3 = vlistDuplicate(vlistID1);

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( streamID3 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(2));

  vlistDefTaxis(vlistID3, taxisID3);
  streamDefVlist(streamID3, vlistID3);

  gridsize = vlistGridsizeMax(vlistID1);
  array = (double *) malloc(gridsize*sizeof(double));

  nts2 = streamNtsteps(streamID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID3, taxisID1);

      if ( tsID == 0 || (nts2 != 0 && nts2 != 1) )
	{
	  nrecs2 = streamInqTimestep(streamID2, tsID);
	  if ( nrecs2 == 0 )
	    cdoAbort("Input streams have different number of timesteps!");

	  for ( recID = 0; recID < nrecs2; recID++ )
	    {
	      streamInqRecord(streamID2, &varID, &levelID);
	      
	      for ( index = 0; index < nchvars; index++ )
		if ( varlist2[index] == varID )
		  {
		    gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
		    offset   = gridsize*levelID;
		    parray   = vardata2[index]+offset;
		    streamReadRecord(streamID2, parray, &nmiss);
		    varnmiss2[index][levelID] = nmiss;
		    break;
		  }
	    }
	}

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  parray = array;

	  for ( index = 0; index < nchvars; index++ )
	    if ( varlist1[index] == varID )
	      {
		gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
		offset   = gridsize*levelID;
		parray   = vardata2[index]+offset;
		nmiss    = varnmiss2[index][levelID];
		break;
	      }

	  if ( index == nchvars )
	    {
	      streamReadRecord(streamID1, parray, &nmiss);
	    }

	  streamDefRecord(streamID3, varID, levelID);
	  streamWriteRecord(streamID3, parray, nmiss);
	}

      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);
 
  if ( vardata2 )
    {
      for ( varID = 0; varID < nchvars; varID++ )
	{
	  free(vardata2[varID]);
	  free(varnmiss2[varID]);
	}

      free(vardata2);
      free(varnmiss2);
    }

  if ( array ) free(array);

  cdoFinish();

  return (0);
}
