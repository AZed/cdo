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

      Split      splitcode       Split codes
      Split      splitname       Split variables
      Split      splitlevel      Split levels
      Split      splitgrid       Split grids
      Split      splitzaxis      Split zaxis
*/


#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Split(void *argument)
{
  static char func[] = "Split";
  int SPLITCODE, SPLITNAME, SPLITLEVEL, SPLITGRID, SPLITZAXIS;
  int operatorID;
  int nchars;
  int streamID1;
  int varID;
  int code;
  int nrecs, nvars, nzaxis, nlevs;
  int tsID, recID, levelID, zaxisID, levID;
  int varID2, levelID2;
  int vlistID1, vlistID2;
  int  *vlistIDs = NULL, *streamIDs = NULL;
  int  itmp[999];
  double ftmp[999];
  char filesuffix[32];
  char filename[1024];
  int nsplit = 0;
  int index;
  int i;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  double *array = NULL;

  cdoInitialize(argument);

  SPLITCODE  = cdoOperatorAdd("splitcode",  0, 0, NULL);
  SPLITNAME  = cdoOperatorAdd("splitname",  0, 0, NULL);
  SPLITLEVEL = cdoOperatorAdd("splitlevel", 0, 0, NULL);
  SPLITGRID  = cdoOperatorAdd("splitgrid",  0, 0, NULL);
  SPLITZAXIS = cdoOperatorAdd("splitzaxis", 0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  nvars  = vlistNvars(vlistID1);
  nrecs  = vlistNrecs(vlistID1);
  nzaxis = vlistNzaxis(vlistID1);

  strcpy(filename, cdoStreamName(1));
  nchars = strlen(filename);

  filesuffix[0] = 0;
  if ( cdoDisableFilesuffix == FALSE )
    {
      strcat(filesuffix, streamFilesuffix(cdoDefaultFileType));
      if ( cdoDefaultFileType == FILETYPE_GRB )
	if ( vlistIsSzipped(vlistID1) || cdoZtype == COMPRESS_SZIP )
	  strcat(filesuffix, ".sz");
    }

  if ( operatorID == SPLITCODE )
    {
      int *codes = NULL;
      nsplit = 0;
      for ( varID = 0; varID < nvars; varID++ )
	{
	  code = vlistInqVarCode(vlistID1, varID);
	  for ( index = 0; index < varID; index++ )
	    if ( code == vlistInqVarCode(vlistID1, index) ) break;

	  if ( index == varID )
	    {
	      itmp[nsplit] = code;
	      nsplit++;
	    }
	}

      codes     = (int *) malloc(nsplit*sizeof(int));
      vlistIDs  = (int *) malloc(nsplit*sizeof(int));
      streamIDs = (int *) malloc(nsplit*sizeof(int));
      memcpy(codes, itmp, nsplit*sizeof(int));

      for ( index = 0; index < nsplit; index++ )
	{
	  vlistClearFlag(vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      code    = vlistInqVarCode(vlistID1, varID);
	      zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      nlevs   = zaxisInqSize(zaxisID);
	      if ( codes[index] == code )
		{
		  for ( levID = 0; levID < nlevs; levID++ )
		    {
		      vlistDefIndex(vlistID1, varID, levID, index);
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	    }
	  vlistID2 = vlistCreate();
	  vlistCopyFlag(vlistID2, vlistID1);
	  vlistIDs[index] = vlistID2;

	  sprintf(filename+nchars, "%03d", codes[index]);
	  if ( filesuffix[0] )
	    sprintf(filename+nchars+3, "%s", filesuffix);
	  streamIDs[index] = streamOpenWrite(filename, cdoFiletype());
	  if ( streamIDs[index] < 0 ) cdiError(streamIDs[index], "Open failed on %s", filename);

	  streamDefVlist(streamIDs[index], vlistIDs[index]);
	}
      if ( codes ) free(codes);
    }
  else if ( operatorID == SPLITNAME )
    {
      char varname[128];
      nsplit = nvars;

      vlistIDs  = (int *) malloc(nsplit*sizeof(int));
      streamIDs = (int *) malloc(nsplit*sizeof(int));

      for ( index = 0; index < nsplit; index++ )
	{
	  vlistClearFlag(vlistID1);
	  varID = index;
	  zaxisID = vlistInqVarZaxis(vlistID1, varID);
	  nlevs   = zaxisInqSize(zaxisID);
	  for ( levID = 0; levID < nlevs; levID++ )
	    {
	      vlistDefIndex(vlistID1, varID, levID, index);
	      vlistDefFlag(vlistID1, varID, levID, TRUE);
	    }

	  vlistID2 = vlistCreate();
	  vlistCopyFlag(vlistID2, vlistID1);
	  vlistIDs[index] = vlistID2;

	  filename[nchars] = '\0';
	  vlistInqVarName(vlistID1, varID, varname);
	  strcat(filename, varname);
	  if ( filesuffix[0] )
	    strcat(filename, filesuffix);
	  streamIDs[index] = streamOpenWrite(filename, cdoFiletype());
	  if ( streamIDs[index] < 0 ) cdiError(streamIDs[index], "Open failed on %s", filename);

	  streamDefVlist(streamIDs[index], vlistID2);
	}
    }
  else if ( operatorID == SPLITLEVEL )
    {
      double level, *levels = NULL;
      nzaxis = vlistNzaxis(vlistID1);
      nsplit = 0;
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID = vlistZaxis(vlistID1, index);
	  nlevs   = zaxisInqSize(zaxisID);
	  for ( levID = 0; levID < nlevs; levID++ )
	    {
	      level = zaxisInqLevel(zaxisID, levID);
	      for ( i = 0; i < nsplit; i++ )
		if ( IS_EQUAL(level, ftmp[i]) ) break;
	      if ( i == nsplit )
		ftmp[nsplit++] = level;
	    }
	}

      levels    = (double *) malloc(nsplit*sizeof(double));
      vlistIDs  = (int *) malloc(nsplit*sizeof(int));
      streamIDs = (int *) malloc(nsplit*sizeof(int));
      memcpy(levels, ftmp, nsplit*sizeof(double));

      for ( index = 0; index < nsplit; index++ )
	{
	  vlistClearFlag(vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      nlevs   = zaxisInqSize(zaxisID);
	      for ( levID = 0; levID < nlevs; levID++ )
		{
		  level = zaxisInqLevel(zaxisID, levID);
		  if ( IS_EQUAL(levels[index], level) )
		    {
		      vlistDefIndex(vlistID1, varID, levID, index);
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	    }
	  vlistID2 = vlistCreate();
	  vlistCopyFlag(vlistID2, vlistID1);
	  vlistIDs[index] = vlistID2;

	  sprintf(filename+nchars, "%06g", levels[index]);
	  if ( filesuffix[0] )
	    sprintf(filename+nchars+6, "%s", filesuffix);
	  streamIDs[index] = streamOpenWrite(filename, cdoFiletype());
	  if ( streamIDs[index] < 0 ) cdiError(streamIDs[index], "Open failed on %s", filename);

	  streamDefVlist(streamIDs[index], vlistID2);
	}
      if ( levels ) free(levels);
    }
  else if ( operatorID == SPLITGRID )
    {
      int gridID, *gridIDs = NULL;

      nsplit = vlistNgrids(vlistID1);

      gridIDs   = (int *) malloc(nsplit*sizeof(int));
      vlistIDs  = (int *) malloc(nsplit*sizeof(int));
      streamIDs = (int *) malloc(nsplit*sizeof(int));

      for ( index = 0; index < nsplit; index++ )
	gridIDs[index] = vlistGrid(vlistID1, index);

      for ( index = 0; index < nsplit; index++ )
	{
	  vlistClearFlag(vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID  = vlistInqVarGrid(vlistID1, varID);
	      zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      nlevs   = zaxisInqSize(zaxisID);
	      if ( gridIDs[index] == gridID )
		{
		  for ( levID = 0; levID < nlevs; levID++ )
		    {
		      vlistDefIndex(vlistID1, varID, levID, index);
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	    }
	  vlistID2 = vlistCreate();
	  vlistCopyFlag(vlistID2, vlistID1);
	  vlistIDs[index] = vlistID2;

	  sprintf(filename+nchars, "%02d", gridIDs[index]+1);
	  if ( filesuffix[0] )
	    sprintf(filename+nchars+2, "%s", filesuffix);
	  streamIDs[index] = streamOpenWrite(filename, cdoFiletype());
	  if ( streamIDs[index] < 0 ) cdiError(streamIDs[index], "Open failed on %s", filename);

	  streamDefVlist(streamIDs[index], vlistID2);
	}
      if ( gridIDs ) free(gridIDs);
    }
  else if ( operatorID == SPLITZAXIS )
    {
      int zaxisID, *zaxisIDs = NULL;

      nsplit = vlistNzaxis(vlistID1);

      zaxisIDs  = (int *) malloc(nsplit*sizeof(int));
      vlistIDs  = (int *) malloc(nsplit*sizeof(int));
      streamIDs = (int *) malloc(nsplit*sizeof(int));

      for ( index = 0; index < nsplit; index++ )
	zaxisIDs[index] = vlistZaxis(vlistID1, index);

      for ( index = 0; index < nsplit; index++ )
	{
	  vlistClearFlag(vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      nlevs   = zaxisInqSize(zaxisID);
	      if ( zaxisIDs[index] == zaxisID )
		{
		  for ( levID = 0; levID < nlevs; levID++ )
		    {
		      vlistDefIndex(vlistID1, varID, levID, index);
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	    }
	  vlistID2 = vlistCreate();
	  vlistCopyFlag(vlistID2, vlistID1);
	  vlistIDs[index] = vlistID2;

	  sprintf(filename+nchars, "%02d", zaxisIDs[index]+1);
	  if ( filesuffix[0] )
	    sprintf(filename+nchars+2, "%s", filesuffix);
	  streamIDs[index] = streamOpenWrite(filename, cdoFiletype());
	  if ( streamIDs[index] < 0 ) cdiError(streamIDs[index], "Open failed on %s", filename);

	  streamDefVlist(streamIDs[index], vlistID2);
	}
      if ( zaxisIDs ) free(zaxisIDs);
    }
  else
    {
      cdoAbort("not implemented!");
    }

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      array = (double *) malloc(gridsize*sizeof(double));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      for ( index = 0; index < nsplit; index++ )
	streamDefTimestep(streamIDs[index], tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  index    = vlistInqIndex(vlistID1, varID, levelID);
	  vlistID2 = vlistIDs[index];
	  varID2   = vlistFindVar(vlistID2, varID);
	  levelID2 = vlistFindLevel(vlistID2, varID, levelID);
	  /*
	    printf("%d %d %d %d %d %d\n", index, vlistID2, varID, levelID, varID2, levelID2);
	  */

	  streamDefRecord(streamIDs[index], varID2, levelID2);
	  if ( lcopy )
	    {
	      streamCopyRecord(streamIDs[index], streamID1);
	    }
	  else
	    {
	      streamReadRecord(streamID1, array, &nmiss);
	      streamWriteRecord(streamIDs[index], array, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID1);

  for ( index = 0; index < nsplit; index++ )
    {
      vlistDestroy(vlistIDs[index]);
      streamClose(streamIDs[index]);
    }
 
  if ( ! lcopy )
    if ( array ) free(array);

  if ( vlistIDs ) free(vlistIDs);
  if ( streamIDs ) free(streamIDs);

  cdoFinish();

  return (0);
}
