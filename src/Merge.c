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

      Merge      merge           Merge datasets with different fields
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Merge(void *argument)
{
  static char func[] = "Merge";
  int streamID1 = -1, streamID2 = -1;
  int varID, varID2;
  int nrecs = 0;
  int tsID, recID, levelID, levelID2;
  int index;
  int streamCnt;
  int *streamIDs;
  int *vlistIDs;
  int vlistID1 = -1, vlistID2;
  int recID2;
  int nmerge;
  int idum = -4711;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  int taxisID1, taxisID2;
  double *array = NULL;

  cdoInitialize(argument);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  streamCnt = cdoStreamCnt();
  nmerge    = streamCnt - 1;
  streamIDs = (int *) malloc(nmerge*sizeof(int));
  vlistIDs  = (int *) malloc(nmerge*sizeof(int));

  for ( index = 0; index < nmerge; index++ )
    {
      streamID1 = streamOpenRead(cdoStreamName(index));
      if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(index));

      streamIDs[index] = streamID1;

      vlistID1 = streamInqVlist(streamID1);
      vlistIDs[index] = vlistID1;
    }

  vlistID1 = vlistIDs[0];
  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);

  vlistID2 = vlistCreate();
  vlistCopy(vlistID2, vlistIDs[0]);
  /*  for ( index = 1; index < nmerge; index++ ) vlistCat(vlistID2, vlistIDs[index]); */
  for ( index = 1; index < nmerge; index++ ) vlistMerge(vlistID2, vlistIDs[index]);

  if ( cdoVerbose ) 
    {
      for ( index = 0; index < nmerge; index++ ) vlistPrint(vlistIDs[index]);
      vlistPrint(vlistID2);
    }

  streamID2 = streamOpenWrite(cdoStreamName(streamCnt-1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(streamCnt-1));

  vlistDefTaxis(vlistID2, taxisID2);
  streamDefVlist(streamID2, vlistID2);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      array = (double *) malloc(gridsize*sizeof(double));
    }

  tsID = 0;
  while ( tsID >= 0 )
    {
      recID2 = 0;
      for ( index = 0; index < nmerge; index++ )
	{
	  streamID1 = streamIDs[index];
	  vlistID1  = vlistIDs[index];

	  if ( vlistID1 == idum ) continue;

	  nrecs = streamInqTimestep(streamID1, tsID);

	  if ( nrecs == 0 )
	    {
	      if ( tsID == 1 )
		{
		  vlistIDs[index] = idum;
		  continue;
		}
	      else
		{
		  tsID = idum;
		  break;
		}
	    }

	  if ( index == 0 )
	    {
	      taxisCopyTimestep(taxisID2, taxisID1);

	      streamDefTimestep(streamID2, tsID);
	    }

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);

	      varID2   = vlistMergedVar(vlistID1, varID);
	      levelID2 = vlistMergedLevel(vlistID1, varID, levelID);

	      if ( cdoVerbose )	cdoPrint("var %d %d %d %d", varID, levelID, varID2, levelID2);

	      streamDefRecord(streamID2, varID2, levelID2);
	      if ( lcopy )
		{
		  streamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  streamReadRecord(streamID1, array, &nmiss);
		  streamWriteRecord(streamID2, array, nmiss);
		}

	      recID2++;
	    }
	}
      tsID++;

      for ( index = 0; index < nmerge; index++ )
	if ( vlistIDs[index] != idum ) break;

      if ( index == nmerge ) tsID = idum;
    }

  for ( index = 0; index < nmerge; index++ )
    streamClose(streamIDs[index]);

  streamClose(streamID2);

  vlistDestroy(vlistID2);

  if ( streamIDs ) free(streamIDs);
  if ( vlistIDs  ) free(vlistIDs);
 
  if ( ! lcopy )
    if ( array ) free(array);

  cdoFinish();

  return (0);
}
