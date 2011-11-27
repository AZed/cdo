/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2011 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Copy       cat             Concatenate datasets
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"

void    vlistDefVarTime(int vlistID, int varID, int timeID);

void *Cat(void *argument)
{
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int tsID1, tsID2 = 0, recID, varID, levelID;
  int vlistID1, vlistID2 = CDI_UNDEFID;
  int streamCnt, nfiles, indf;
  int taxisID1, taxisID2 = CDI_UNDEFID;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  int ntsteps, nvars;
  double *array = NULL;

  cdoInitialize(argument);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  streamCnt = cdoStreamCnt();
  nfiles = streamCnt - 1;

  for ( indf = 0; indf < nfiles; indf++ )
    {
      if ( cdoVerbose ) cdoPrint("Process file: %s", cdoStreamName(indf));

      streamID1 = streamOpenRead(cdoStreamName(indf));

      vlistID1 = streamInqVlist(streamID1);
      taxisID1 = vlistInqTaxis(vlistID1);

      if ( indf == 0 )
	{
	  if ( fileExist(cdoStreamName(nfiles)) )
	    {
	      streamID2 = streamOpenAppend(cdoStreamName(nfiles));

	      vlistID2 = streamInqVlist(streamID2);
	      taxisID2 = vlistInqTaxis(vlistID2);

	      vlistCompare(vlistID1, vlistID2, CMP_ALL);

	      tsID2 = vlistNtsteps(vlistID2);
	      if ( tsID2 == 0 ) tsID2 = 1; /* bug fix for time constant data only */
	    }
	  else
	    {
	      if ( cdoVerbose )
		cdoPrint("Output file doesn't exist, creating: %s", cdoStreamName(nfiles));

	      streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());

	      vlistID2 = vlistDuplicate(vlistID1);
	      taxisID2 = taxisDuplicate(taxisID1);
	      vlistDefTaxis(vlistID2, taxisID2);
	  
	      ntsteps = vlistNtsteps(vlistID1);
	      nvars   = vlistNvars(vlistID1);
	      
	      if ( ntsteps == 1 )
		{
		  for ( varID = 0; varID < nvars; ++varID )
		    if ( vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE ) break;
		  
		  if ( varID == nvars ) ntsteps = 0;
		}

	      if ( ntsteps == 0 && nfiles > 1 )
		{		  
		  for ( varID = 0; varID < nvars; ++varID )
		    vlistDefVarTime(vlistID2, varID, TIME_VARIABLE);
		}

	      streamDefVlist(streamID2, vlistID2);
	    }

	  if ( ! lcopy )
	    {
	      gridsize = vlistGridsizeMax(vlistID1);
	      array = (double *) malloc(gridsize*sizeof(double));
	    }
	}
      else
	{
	  vlistCompare(vlistID1, vlistID2, CMP_ALL);
	}

      tsID1 = 0;
      while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);

	  streamDefTimestep(streamID2, tsID2);
	       
	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);
	      streamDefRecord(streamID2,  varID,  levelID);

	      if ( lcopy )
		{
		  streamCopyRecord(streamID2, streamID1); 
		}
	      else
		{
		  streamReadRecord(streamID1, array, &nmiss);
		  streamWriteRecord(streamID2, array, nmiss);
		}
	    }
	  tsID1++;
	  tsID2++;
	}
      streamClose(streamID1);
    }

  streamClose(streamID2);
 
  if ( array ) free(array);

  cdoFinish();

  return (0);
}
