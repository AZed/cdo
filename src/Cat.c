/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2008 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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


#if  defined  (HAVE_CONFIG_H)
#  include "config.h" /* large file */
#endif

#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

void    vlistDefVarTime(int vlistID, int varID, int timeID);

void *Cat(void *argument)
{
  static char func[] = "Cat";
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int tsID1, tsID2 = 0, recID, varID, levelID;
  int vlistID1, vlistID2;
  int streamCnt, nfiles, indf;
  int taxisID1, taxisID2 = CDI_UNDEFID;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  int ntsteps;
  double *array = NULL;

  cdoInitialize(argument);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  streamCnt = cdoStreamCnt();
  nfiles = streamCnt - 1;

  for ( indf = 0; indf < nfiles; indf++ )
    {
      if ( cdoVerbose ) cdoPrint("Process file: %s", cdoStreamName(indf));

      streamID1 = streamOpenRead(cdoStreamName(indf));
      if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(indf));

      vlistID1 = streamInqVlist(streamID1);
      taxisID1 = vlistInqTaxis(vlistID1);

      if ( indf == 0 )
	{
	  int fileExist = FALSE;
	  FILE *fp;

	  fp = fopen(cdoStreamName(nfiles), "r");
	  if ( fp )
	    {
	      fclose(fp);
	      fileExist = TRUE;
	    }

	  if ( fileExist )
	    {
	      streamID2 = streamOpenAppend(cdoStreamName(nfiles));
	      if ( streamID2 < 0 )
		cdiError(streamID2, "Open failed on %s", cdoStreamName(nfiles));

	      vlistID2 = streamInqVlist(streamID2);
	      taxisID2 = vlistInqTaxis(vlistID2);

	      vlistCompare(vlistID1, vlistID2, func_sft);

	      tsID2 = vlistNtsteps(vlistID2);
	      if ( tsID2 == 0 ) tsID2 = 1; /* bug fix for constant data only */
	    }
	  else
	    {
	      if ( cdoVerbose )
		cdoPrint("Output file doesn't exist, creating: %s", cdoStreamName(nfiles));

	      streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
	      if ( streamID2 < 0 )
		cdiError(streamID2, "Open failed on %s", cdoStreamName(nfiles));

	      vlistID2 = vlistDuplicate(vlistID1);
	      taxisID2 = taxisDuplicate(taxisID1);
	      vlistDefTaxis(vlistID2, taxisID2);
	  
	      ntsteps = vlistNtsteps(vlistID1);
	      if ( ntsteps == 0 && nfiles > 1 )
		{
		  int nvars = vlistNvars(vlistID1);
		  
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
