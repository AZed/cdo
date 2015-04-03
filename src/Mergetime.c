/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2007 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Merge      mergetime       Merge datasets sorted by date and time
*/


#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Mergetime(void *argument)
{
  static char func[] = "Mergetime";
  int streamID1, streamID2 = CDI_UNDEFID;
  int tsID2 = 0, recID, varID, levelID;
  int vlistID1, vlistID2;
  int streamCnt, nfiles, fileID;
  int taxisID1, taxisID2 = CDI_UNDEFID;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  int vdate, vtime;
  int next_fileID;
  double *array = NULL;
  typedef struct
  {
    int streamID;
    int vlistID;
    int taxisID;
    int tsID;
    int vdate;
    int vtime;
    int nrecs;
  } sfile_t;
  sfile_t *sf = NULL;

  cdoInitialize(argument);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  streamCnt = cdoStreamCnt();
  nfiles = streamCnt - 1;

  sf = (sfile_t *) malloc(nfiles*sizeof(sfile_t));

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID1 = streamOpenRead(cdoStreamName(fileID));
      if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(fileID));

      vlistID1 = streamInqVlist(streamID1);
      taxisID1 = vlistInqTaxis(vlistID1);

      sf[fileID].streamID = streamID1;
      sf[fileID].vlistID  = vlistID1;
      sf[fileID].taxisID  = taxisID1;
    }

  
  /* check that the contents is always the same */
  for ( fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(sf[0].vlistID, sf[fileID].vlistID, func_hrd);

  /* read the first time step */
  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      sf[fileID].tsID = 0;
      sf[fileID].nrecs = streamInqTimestep(sf[fileID].streamID, sf[fileID].tsID);
      if ( sf[fileID].nrecs == 0 )
	{
	  streamClose(sf[fileID].streamID);
	  sf[fileID].streamID = -1;
	}
      else
	{
	  sf[fileID].vdate = taxisInqVdate(sf[fileID].taxisID);
	  sf[fileID].vtime = taxisInqVtime(sf[fileID].taxisID);
	}
    }

  streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
  if ( streamID2 < 0 )
    cdiError(streamID2, "Open failed on %s", cdoStreamName(nfiles));

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(sf[0].vlistID);
      array = (double *) malloc(gridsize*sizeof(double));
    }

  while ( TRUE )
    {
      next_fileID = -1;
      vdate = 0;
      vtime = 0;
      for ( fileID = 0; fileID < nfiles; fileID++ )
	{
	  if ( sf[fileID].streamID != -1 )
	    if ( next_fileID == -1 || sf[fileID].vdate < vdate ||
		 (sf[fileID].vdate == vdate && sf[fileID].vtime < vtime) )
	      {
		next_fileID = fileID;
		vdate = sf[fileID].vdate;
		vtime = sf[fileID].vtime;
	      }
	}

      if ( cdoVerbose )
	cdoPrint("nextstep = %d  vdate = %d  vtime = %d", next_fileID, vdate, vtime);

      if ( next_fileID == -1 ) break;

      fileID = next_fileID;

      if ( tsID2 == 0 )
	{
	  vlistID1 = sf[0].vlistID;
	  vlistID2 = vlistDuplicate(vlistID1);
	  taxisID1 = vlistInqTaxis(vlistID1);
	  taxisID2 = taxisDuplicate(taxisID1);
	  vlistDefTaxis(vlistID2, taxisID2);
	      
	  streamDefVlist(streamID2, vlistID2);
	}

      taxisCopyTimestep(taxisID2, sf[fileID].taxisID);

      streamDefTimestep(streamID2, tsID2);
	       
      for ( recID = 0; recID < sf[fileID].nrecs; recID++ )
	{
	  streamInqRecord(sf[fileID].streamID, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  if ( lcopy )
	    {
	      streamCopyRecord(streamID2, sf[fileID].streamID); 
	    }
	  else
	    {
	      streamReadRecord(sf[fileID].streamID, array, &nmiss);
	      streamWriteRecord(streamID2, array, nmiss);
	    }
	}
      tsID2++;

      sf[fileID].nrecs = streamInqTimestep(sf[fileID].streamID, ++sf[fileID].tsID);
      if ( sf[fileID].nrecs == 0 )
	{
	  streamClose(sf[fileID].streamID);
	  sf[fileID].streamID = -1;
	}
      else
	{
	  sf[fileID].vdate = taxisInqVdate(sf[fileID].taxisID);
	  sf[fileID].vtime = taxisInqVtime(sf[fileID].taxisID);
	}
    }

  streamClose(streamID2);

  if ( ! lcopy )
    if ( array ) free(array);

  if ( sf ) free(sf);

  cdoFinish();

  return (0);
}
