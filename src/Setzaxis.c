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

      Setzaxis   setzaxis        Set zaxis
*/

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Setzaxis(void *argument)
{
  static char func[] = "Setzaxis";
  int SETZAXIS;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int tsID, recID, varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int zaxisID1, zaxisID2 = -1;
  int nzaxis, index;
  int gridsize;
  int nmiss;
  int found;
  double *array = NULL;

  cdoInitialize(argument);

  SETZAXIS = cdoOperatorAdd("setzaxis",     0, 0, "zaxis description file");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));  

  if ( operatorID == SETZAXIS )
    {
      zaxisID2 = cdoDefineZaxis(operatorArgv()[0]);
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

  if ( operatorID == SETZAXIS )
    {
      found = 0;
      nzaxis = vlistNzaxis(vlistID1);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID1, index);

	  if ( zaxisInqSize(zaxisID1) == zaxisInqSize(zaxisID2) )
	    {
	      vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
	      found++;
	    }
	}
      if ( ! found ) cdoWarning("No zaxis with %d levels found!", zaxisInqSize(zaxisID2));
    }

  streamDefVlist(streamID2, vlistID2);

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
	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID++;
    }

  streamClose(streamID1);
  streamClose(streamID2);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}
