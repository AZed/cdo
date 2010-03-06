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

      Setrtoc     setrtoc       Set range to new value
      Setrtoc2    setrtoc2      Set range to new value others to value2
*/


#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


double arg2val(char *arg)
{
  if (strcmp(arg,"inf")==0) return DBL_MAX;
  else if (strcmp(arg,"-inf")==0) return -DBL_MAX;
  else return atof(arg);
}

void *Setrange(void *argument)
{
  static char func[] = "Setrange";
  int SETRTOC, SETRTOC2;
  int operatorID;
  int streamID1, streamID2;
  int gridsize;
  int nrecs, recID;
  int tsID;
  int varID, levelID;
  int vlistID1, vlistID2;
  int nmiss;
  int i;
  double missval;
  double rmin = 0, rmax = 0;
  double *array;
  int taxisID1, taxisID2;
  double newval = 0, newval2 = 0;

  cdoInitialize(argument);

  SETRTOC = cdoOperatorAdd("setrtoc", 0, 0, "range (min, max), value");
  SETRTOC2 = cdoOperatorAdd("setrtoc2", 0, 0, "range (min, max), value1, value2");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorID == SETRTOC )
    {
      operatorCheckArgc(3);
      rmin = arg2val(operatorArgv()[0]);
      rmax = arg2val(operatorArgv()[1]);
      newval = arg2val(operatorArgv()[2]);
    }
  else if ( operatorID == SETRTOC2 )
    {
      operatorCheckArgc(4);
      rmin = arg2val(operatorArgv()[0]);
      rmax = arg2val(operatorArgv()[1]);
      newval = arg2val(operatorArgv()[2]);
      newval2 = arg2val(operatorArgv()[3]);
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
	  streamReadRecord(streamID1, array, &nmiss);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  missval = vlistInqVarMissval(vlistID1, varID);

	  if ( operatorID == SETRTOC )
	    {
	      for ( i = 0; i < gridsize; i++ )
		if ( !DBL_IS_EQUAL(array[i], missval) )
		  {
		    if ( array[i] >= rmin && array[i] <= rmax)
		      {
			array[i] = newval;
		      }
		  }
	    }
	  else if ( operatorID == SETRTOC2 )
	    {
	      for ( i = 0; i < gridsize; i++ )
		if ( !DBL_IS_EQUAL(array[i], missval) )
		  {
		    if ( array[i] >= rmin && array[i] <= rmax )
		      {
			array[i] = newval;
		      }
		    else
		      {
			array[i] = newval2;
		      }
		  }
	    }

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}