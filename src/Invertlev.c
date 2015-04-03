/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Invertlev     invertlev       Invert level
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "error.h"



static
void invertLevDes(int vlistID)
{
  int index, nzaxis;
  int zaxisID1, zaxisID2;
  int nlev;
  int ilev;
  int zaxistype;
  double *yv1, *yv2;
  double *ylb1, *ylb2;
  double *yub1, *yub2;

  nzaxis = vlistNzaxis(vlistID);
  for ( index = 0; index < nzaxis; index++ )
    {
      zaxisID1 = vlistZaxis(vlistID, index);
      zaxisID2 = zaxisDuplicate(zaxisID1);

      zaxistype = zaxisInqType(zaxisID1);

      nlev = zaxisInqSize(zaxisID1);

      if ( nlev < 2 || zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF ) continue;

      /* if ( zaxisInqLevels(zaxisID1, NULL) ) */
	{

	  yv1 = (double *) malloc(nlev*sizeof(double));
	  yv2 = (double *) malloc(nlev*sizeof(double));

	  zaxisInqLevels(zaxisID1, yv1);

	  for ( ilev = 0; ilev < nlev; ilev++ )
	    yv2[nlev-ilev-1] = yv1[ilev];

	  zaxisDefLevels(zaxisID2, yv2);

	  if ( yv2 ) free(yv2);
	  if ( yv1 ) free(yv1);
	}

      if ( zaxisInqLbounds(zaxisID1, NULL) && zaxisInqUbounds(zaxisID1, NULL) )
	{
	  ylb1 = (double *) malloc(nlev*sizeof(double));
	  ylb2 = (double *) malloc(nlev*sizeof(double));

	  zaxisInqLbounds(zaxisID1, ylb1);

	  yub1 = (double *) malloc(nlev*sizeof(double));
	  yub2 = (double *) malloc(nlev*sizeof(double));

	  zaxisInqUbounds(zaxisID1, yub1);

	  for ( ilev = 0; ilev < nlev; ilev++ )
	    {
	      ylb2[nlev-ilev-1] = ylb1[ilev];
	      yub2[nlev-ilev-1] = yub1[ilev];
	    }

	  zaxisDefLbounds(zaxisID2, ylb2);
	  zaxisDefUbounds(zaxisID2, yub2);

	  if ( ylb2 ) free(ylb2);
	  if ( ylb1 ) free(ylb1);
	  if ( yub2 ) free(yub2);
	  if ( yub1 ) free(yub1);
	}

      vlistChangeZaxis(vlistID, zaxisID1, zaxisID2);
    }
}


void *Invertlev(void *argument)
{
  int INVERTLEV;
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int nrecs, nvars;
  int tsID, recID, varID, levelID;
  int gridsize;
  int vlistID1, vlistID2;
  int nmiss;
  int **varnmiss;
  double *array;
  double **vardata;
  int taxisID1, taxisID2;
  int lcopy = FALSE;
  int nlev, nlevel;
  int gridID, zaxisID, zaxistype, offset;
  int linvert = FALSE;

  cdoInitialize(argument);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  INVERTLEV     = cdoOperatorAdd("invertlev",     func_all, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc   = cdoOperatorF1(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operfunc == func_all || operfunc == func_hrd )
    {
      invertLevDes(vlistID2);
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);

  array = (double *) malloc(gridsize*sizeof(double));

  nvars = vlistNvars(vlistID1);

  vardata  = (double **) malloc(nvars*sizeof(double*));
  varnmiss = (int **) malloc(nvars*sizeof(int*));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID    = vlistInqVarGrid(vlistID1, varID);
      zaxisID   = vlistInqVarZaxis(vlistID1, varID);
      gridsize  = gridInqSize(gridID);
      zaxistype = zaxisInqType(zaxisID);
      nlev      = zaxisInqSize(zaxisID);

      if ( nlev < 2 || zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF )
	{
	  vardata[varID]  = NULL;
	  varnmiss[varID] = NULL;
	}
      else
	{
	  linvert = TRUE;
	  vardata[varID]  = (double *) malloc(gridsize*nlev*sizeof(double));
	  varnmiss[varID] = (int *) malloc(nlev*sizeof(int));
	}
    }

  if ( linvert == FALSE ) cdoWarning("No variables with invertable levels found!");

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( vardata[varID] )
	    {    
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevel   = zaxisInqSize(zaxisID);
	      offset   = gridsize*levelID;

	      streamReadRecord(streamID1, vardata[varID]+offset, &nmiss);
	      varnmiss[varID][levelID] = nmiss;
	    }
	  else
	    {
	      streamDefRecord(streamID2, varID, levelID);
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
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vardata[varID] )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevel   = zaxisInqSize(zaxisID);
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  streamDefRecord(streamID2, varID, levelID);

		  offset   = gridsize*(nlevel-levelID-1);

		  nmiss = varnmiss[varID][nlevel-levelID-1];

		  streamWriteRecord(streamID2, vardata[varID]+offset, nmiss);
		}   
	    }
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vardata[varID] )
	{
	  free(varnmiss[varID]);
	  free(vardata[varID]);
	}
    }

  free(varnmiss);
  free(vardata);

  cdoFinish();

  return (0);
}
