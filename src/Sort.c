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

      Sort sortcode  Sort by code number
*/


#include <string.h>
#include <stdlib.h>  /* qsort */

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


typedef struct
{
  int      recID;
  int      varID;
  int      levelID;
  int      code;
  double   level;
  char     name[128];
}
RecInfo;

static
int cmpreccode(const void *s1, const void *s2)
{
  int cmp = 0;
  RecInfo *x = (RecInfo *) s1;
  RecInfo *y = (RecInfo *) s2;
  /*
  printf("%d %d  %d %d\n", x->code, y->code, x, y);
  */
  if      ( x->code < y->code ) cmp = -1;
  else if ( x->code > y->code ) cmp =  1;

  return (cmp);
}

static
int cmpreclevel(const void *s1, const void *s2)
{
  int cmp = 0;
  RecInfo *x = (RecInfo *) s1;
  RecInfo *y = (RecInfo *) s2;
  /*
  printf("%g %g  %d %d\n", x->level, y->level, x, y);
  */
  if      ( x->level < y->level ) cmp = -1;
  else if ( x->level > y->level ) cmp =  1;

  return (cmp);
}

static
int cmprecname(const void *s1, const void *s2)
{
  RecInfo *x = (RecInfo *) s1;
  RecInfo *y = (RecInfo *) s2;

  return (strcmp(x->name, y->name));
}

static
int findrec(RecInfo *recInfo[], int nrecords, int varID, int levelID)
{
  int index;

  for ( index = 0; index < nrecords; index++ )
    if ( recInfo[index]->varID == varID && recInfo[index]->levelID == levelID )
      break;

  if ( index == nrecords )
    cdoAbort("Internal problem! Record not found.");

  return (index);
}


void *Sort(void *argument)
{
  static char func[] = "Sort";
  int SORTCODE, SORTNAME, SORTLEVEL;
  int operatorID;
  int streamID1, streamID2;
  int nrecs;
  int tsID, recID, varID, levelID, zaxisID;
  int nvars, nrecords, nlevel, offset, index;
  int vlistID1, vlistID2;
  int gridsize;
  int nmiss;
  int *recNmiss;
  double *single;
  double **vardata = NULL;
  RecInfo **recInfo;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  SORTCODE  = cdoOperatorAdd("sortcode",  0, 0, NULL);
  SORTNAME  = cdoOperatorAdd("sortname",  0, 0, NULL);
  SORTLEVEL = cdoOperatorAdd("sortlevel", 0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  /*
  if ( operatorID == SORTCODE )
      vlistSortCode(vlistID2);
   else if ( operatorID == SORTNAME )
      ;
   else if ( operatorID == SORTLEVEL )
      ;
  */

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recNmiss   = (int *) malloc(nrecords*sizeof(int));

  recInfo    = (RecInfo **) malloc(nrecords*sizeof(RecInfo *));
  recInfo[0] = (RecInfo *) malloc(nrecords*sizeof(RecInfo));

  for ( index = 1; index < nrecords; index++ )
    recInfo[index] = recInfo[0] + index;

  vardata = (double **) malloc(nvars*sizeof(double*));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      vardata[varID] = (double *) malloc(gridsize*nlevel*sizeof(double));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recInfo[recID]->recID   = recID;
	      recInfo[recID]->varID   = varID;
	      recInfo[recID]->levelID = levelID;
	      recInfo[recID]->code    = vlistInqVarCode(vlistID1, varID);
	      zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      recInfo[recID]->level   = zaxisInqLevel(zaxisID, levelID);
	      vlistInqVarName(vlistID1, varID, recInfo[recID]->name);
	    }

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single   = vardata[varID] + offset;

	  streamReadRecord(streamID1, single, &nmiss);

	  index = findrec(recInfo, nrecords, varID, levelID);
	  recNmiss[index] = nmiss;
	}

      if ( tsID == 0 )
	{
	  if      ( operatorID == SORTCODE )
	    qsort(recInfo[0], nrecords, sizeof(RecInfo), cmpreccode);
	  else if ( operatorID == SORTNAME )
	    qsort(recInfo[0], nrecords, sizeof(RecInfo), cmprecname);
	  else if ( operatorID == SORTLEVEL )
	    qsort(recInfo[0], nrecords, sizeof(RecInfo), cmpreclevel);
	}

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  /*
	  printf("recID, recID %d %d\n", recID, recInfo[recID]->recID);
	  */
	  varID   = recInfo[recID]->varID;
	  levelID = recInfo[recID]->levelID;
	  nmiss   = recNmiss[recID];

	  if ( tsID == 0 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	      offset   = gridsize*levelID;
	      single   = vardata[varID] + offset;

	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, single, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID1);
  streamClose(streamID2);

  for ( varID = 0; varID < nvars; varID++ )
    free(vardata[varID]);

  free(vardata);

  free(recInfo[0]);
  free(recInfo);

  if ( recNmiss ) free(recNmiss);

  cdoFinish();

  return (0);
}
