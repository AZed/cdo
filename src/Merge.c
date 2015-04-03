/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


static
void checkDupEntry(int vlistID1, int vlistID2, const char *filename)
{
  char vname1[CDI_MAX_NAME], vname2[CDI_MAX_NAME];
  int k;
  int gridID1, gridID2;
  int zaxisID1, zaxisID2;
  int varID1, varID2;
  int nvars1, nvars2;
  int param1, param2;
  int ztype1, ztype2;
  int gtype1, gtype2;
  int nlev1,  nlev2;
  int gsize1, gsize2;
  int mlev1 = 0, mlev2 = 0;
  double *lev1 = NULL, *lev2 = NULL;

  nvars1 = vlistNvars(vlistID1);
  nvars2 = vlistNvars(vlistID2);

  for ( varID1 = 0; varID1 < nvars1; ++varID1 )
    {
      vlistInqVarName(vlistID1, varID1, vname1);
      param1   = vlistInqVarParam(vlistID1, varID1);
      gridID1  = vlistInqVarGrid(vlistID1, varID1);
      zaxisID1 = vlistInqVarZaxis(vlistID1, varID1);
      gtype1   = gridInqType(gridID1);
      gsize1   = gridInqSize(gridID1);
      ztype1   = zaxisInqType(zaxisID1);
      nlev1    = zaxisInqSize(zaxisID1);
      if ( nlev1 > mlev1 )
	{
	  mlev1 = nlev1;
	  lev1 = realloc(lev1, mlev1*sizeof(double));
	}
      zaxisInqLevels(zaxisID1, lev1);

      for ( varID2 = 0; varID2 < nvars2; ++varID2 )
	{
	  vlistInqVarName(vlistID2, varID2, vname2);
	  param2   = vlistInqVarParam(vlistID2, varID2);
	  gridID2  = vlistInqVarGrid(vlistID2, varID2);
	  zaxisID2 = vlistInqVarZaxis(vlistID2, varID2);
	  gtype2   = gridInqType(gridID2);
	  gsize2   = gridInqSize(gridID2);
	  ztype2   = zaxisInqType(zaxisID2);
	  nlev2    = zaxisInqSize(zaxisID2);
	  if ( gtype1 == gtype2 && gsize1 == gsize2 && ztype1 == ztype2 && nlev1 == nlev2 )
	    {
	      if ( nlev2 > mlev2 )
		{
		  mlev2 = nlev2;
		  lev2 = realloc(lev2, mlev2*sizeof(double));
		}
	      zaxisInqLevels(zaxisID2, lev2);

	      for ( k = 0; k < nlev2; ++k )
		if ( !IS_EQUAL(lev1[k], lev2[k]) ) break;

	      if ( k == nlev2 )
		{
		  if ( param1 < 0 || param2 < 0 )
		    {
		      if ( strcmp(vname1, vname2) == 0 )
			{
			  cdoWarning("Duplicate entry of parameter %s in %s!", vname2, filename);
			}
		    }
		  else
		    {
		      if ( param1 == param2 )
			{
			  char paramstr[32];
			  cdiParamToString(param2, paramstr, sizeof(paramstr));
			  cdoWarning("Duplicate entry of parameter %s in %s!", paramstr, filename);
			}
		    }
		}
	    }
	}
    }

  if ( lev1 ) free(lev1);
  if ( lev2 ) free(lev2);
}


void *Merge(void *argument)
{
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
  int skip_same_var = FALSE;
  const char *ofilename;
  double *array = NULL;

  cdoInitialize(argument);

  {
    char *envstr;
    envstr = getenv("SKIP_SAME_VAR");
    if ( envstr )
      {
	int ival;
	ival = atoi(envstr);
	if ( ival == 1 )
	  {
	    skip_same_var = TRUE;
	    if ( cdoVerbose )
	      cdoPrint("Set SKIP_SAME_VAR to %d", ival);
	  }
      }
  }

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  streamCnt = cdoStreamCnt();
  nmerge    = streamCnt - 1;

  ofilename = cdoStreamName(streamCnt-1)->args;

  if ( !cdoSilentMode && !cdoOverwriteMode )
    if ( fileExists(ofilename) )
      if ( !userFileOverwrite(ofilename) )
	cdoAbort("Outputfile %s already exists!", ofilename);

  streamIDs = malloc(nmerge*sizeof(int));
  vlistIDs  = malloc(nmerge*sizeof(int));

  for ( index = 0; index < nmerge; index++ )
    {
      streamID1 = streamOpenRead(cdoStreamName(index));

      streamIDs[index] = streamID1;

      vlistID1 = streamInqVlist(streamID1);
      vlistIDs[index] = vlistID1;
    }

  vlistID1 = vlistIDs[0];
  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);

  vlistID2 = vlistCreate();
  vlistCopy(vlistID2, vlistIDs[0]);
  for ( index = 1; index < nmerge; index++ )
    {
      checkDupEntry(vlistID2, vlistIDs[index], cdoStreamName(index)->args);
      /* vlistCat(vlistID2, vlistIDs[index]); */
      vlistMerge(vlistID2, vlistIDs[index]);
    }

  if ( cdoVerbose ) 
    {
      for ( index = 0; index < nmerge; index++ ) vlistPrint(vlistIDs[index]);
      vlistPrint(vlistID2);
    }
       
  streamID2 = streamOpenWrite(cdoStreamName(streamCnt-1), cdoFiletype());

  vlistDefTaxis(vlistID2, taxisID2);
  streamDefVlist(streamID2, vlistID2);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID2);
      array = malloc(gridsize*sizeof(double));
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
