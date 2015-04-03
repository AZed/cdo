/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2007 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

        Timstat2        timcor      correlates two data files on the same grid
*/


#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#define  NWORK  5

void *Timstat2(void *argument)
{
  static char func[] = "Timstat2";
  int operatorID;
  int operfunc;
  int streamID1, streamID2, streamID3;
  int gridsize;
  int vdate = 0, vtime = 0;
  int nrecs, nrecs2, nrecs3, nvars, *nlevs ;
  int i;
  int tsID;
  int varID, recID, levelID, gridID;
  int nmiss1, nmiss2;
  int *recVarID, *recLevelID;
  int vlistID1, vlistID2, vlistID3;
  int taxisID1, taxisID2, taxisID3;
  int nvals2;
  double missval, missval1, missval2;
  FIELD **vars1 = NULL;
  FIELD **vars2 = NULL, ***work;
  double temp0, temp1, temp2, temp3, temp4, temp5, temp6;
  int ***nofvals;
  FIELD field1, field2;


  cdoInitialize(argument);

  cdoOperatorAdd("timcor", 2, 1, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, func_sft);
 
  gridsize = vlistGridsizeMax(vlistID1);
  nvars = vlistNvars(vlistID1);
  nrecs = vlistNrecs(vlistID1);
  nrecs3 = nrecs;
  recVarID   = (int *) malloc(nrecs*sizeof(int));
  recLevelID = (int *) malloc(nrecs*sizeof(int));
  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisID3 = taxisDuplicate(taxisID1);
 
  vlistDefTaxis(vlistID3, taxisID3);
  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));

  streamDefVlist(streamID3, vlistID3);
 
  field1.ptr = (double *) malloc(gridsize*sizeof(double));
  memset(field1.ptr, 0, gridsize*sizeof(double));
  field2.ptr = (double *) malloc(gridsize*sizeof(double));
  				 

  vars1  = (FIELD **)  malloc(nvars*sizeof(FIELD *));
  vars2  = (FIELD **)  malloc(nvars*sizeof(FIELD *));
  work   = (FIELD ***) malloc(nvars*sizeof(FIELD **));
  nofvals= (int ***)   malloc(nvars*sizeof(int **));
  nlevs  = (int *)     malloc(nvars*sizeof(int)); 

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, 0);      
      nlevs[varID]= zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval = missval1 = vlistInqVarMissval(vlistID1, varID);
      missval2 = vlistInqVarMissval(vlistID1, varID); 

      vars1[varID] = (FIELD *)   malloc(nlevs[varID]*sizeof(FIELD));
      vars2[varID] = (FIELD *)   malloc(nlevs[varID]*sizeof(FIELD));
      work[varID]  = (FIELD **)  malloc(nlevs[varID]*sizeof(FIELD*));
      nofvals[varID] = (int **)  malloc(nlevs[varID]*sizeof(int*));  

      for ( levelID = 0; levelID < nlevs[varID]; levelID++ )
	{
	  nofvals[varID][levelID] = (int *) malloc(gridsize*sizeof(int));
	  memset(nofvals[varID][levelID], 0, gridsize*sizeof(int));
      
	  vars1[varID][levelID].grid    = gridID;
	  vars1[varID][levelID].nmiss   = 0;
	  vars1[varID][levelID].missval = missval1;
	  vars1[varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));

	  vars2[varID][levelID].grid    = gridID;
	  vars2[varID][levelID].nmiss   = 0;
	  vars2[varID][levelID].missval = missval2;
	  vars2[varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));

	  work[varID][levelID] = (FIELD *) malloc(7*sizeof(FIELD));
	  for ( i = 0; i < NWORK; i++ )
	    {
	      work[varID][levelID][i].grid    = gridID;
	      work[varID][levelID][i].nmiss   = 0;
	      work[varID][levelID][i].missval = missval;
	      work[varID][levelID][i].ptr     = (double *) malloc(gridsize*sizeof(double));
	      memset(work[varID][levelID][i].ptr, 0, gridsize*sizeof(double));
	    }
	}
    }
 
  tsID=0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      nrecs2 = streamInqTimestep(streamID2, tsID);

      for ( recID = 0; recID<nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamInqRecord(streamID2, &varID, &levelID);

	  if (tsID == 0)
	    {
	      recVarID[recID] = varID;
	      recLevelID[recID] = levelID;	     	     
	    }	 

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  missval1 = vlistInqVarMissval(vlistID1, varID);
	  missval2 = vlistInqVarMissval(vlistID2, varID);

	  streamReadRecord(streamID1, field1.ptr, &field1.nmiss);
	  streamReadRecord(streamID2, field2.ptr, &field2.nmiss);
	      	     
	  for ( i = 0; i < gridsize; i++)
	    {
	      if ( ( ! DBL_IS_EQUAL(field1.ptr[i], missval1) ) && 
		   ( ! DBL_IS_EQUAL(field2.ptr[i], missval2) ) )
		{
		  work[varID][levelID][0].ptr[i] += field1.ptr[i];
		  work[varID][levelID][1].ptr[i] += field2.ptr[i];
		  work[varID][levelID][2].ptr[i] += ( field1.ptr[i]*field1.ptr[i] );
		  work[varID][levelID][3].ptr[i] += ( field2.ptr[i]*field2.ptr[i] );
		  work[varID][levelID][4].ptr[i] += ( field1.ptr[i]*field2.ptr[i] );
		  nofvals[varID][levelID][i]++;
		}
	    }	 
	}

      tsID++;
    }

  tsID = 0;
  taxisDefVdate(taxisID3, vdate);
  taxisDefVtime(taxisID3, vtime);
  streamDefTimestep(streamID3, tsID);

  for ( recID = 0; recID < nrecs3; recID++ )
    {
      varID    = recVarID[recID];
      levelID  = recLevelID[recID];     
      missval1 = vlistInqVarMissval(vlistID1, varID);
      missval2 = missval1;

      for ( i = 0; i < gridsize; i++ )
	{	  
	  nvals2 = nofvals[varID][levelID][i];
 	  if ( nvals2 > 0 )
	    {
	      temp0 = MUL(work[varID][levelID][0].ptr[i], work[varID][levelID][1].ptr[i]);
	      temp1 = SUB(work[varID][levelID][4].ptr[i], DIV(temp0, nvals2));
	      temp2 = MUL(work[varID][levelID][0].ptr[i], work[varID][levelID][0].ptr[i]);
	      temp3 = MUL(work[varID][levelID][1].ptr[i], work[varID][levelID][1].ptr[i]);
	      temp4 = SUB(work[varID][levelID][2].ptr[i], DIV(temp2, nvals2));
	      temp5 = SUB(work[varID][levelID][3].ptr[i], DIV(temp3, nvals2));
	      temp6 = MUL(temp4, temp5);

	      work[varID][levelID][0].ptr[i] = DIV(temp1, ROOT(temp6));
	      /*
		if ( work[varID][levelID][0].ptr[i] < -1)
		  work[varID][levelID][0].ptr[i] = -1;
		else if ( work[varID][levelID][0].ptr[i] > 1)
	          work[varID][levelID][0].ptr[i] = 1;
	      */
	      if ( DBL_IS_EQUAL(work[varID][levelID][0].ptr[i], missval1) ) work[varID][levelID][0].nmiss++;
	    }
	  else if ( nvals2 <= 0 ) 
	    {
	      work[varID][levelID][0].nmiss++;
	      work[varID][levelID][0].ptr[i] = missval1;
	    }
	}									       

      streamDefRecord(streamID3, varID, levelID);
      streamWriteRecord(streamID3, work[varID][levelID][0].ptr, work[varID][levelID][0].nmiss);
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      nlevs[varID] = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevs[varID]; levelID++ )
	{
	  free(vars1[varID][levelID].ptr);
	  free(vars2[varID][levelID].ptr);
	  for ( i = 0; i < NWORK; i++ )
	    free(work[varID][levelID][i].ptr);
	}
    
      free(vars1[varID]);
      free(vars2[varID]);
      free(work[varID]);
    }


  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);
    
  free(vars1);
  free(vars2);
  free(work);
  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);
  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);
    
  cdoFinish();   
 
  return (0);
}
