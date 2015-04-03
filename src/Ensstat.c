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

      Ensstat    ensmin          Ensemble minimum
      Ensstat    ensmax          Ensemble maximum
      Ensstat    enssum          Ensemble sum
      Ensstat    ensmean         Ensemble mean
      Ensstat    ensavg          Ensemble average
      Ensstat    ensstd          Ensemble standard deviation
      Ensstat    ensvar          Ensemble variance
      Ensstat    enspctl         Ensemble percentiles
*/


#if defined (_OPENMP)
#  include <omp.h>
#endif

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Ensstat(void *argument)
{
  static char func[] = "Ensstat";
  int operatorID;
  int operfunc;
  int i;
  int varID, recID;
  int gridsize = 0;
  int gridID;
  int nrecs, nrecs0;
  int levelID;
  int tsID;
  int streamID = 0, streamID2;
  int vlistID, vlistID1, vlistID2;
  int nmiss;
  int taxisID1, taxisID2;
  int ompthID;
  double missval;
  double *array2 = NULL;
  FIELD *field;
  int fileID, nfiles;
  typedef struct
  {
    int streamID;
    int vlistID;
    double *array;
  } ens_file_t;
  ens_file_t *ef = NULL;
  /* RQ */
  int pn = 0;
  /* QR */

  cdoInitialize(argument);

  cdoOperatorAdd("ensmin",  func_min,  0, NULL);
  cdoOperatorAdd("ensmax",  func_max,  0, NULL);
  cdoOperatorAdd("enssum",  func_sum,  0, NULL);
  cdoOperatorAdd("ensmean", func_mean, 0, NULL);
  cdoOperatorAdd("ensavg",  func_avg,  0, NULL);
  cdoOperatorAdd("ensvar",  func_var,  0, NULL);
  cdoOperatorAdd("ensstd",  func_std,  0, NULL);
  /* RQ */
  cdoOperatorAdd("enspctl", func_pctl,  0, NULL);
  /* QR */

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  /* RQ */
  if ( operfunc == func_pctl )
    {
      operatorInputArg("percentile number");
      pn = atoi(operatorArgv()[0]);
      
      if ( pn < 1 || pn > 99 )
        cdoAbort("Illegal argument: percentile number %d is not in the range 1..99!", pn);
    }
  /* QR */
    
  nfiles = cdoStreamCnt() - 1;

  if ( cdoVerbose )
    cdoPrint("Ensemble over %d files.", nfiles);

  ef = (ens_file_t *) malloc(nfiles*sizeof(ens_file_t));

  field = (FIELD *) malloc(ompNumThreads*sizeof(FIELD));
  for ( i = 0; i < ompNumThreads; i++ )
    {
      field[i].size   = nfiles;
      field[i].ptr    = (double *) malloc(nfiles*sizeof(double));
      field[i].weight = (double *) malloc(nfiles*sizeof(double));
      for ( fileID = 0; fileID < nfiles; fileID++ )
	field[i].weight[fileID] = 1;
    }

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID = streamOpenRead(cdoStreamName(fileID));
      if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(fileID));

      vlistID = streamInqVlist(streamID);

      ef[fileID].streamID = streamID;
      ef[fileID].vlistID = vlistID;
    }

  /* check that the contents is always the same */
  for ( fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(ef[0].vlistID, ef[fileID].vlistID, func_hrd);

  vlistID1 = ef[0].vlistID;
  vlistID2 = vlistDuplicate(vlistID1);
  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(nfiles));

  streamDefVlist(streamID2, vlistID2);
	  
  gridsize = vlistGridsizeMax(vlistID1);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    ef[fileID].array = (double *) malloc(gridsize*sizeof(double));

  array2 = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  do
    {
      nrecs0 = streamInqTimestep(ef[0].streamID, tsID);
      for ( fileID = 1; fileID < nfiles; fileID++ )
	{
	  streamID = ef[fileID].streamID;
	  nrecs = streamInqTimestep(streamID, tsID);
	  if ( nrecs != nrecs0 )
	    cdoAbort("Number of records changed from %d to %d", nrecs0, nrecs);
	}

      taxisCopyTimestep(taxisID2, taxisID1);

      if ( nrecs0 > 0 ) streamDefTimestep(streamID2, tsID);
      
      for ( recID = 0; recID < nrecs0; recID++ )
	{
	  for ( fileID = 0; fileID < nfiles; fileID++ )
	    {
	      streamID = ef[fileID].streamID;
	      streamInqRecord(streamID, &varID, &levelID);
	      streamReadRecord(streamID, ef[fileID].array, &nmiss);
	    }

	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  missval  = vlistInqVarMissval(vlistID1, varID);

	  nmiss = 0;
#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(i, ompthID, fileID)
#endif
	  for ( i = 0; i < gridsize; i++ )
	    {
#if defined (_OPENMP)
	      ompthID = omp_get_thread_num();
#else
	      ompthID = 0;
#endif
	      field[ompthID].missval = missval;
	      field[ompthID].nmiss = 0;
	      for ( fileID = 0; fileID < nfiles; fileID++ )
		{
		  field[ompthID].ptr[fileID] = ef[fileID].array[i];
		  if ( DBL_IS_EQUAL(field[ompthID].ptr[fileID], missval) )
		    field[ompthID].nmiss++;
		}

	      /* RQ */
	      if ( operfunc == func_pctl )
	        array2[i] = fldpctl(field[ompthID], pn);
	      else  
	        array2[i] = fldfun(field[ompthID], operfunc);
	      /* QR */

	      if ( DBL_IS_EQUAL(array2[i], field[ompthID].missval) ) nmiss++;
	    }

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss);
	}

      tsID++;
    }
  while ( nrecs0 > 0 );

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID = ef[fileID].streamID;
      streamClose(streamID);
    }

  streamClose(streamID2);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    if ( ef[fileID].array ) free(ef[fileID].array);

  if ( ef ) free(ef);
  if ( array2 ) free(array2);

  for ( i = 0; i < ompNumThreads; i++ )
    {
      if ( field[i].ptr    ) free(field[i].ptr);
      if ( field[i].weight ) free(field[i].weight);
    }

  if ( field ) free(field);

  cdoFinish();

  return (0);
}
