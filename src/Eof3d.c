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

     EOF3d        eof3d             3D-EOF in spatial or time space
     EOF3d        eof3dspatial      3D-EOF in spatial space
     EOF3d        eof3dtime         3D-EOF in time space
*/
/*
 * TODO: 
 * Role of the weights for eofs. Should not be mixed up with division with
 * number of contributing values during summation.
 */

enum T_EIGEN_MODE {JACOBI, DANIELSON_LANCZOS};

#define WEIGHTS 1

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "statistic.h"
#if defined (_OPENMP)
#include "omp.h"
#endif

// NO MISSING VALUE SUPPORT ADDED SO FAR

void *EOF3d(void * argument)
{
  char *envstr;

  enum {EOF3D_, EOF3D_TIME, EOF3D_SPATIAL};

  int **datacounts;
  int gridsize, temp_size = 0;
  int gridID1, gridID2, gridID3;
  int i, i2, j, j1, j2, eofID, varID, recID, levelID, tsID;
  int *miss;
  int missval_warning=0;
  int nmiss,ngrids,n_eig,nrecs,nvars,n=0,nlevs=0,npack=0,nts=0;
  int offset;
  int operatorID, operfunc;
  int *pack;
  int reached_eof;
  int streamID1, streamID2, streamID3;
  int taxisID1, taxisID2, taxisID3;
  int timer_init = 0, timer_alloc = 0, timer_read = 0, timer_cov = 0;
  int timer_eig = 0, timer_post = 0, timer_write = 0, timer_finish = 0;
  int *varID2;
  int vdate=0, vtime=0;
  int vlistID1, vlistID2=-1, vlistID3=-1;
  int zaxisID2;

  int calendar = CALENDAR_STANDARD;
  juldate_t juldate;

  double missval=0;
  double sum_w, sum;
  double **cov = NULL;                                /* TODO: covariance matrix / eigenvectors after solving */
  double *eigv;
  double *weight;
  double *xvals, *yvals, *zvals;
  double *df1p, *df2p;

  field_t **datafields;
  field_t **eigenvectors, **eigenvalues;
  field_t in;

  enum T_EIGEN_MODE eigen_mode = JACOBI;


  if ( cdoTimer ) {
    timer_init = timer_new("Timeof init");
    timer_alloc= timer_new("Timeof alloc");
    timer_read = timer_new("Timeof read");
    timer_cov  = timer_new("Timeof cov");
    timer_eig  = timer_new("Timeof eig");
    timer_post = timer_new("Timeof post");
    timer_write= timer_new("Timeof write");
    timer_finish=timer_new("Timeof finish");

    timer_start(timer_init);
  }

  cdoInitialize(argument);
  cdoOperatorAdd("eof3d",       EOF3D_,       0, NULL);
  cdoOperatorAdd("eof3dtime",   EOF3D_TIME,   0, NULL);
  cdoOperatorAdd("eof3dspatial",EOF3D_SPATIAL,0, NULL);

  operatorID  = cdoOperatorID();
  operfunc    = cdoOperatorF1(operatorID);

  operatorInputArg("Number of eigen functions to write out");
  n_eig       = atoi(operatorArgv()[0]);

  envstr = getenv("CDO_SVD_MODE");

  if ( envstr &&! strncmp(envstr,"danielson_lanczos",17) )
    eigen_mode = DANIELSON_LANCZOS;
  else if ( envstr && ! strncmp(envstr,"jacobi",6) )
    eigen_mode = JACOBI;
  else if ( envstr ) {
    cdoWarning("Unknown environmental setting %s for CDO_SVD_MODE. Available options are",envstr);
    cdoWarning("  - 'jacobi' for a one-sided parallelized jacobi algorithm");
    cdoWarning("  - 'danielson_lanzcos' for the D/L algorithm");
  }

  if ( cdoVerbose ) 
    cdoPrint("Set eigen_mode to %s\n",eigen_mode == JACOBI? "jacobi" : "danielson_lanczos");

#if defined (_OPENMP)
  if ( omp_get_max_threads() > 1 && eigen_mode == DANIELSON_LANCZOS )  {
    cdoWarning("Requested parallel computation with %i Threads ",omp_get_max_threads());
    cdoWarning("  but environmental setting CDO_SVD_MODE causes sequential ");
    cdoWarning("  Singular value decomposition");
  }
#endif 

  streamID1   = streamOpenRead(cdoStreamName(0));
  vlistID1    = streamInqVlist(streamID1);
  taxisID1    = vlistInqTaxis(vlistID1);
  gridID1     = vlistInqVarGrid(vlistID1, 0);
  gridsize    = vlistGridsizeMax(vlistID1);
  nvars       = vlistNvars(vlistID1);
  nrecs       = vlistNrecs(vlistID1);
  taxisID1    = vlistInqTaxis(vlistID1);
  weight      = (double*) malloc(gridsize*sizeof(double));
  if ( WEIGHTS )
      gridWeights(gridID1, &weight[0]);
  else
    for(i=0;i<gridsize;i++)
      weight[i]=1;


  /*  eigenvalues */
  streamID2   = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  taxisID2    = taxisDuplicate(taxisID1);

  gridID2     = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  xvals       = (double*) malloc(1*sizeof(double));
  yvals       = (double*) malloc(1*sizeof(double));
  zvals       = (double*) malloc(1*sizeof(double));
  xvals[0]    = 0;
  yvals[0]    = 0;
  zvals[0]    = 0;
  gridDefXvals(gridID2, xvals);
  gridDefYvals(gridID2, yvals);

  zaxisID2 = zaxisCreate(ZAXIS_GENERIC,1);
  zaxisDefLevels(zaxisID2,zvals);
  zaxisDefName(zaxisID2,"zaxis_Reduced");
  zaxisDefLongname(zaxisID2,"Reduced zaxis from EOF3D - only one eigen value per 3D eigen vector");

  vlistID2 = vlistCreate();
  taxisDefRdate(taxisID2, 0);
  taxisDefRtime(taxisID2, 0);
  vlistDefTaxis(vlistID2, taxisID2);

  varID2 = (int *) malloc (nvars*sizeof(int));
  for ( varID=0; varID<nvars; varID++ )
    varID2[varID] = vlistDefVar(vlistID2, gridID2, zaxisID2, TIME_VARIABLE);
  ngrids      = vlistNgrids(vlistID2);
  for ( i = 0; i < ngrids; i++ )
    vlistChangeGridIndex(vlistID2, i, gridID2);

  /*  eigenvectors */
  streamID3   = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  vlistID3    = vlistDuplicate(vlistID1);
  taxisID3    = taxisDuplicate(taxisID1);
  gridID3     = gridDuplicate(gridID1);
  taxisDefRdate(taxisID3, 0);
  taxisDefRtime(taxisID3, 0);
  vlistDefTaxis(vlistID3, taxisID3);
  

  if ( cdoVerbose )
    cdoPrint("Initialized streams");

  /*  eigenvalues */

  reached_eof = 0;
  tsID        = 0;

  if ( operfunc == EOF3D_SPATIAL )
    cdoAbort("Operator not Implemented - use eof3d or eof3dtime instead");


  /* COUNT NUMBER OF TIMESTEPS if EOF3D_ or EOF3D_TIME */
  while ( TRUE )
    {
      if ( reached_eof ) continue;
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) {
	reached_eof = 1;
	break;
      }
      tsID++;
    }

  nts         = tsID;
  reached_eof = 0;
  streamID1   = streamOpenRead(cdoStreamName(0));

  /* reset the requested number of eigen-function to the maximum if neccessary */
  if ( n_eig > nts )
    {
      cdoWarning("Solving in time-space:");
      cdoWarning("Number of eigen-functions to write out is bigger than number of time-steps.");
      cdoWarning("Setting n_eig to %i.", nts);
      n_eig = nts;
    }
  n = nts;

  if ( cdoVerbose ) 
    cdoPrint("counted %i timesteps",n);

  if ( cdoTimer ) timer_stop(timer_init);
  if ( cdoTimer ) timer_start(timer_alloc);

  /* allocation of temporary fields and output structures */
  in.ptr       = (double *)   malloc(gridsize*sizeof(double));
  datafields   = (field_t **) malloc(nvars*sizeof(field_t*));
  datacounts   = (int     **) malloc(nvars*sizeof(int *));
  eigenvectors = (field_t **) malloc(nvars*sizeof(field_t*));
  eigenvalues  = (field_t **) malloc(nvars*sizeof(field_t*));

  for ( varID = 0; varID < nvars; ++varID )
    {
      gridID1             = vlistInqVarGrid(vlistID1, varID);
      gridsize            = vlistGridsizeMax(vlistID1);
      nlevs               = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      temp_size           = gridsize * nlevs;
      missval             = vlistInqVarMissval(vlistID1, varID);

      datafields[varID]   = (field_t *) malloc(nlevs*sizeof(field_t*));
      datacounts[varID]   = (int *)     malloc(nlevs*sizeof(int* ));
      eigenvectors[varID] = (field_t *) malloc(nlevs*sizeof(field_t*));

      datafields[varID] = (field_t *) malloc(nts*sizeof(field_t));
      for ( tsID = 0; tsID < nts; tsID++ )
	{
	  datafields[varID][tsID].grid    = gridID1;
	  datafields[varID][tsID].nmiss   = 0;
	  datafields[varID][tsID].missval = missval;
	  datafields[varID][tsID].ptr     = (double *) malloc(temp_size*sizeof(double));
	  for ( i = 0; i < temp_size; ++i )
	    datafields[varID][tsID].ptr[i] = 0;
	}
      datacounts[varID] = (int *) malloc(temp_size*sizeof(int));	      
      for(i=0;i<temp_size;i++)
	datacounts[varID][i] = 0;
      
      eigenvectors[varID] = (field_t *) malloc(n_eig*sizeof(field_t));
      eigenvalues[varID]  = (field_t *) malloc(nts*sizeof(field_t));

      for ( i = 0; i < n; i++ )
	{
	  if ( i < n_eig )
	    {
	      eigenvectors[varID][i].grid    = gridID2;
	      eigenvectors[varID][i].nmiss   = 0;
	      eigenvectors[varID][i].missval = missval;
	      eigenvectors[varID][i].ptr     = (double *) malloc(temp_size*sizeof(double));
	      for ( i2 = 0; i2 < temp_size; ++i2 )
		eigenvectors[varID][i].ptr[i2] = missval;
	    }
	  
	  eigenvalues[varID][i].grid    = gridID3;
	  eigenvalues[varID][i].nmiss   = 0;
	  eigenvalues[varID][i].missval = missval;
	  eigenvalues[varID][i].ptr     = (double *) malloc(1*sizeof(double));
	  eigenvalues[varID][i].ptr[0]  = missval;
	}
    }

  if ( cdoVerbose)
    cdoPrint("allocated eigenvalue/eigenvector with nts=%i, n=%i, gridsize=%i for processing in %s",
	     nts,n,gridsize,"time_space");
  
  if ( cdoTimer ) timer_stop(timer_alloc);

  if ( cdoTimer ) timer_start(timer_read);
  tsID = 0;

  /* read the data and create covariance matrices for each var & level */
  while ( TRUE )
    {
      if ( reached_eof ) continue;
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 )
        {
          reached_eof = 1;
          break;
        }

      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID1, &varID, &levelID);
          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

          missval  = in.missval = vlistInqVarMissval(vlistID1, varID);
          streamReadRecord(streamID1, in.ptr, &in.nmiss);

	  offset = gridsize * levelID;
	  for ( i=0; i<gridsize; ++i )
	    {
	      if ( ! DBL_IS_EQUAL(in.ptr[i], missval ) )
		{
		  datafields[varID][tsID].ptr[offset + i] = in.ptr[i];
		  datacounts[varID][offset + i]++;
		}
	      else
		{
		  if ( missval_warning == 0 )
		    {
		      cdoWarning("Missing Value Support not Checked for this Operator!");
		      cdoWarning("Does not work with changing locations of missing values in time.");
		      missval_warning = 1;
		    }
		  datafields[varID][tsID].ptr[i+offset] = 0;
		}
	    }
        }
      tsID++;
    }

  if ( cdoVerbose ) 
    cdoPrint("Read data for %i variables",nvars);
  
  pack = (int *) malloc(temp_size*sizeof(int)); //TODO
  miss = (int *) malloc(temp_size*sizeof(int));

  if ( cdoTimer ) timer_stop(timer_read);

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevs               = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      temp_size = gridsize * nlevs;

      if ( cdoVerbose )  {
	char vname[64];
	vlistInqVarName(vlistID1,varID,&vname[0]);
	cdoPrint("============================================================================");
	cdoPrint("Calculating covariance matrix and SVD for var%i (%s)",varID,vname);
      }

      npack        = 0;    // TODO already set to 0
      sum_w        = 0;

      if ( cdoTimer ) timer_start(timer_cov);
      
      sum_w = 0;
      for ( i = 0; i < temp_size ; i++ )
	{
	  if ( datacounts[varID][i] > 1)
	    {
	      pack[npack] = i;
	      npack++;
	      sum_w += weight[i%gridsize];
	    }
	}

      if ( npack < 1 ) {
	char vname[64];
	vlistInqVarName(vlistID1,varID,&vname[0]);
	cdoWarning("Refusing to calculate EOF from a single time step for var%i (%s)",varID,&vname[0]);
	continue;
      }

	  
      cov = (double **) malloc (nts*sizeof(double*));
      for ( j1 = 0; j1 < nts; j1++)
	cov[j1] = (double*) malloc(nts*sizeof(double));
      eigv = (double *) malloc(n*sizeof(double));

      if ( cdoVerbose )  {
	cdoPrint("varID %i allocated eigv and cov with nts=%i and n=%i",varID,nts,n);
	cdoPrint("   npack=%i, nts=%i temp_size=%i",npack,nts,temp_size);
      }


#if defined (_OPENMP)
#pragma omp parallel for private(j1,j2,sum,df1p,df2p) default(shared) schedule(static,2000)
#endif 
      for ( j1 = 0; j1 < nts; j1++)
	for ( j2 = j1; j2 < nts; j2++ )
	  {
	    sum = 0;
	    df1p = datafields[varID][j1].ptr;
	    df2p = datafields[varID][j2].ptr;
	    for ( i = 0; i < npack; i++ )
	      sum += weight[pack[i]%gridsize]*df1p[pack[i]]*df2p[pack[i]];
	    cov[j2][j1] = cov[j1][j2] = sum / sum_w / nts;
	  }
      if ( cdoVerbose ) 
	cdoPrint("calculated cov-matrix");

      /* SOLVE THE EIGEN PROBLEM */
      if ( cdoTimer ) timer_stop(timer_cov);


      if ( cdoTimer ) timer_start(timer_eig);

      if ( cdoVerbose ) 
	cdoPrint("Processed correlation matrix for var %2i | npack: %4i",varID,n);

      if ( eigen_mode == JACOBI ) 
	parallel_eigen_solution_of_symmetric_matrix(&cov[0],&eigv[0],n,n,__func__);
      else 
	eigen_solution_of_symmetric_matrix(&cov[0],&eigv[0],n,n,__func__);
      /* NOW: cov contains the eigenvectors, eigv the eigenvalues */

      if ( cdoVerbose ) 
	cdoPrint("Processed SVD decomposition for var %i from %i x %i matrix",varID,n,n);

      for( eofID=0; eofID<n; eofID++ )
	eigenvalues[varID][eofID].ptr[0] = eigv[eofID];
      
      if ( cdoTimer ) timer_stop(timer_eig);
      if ( cdoTimer ) timer_start(timer_post);

      for ( eofID = 0; eofID < n_eig; eofID++ )
	{
#if defined (_OPENMP)
#pragma omp parallel for private(i,j,sum) shared(datafields, eigenvectors)
#endif 
	  for ( i = 0; i < npack; i++ )
	    {
	      sum = 0;
	      for ( j = 0; j < nts; j++ )
		sum += datafields[varID][j].ptr[pack[i]] * cov[eofID][j];
	      eigenvectors[varID][eofID].ptr[pack[i]] = sum;
	    }
	  // NORMALIZING
	  sum = 0;

#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) reduction(+:sum) \
  shared(eigenvectors,weight,pack,varID,eofID,npack,gridsize)
#endif 
	  for ( i = 0; i < npack; i++ )
	    sum +=  weight[pack[i]%gridsize] *
	      eigenvectors[varID][eofID].ptr[pack[i]] *
	      eigenvectors[varID][eofID].ptr[pack[i]];

	  if ( sum > 0 ) {
	    sum = sqrt(sum);
#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) \
  shared(sum,npack,eigenvectors,varID,eofID,pack)
#endif
	    for( i = 0; i < npack; i++ )
	      eigenvectors[varID][eofID].ptr[pack[i]] /= sum;
	  }
	  else
#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) \
  shared(eigenvectors,varID,eofID,pack,missval,npack)
#endif
	    for( i = 0; i < npack; i++ )
	      eigenvectors[varID][eofID].ptr[pack[i]] = missval;
	}     /* for ( eofID = 0; eofID < n_eig; eofID++ )     */

      if ( cdoTimer ) timer_stop(timer_post);

      if ( eigv ) free(eigv);
      for ( i=0; i<n; i++ )
	if ( cov[i] ) 
	  free(cov[i]);
    }         /* for ( varID = 0; varID < nvars; varID++ )    */

  /* write files with eigenvalues (ID3) and eigenvectors (ID2) */


  if ( cdoTimer ) timer_start(timer_write);

  cdoPrint("Started writing");
  streamDefVlist(streamID2, vlistID2);
  streamDefVlist(streamID3, vlistID3);

  vdate = 10101;
  vtime = 0;
  juldate = juldate_encode(calendar, vdate, vtime);
  for ( tsID = 0; tsID < n; tsID++ )
    {
      juldate = juldate_add_seconds(60, juldate);
      juldate_decode(calendar, juldate, &vdate, &vtime);

      taxisDefVdate(taxisID2, vdate);
      taxisDefVtime(taxisID2, vtime);
      streamDefTimestep(streamID2, tsID);

      if ( tsID < n_eig )
        {
          taxisDefVdate(taxisID3, vdate);
          taxisDefVtime(taxisID3, vtime);
          streamDefTimestep(streamID3, tsID);
        }

      for ( varID = 0; varID < nvars; varID++ )
        {
          nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for ( levelID = 0; levelID < nlevs; levelID++ )
            {
	      offset = levelID * gridsize;
              if ( tsID < n_eig )
                {
                  nmiss = 0;
                  for ( i = 0; i < gridsize; i++ )
                    if ( DBL_IS_EQUAL(eigenvectors[varID][tsID].ptr[offset + i], missval) ) nmiss++;

                  streamDefRecord(streamID3, varID, levelID);
                  streamWriteRecord(streamID3, &eigenvectors[varID][tsID].ptr[offset], nmiss);
                }
	    }
	  if ( DBL_IS_EQUAL(eigenvalues[varID][tsID].ptr[i], missval) ) nmiss = 1;
	  else nmiss = 0;
	  streamDefRecord(streamID2, varID, 0);
	  streamWriteRecord(streamID2, eigenvalues[varID][tsID].ptr,nmiss);
        } // for ( varID = 0; ... )
    } // for ( tsID = 0; ... )

  if ( cdoTimer ) timer_stop(timer_write);
  
  if ( cdoTimer ) timer_start(timer_finish);

  cdoPrint("Started cleanup in eof3d");
  
  for ( varID = 0; varID < nvars; varID++)
    {
      for(i = 0; i < nts; i++)
	{
	  if ( i < n_eig )
	    free(eigenvectors[varID][i].ptr);
	  free(eigenvalues[varID][i].ptr);
	}
      free(eigenvectors[varID]);
      free(eigenvalues[varID]);
      free(datacounts[varID]);
    }

  free(eigenvectors);
  free(eigenvalues);
  free(datafields);
  free(datacounts);
  free(in.ptr);

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  if ( cdoTimer ) timer_stop(timer_finish);
  
  cdoFinish();
 
  return (0);
}

