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
     
     Timeof        eof             EOF in spatial or time space
     Timeof        eofspatial      EOF in spatial space
     Timeof        eoftime         EOF in time space 
*/
#define WEIGHTS 1

#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "statistic.h"

// NO MISSING VALUE SUPPORT ADDED SO FAR

void *Timeof(void * argument)
{
  static char func[] = "Timeof";
  enum {EOF_, EOF_TIME, EOF_SPATIAL};
  int operatorID;
  int operfunc;
  int streamID1, streamID2, streamID3;
  int gridsize;
  int vdate = 0, vtime = 0;
  int nrecs, nvars, nlevs=0 ;
  int i, ii, j, i1, i2, j1, j2;
  int nmiss;
  int tsID;
  int varID, recID, levelID;
  int vlistID1, vlistID2 = -1, vlistID3 = -1;
  int taxisID1, taxisID2, taxisID3;
  int gridID1, gridID2, gridID3;
  int ngrids;
  int reached_eof;
  int npack=0, nts=0;
  int *pack;
  int ***iwork;
  int n_eig, n=0;
  int grid_space=0, time_space=0;
  double *w;
  double sum_w;
  double sum;
  double missval=0;
  double *xvals, *yvals;
  FIELD ***fwork;
  FIELD ***o, ***o2; 
  FIELD in;     
  
  cdoInitialize(argument);
  cdoOperatorAdd("eof",       EOF_,       0, NULL);
  cdoOperatorAdd("eoftime",   EOF_TIME,   0, NULL);
  cdoOperatorAdd("eofspatial",EOF_SPATIAL,0, NULL);
  
  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);
  
  operatorInputArg("Number of eigen functions to write out");  
  n_eig = atoi(operatorArgv()[0]);
  
  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));
  vlistID1 = streamInqVlist(streamID1); 
  taxisID1 = vlistInqTaxis(vlistID1);
  
  gridsize = vlistGridsizeMax(vlistID1);
  gridID2 = gridID1 = vlistInqVarGrid(vlistID1, 0);
  nvars = vlistNvars(vlistID1);
  nrecs = vlistNrecs(vlistID1); 
  taxisID1 = vlistInqTaxis(vlistID1);
  w = (double*)malloc(gridsize*sizeof(double));
  gridWeights(gridID1, &w[0]);
  
  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));
  vlistID2 = vlistDuplicate(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1); 
  vlistDefTaxis(vlistID2, taxisID2);  
  
  
  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));
  vlistID3 = vlistDuplicate(vlistID1); 
  taxisID3 = taxisDuplicate(taxisID1);  
  vlistDefTaxis(vlistID3, taxisID3);  
  
  gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  xvals=(double*)malloc(1*sizeof(double));
  yvals=(double*)malloc(1*sizeof(double));
  xvals[0]=0;
  yvals[0]=0;
  gridDefXvals(gridID3, xvals);
  gridDefYvals(gridID3, yvals);
  
  ngrids = vlistNgrids(vlistID3);
  
  for ( i = 0; i < ngrids; i++ )
    vlistChangeGridIndex(vlistID3, i, gridID3);
  
  reached_eof=0;
  tsID = 0;
  
  /* COUNT NUMBER OF TIMESTEPS if EOF_ or EOF_TIME */
  
  if ( operfunc == EOF_ || operfunc == EOF_TIME)
    {
      while ( TRUE )
        {      
          if ( reached_eof ) continue;
          nrecs = streamInqTimestep(streamID1, tsID);
          if ( nrecs == 0 )
            {
              reached_eof = 1;
              break;
            }              
          tsID++;      
        }
      nts = tsID;
      reached_eof = 0;
      streamID1 = streamOpenRead(cdoStreamName(0));
      if ( nts < gridsize || operfunc == EOF_TIME) 
       {
         time_space = 1;
         grid_space = 0;
       }
      else 
        {
          time_space = 0;
          grid_space = 1;
        }
    }
  else if ( operfunc == EOF_SPATIAL )
    {
      time_space = 0;
      grid_space = 1;
    }  
  
  if ( time_space )    
    {
      if ( n_eig > nts )
        {
          cdoWarning("Solving in time-space:");
          cdoWarning("Number of eigen-functions to write out is bigger than number of time-steps.");
          cdoWarning("Setting n_eig to %i.", nts);
          cdoWarning("If You want to force a solution in grid-space use operator eofspatial");
          n_eig = nts;          
        }  
      n=nts;
    }
  else if ( grid_space )
    {
      if ( n_eig > gridsize )
        {
          cdoWarning("Solving in sptial space");
          cdoWarning("Number of eigen-functions to write out is bigger than number of time-steps");
          cdoWarning("Setting n_eig to %i", gridsize);
          cdoWarning("If You want to force a solution in time-space use operator eoftime");  
          n_eig = gridsize;
        }  
      n=gridsize;
    }
  // ALLOCATION OF FIELDS 
  in.ptr = (double *) malloc(gridsize*sizeof(double));
  
  fwork = (FIELD ***) malloc(nvars*sizeof(FIELD**));
  o     = (FIELD ***)malloc(nvars*sizeof(FIELD**));
  o2    = (FIELD ***)malloc(nvars*sizeof(FIELD**));
  iwork = (int ***) malloc(nvars*sizeof(int **));
  
  for ( varID = 0; varID < nvars; ++varID )
    {
      gridID1  = vlistInqVarGrid(vlistID1, varID);      
      gridsize = vlistGridsizeMax(vlistID1);
      nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval  = vlistInqVarMissval(vlistID1, varID);      
      
      fwork[varID] = (FIELD **)  malloc(nlevs*sizeof(FIELD*));
      o[varID]     = (FIELD **)  malloc(nlevs*sizeof(FIELD*));     
      o2[varID]    = (FIELD **)  malloc(nlevs*sizeof(FIELD*));
      iwork[varID] = (int  **)   malloc(nlevs*sizeof(int* ));
      
      for ( levelID = 0; levelID < nlevs; ++levelID )
        { 
          if ( grid_space )
            {
              fwork[varID][levelID] = (FIELD *) malloc (1*sizeof(FIELD));
              fwork[varID][levelID][0].grid    = gridID1;
              fwork[varID][levelID][0].nmiss   = 0;
              fwork[varID][levelID][0].missval = missval;
              fwork[varID][levelID][0].ptr     = (double *)malloc(gridsize*gridsize*sizeof(double));
              for ( i = 0; i < gridsize*gridsize; ++i )
                fwork[varID][levelID][0].ptr[i] = 0;
              iwork[varID][levelID] = (int *) malloc(gridsize*gridsize*sizeof(int));
              memset(iwork[varID][levelID], 0, gridsize*gridsize*sizeof(int)); 
            }
          else if ( time_space ) 
            {
              fwork[varID][levelID] = (FIELD *) malloc (nts*sizeof(FIELD));
              for ( tsID=0; tsID<nts; tsID++ )
                {
                  fwork[varID][levelID][tsID].grid    = gridID1;
                  fwork[varID][levelID][tsID].nmiss   = 0;
                  fwork[varID][levelID][tsID].missval = missval;
                  fwork[varID][levelID][tsID].ptr     = (double *)malloc(gridsize*sizeof(double));
                  for ( i = 0; i < gridsize; ++i )
                    fwork[varID][levelID][tsID].ptr[i] = 0;
                  
                }
              iwork[varID][levelID] = (int *) malloc (gridsize*sizeof(int));
            }                    
          
          o[varID][levelID] = (FIELD *) malloc(n_eig*sizeof(FIELD));
          o2[varID][levelID]= (FIELD *) malloc(gridsize*sizeof(FIELD));
          
          for ( i = 0; i < n; i++ )
            {
              if ( i < n_eig )
                {
                  o[varID][levelID][i].grid   = gridID2;
                  o[varID][levelID][i].nmiss  = 0;
                  o[varID][levelID][i].missval= missval;
                  o[varID][levelID][i].ptr    = (double *)malloc(gridsize*sizeof(double));
                  for ( ii = 0; ii < gridsize; ++ii )
                    o[varID][levelID][i].ptr[ii] = 0;
                }
              
              o2[varID][levelID][i].grid    = gridID3;
              o2[varID][levelID][i].nmiss   = 0;
              o2[varID][levelID][i].missval = missval;
              o2[varID][levelID][i].ptr     = (double *)malloc(1*sizeof(double));
              o2[varID][levelID][i].ptr[0]  = missval;              
            }
        }      
    }    
  
  tsID=0; 
  /* Read the data and create covariance matrices for each var & level */
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
          int i2;
          double tmp;
          streamInqRecord(streamID1, &varID, &levelID); 
          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));          
          missval = in.missval = vlistInqVarMissval(vlistID1, varID);           
          
          streamReadRecord(streamID1, in.ptr, &in.nmiss); 
          if ( grid_space )
            {
              for ( i1 = 0; i1 < gridsize; ++i1 )
                {                    
                  for ( i2 = i1; i2 <= gridsize; i2++ )
                    {                
                      if ( ( ! DBL_IS_EQUAL(in.ptr[i1], missval) ) && 
                          ( ! DBL_IS_EQUAL(in.ptr[i2], missval) ) )                            
                        {
                          tmp = in.ptr[i1]*in.ptr[i2];                    
                          fwork[varID][levelID][0].ptr[i1*gridsize+i2] += tmp;                      
                          iwork[varID][levelID][i1*gridsize+i2]++;                                        
                        }   
                      else
                        cdoWarning("Missing value support not checked for this operator!");
                    }           
                }	 
            }
          else if ( time_space )                    
            for ( i=0; i<gridsize; ++i )         
              {
                fwork[varID][levelID][tsID].ptr[i] = in.ptr[i];                      
                if ( ! DBL_IS_EQUAL(in.ptr[i], missval ) )
                  iwork[varID][levelID][i]++;
                else
                  cdoAbort("No missing value support for EOF in time space");
              }
                  
        }                                  
      tsID++;
    }
  if ( grid_space ) 
    for ( i1=0; i1<gridsize; ++i1 )
      for ( i2=0; i2<i1; ++i2 )                  
        {
          fwork[varID][levelID][0].ptr[i1*gridsize+i2] = 
          fwork[varID][levelID][0].ptr[i2*gridsize+i1];
          iwork[varID][levelID][i1*gridsize+i2] = 
          iwork[varID][levelID][i2*gridsize+i1];
        }
  
  pack = (int *)malloc(gridsize*sizeof(int));   
  
  for ( varID = 0; varID < nvars; varID++ )
    {         
      for ( levelID = 0; levelID < nlevs; levelID++ )
        {
          int i2;
          double **cov = NULL;
          double *eigv;
          npack = 0;
          sum_w = 0;          
          if ( grid_space ) 
            {      
              memset(pack, 0, gridsize*sizeof(int));                
              for ( i1 = 0; i1 < gridsize; i1++ )
                {        
                  for ( i2 = 0; i2 < gridsize; i2++ )              
                    if (!iwork[varID][levelID][i2*gridsize+i])
                      break;
                  
                  if ( i2 == gridsize )
                    pack[npack++] = i1;                                   
                }          
              for(i1=0;i1<npack;i1++)
                sum_w += w[pack[i1]];
              cov = (double **)malloc(npack*sizeof(double *));
              for (i1=0; i1<npack; i1++ )
                cov[i1] = (double*)malloc(npack*sizeof(double));
              
              for (i1=0;i1<npack;i1++) 
                {                                                       
                  for (i2 = i1; i2 < npack; i2++ )                               
                    if ( iwork[varID][levelID][pack[i1]*gridsize+pack[i2]] )
                      {                    
                        cov[i2][i1] = cov[i1][i2] = 
                        fwork[varID][levelID][0].ptr[pack[i1]*gridsize+pack[i2]]* // covariance
                        sqrt(w[pack[i1]]) * sqrt (w[pack[i2]]) / sum_w /          // weights
                        (iwork[varID][levelID][pack[i1]*gridsize+pack[i2]] - 1.); // number of data contributing                        
                        
                      }
                }
            }        
          else if ( time_space ) 
            {  
              sum_w = 0;
              for ( i=0; i < gridsize ;i++ )
                if ( iwork[varID][levelID][i] )                                    
                  {
                    pack[npack] = i;                                          
                    npack++;
                    sum_w += w[i];
                  }
              
              cov = (double **) malloc (nts*sizeof(double*)); 
              for ( j1=0; j1<nts; j1++)
                cov[j1] = (double*)malloc(nts*sizeof(double));
              for ( j1=0; j1<nts; j1++)
                {                  
                  for ( j2=j1; j2<nts; j2++ )
                    {
                      sum = 0;
                      for ( i=0; i<npack; i++ )
                        { 
                          sum   += w[pack[i]]*
                          fwork[varID][levelID][j1].ptr[pack[i]]*
                          fwork[varID][levelID][j2].ptr[pack[i]];    
                        }                      
                      cov[j2][j1] = cov[j1][j2] = sum / sum_w / (nts - 1);                      
                    }
                }             
            } 
          
          /* solve the eigen problem */
          eigv = (double *)malloc(n*sizeof(double));
          eigen_solution_of_symmetric_matrix(&cov[0], &eigv[0], n, n, func);
          /* cov contains the eigenvectors, eigv the eigenvalues */
          
          for (i=0; i<n; i++)            
            o2[varID][levelID][i].ptr[0] = eigv[i];
          
          for (i=0;i<n_eig;i++)
            {                                               
              if ( grid_space )
                {
                  for(j=0;j<npack;j++)
                    {
                      if ( ! WEIGHTS )
                        o[varID][levelID][i].ptr[pack[j]] = cov[i][j];                                
                      else if ( WEIGHTS )
                        o[varID][levelID][i].ptr[pack[j]] = cov[i][j] / sqrt(w[pack[j]] / sum_w);                  
                    }                                             
                }
              else if ( time_space )
                {                  
                  for ( i2=0; i2<npack; i2++ )
                    {
                      sum = 0;
                      for ( j=0; j<nts; j++ )
                        sum += fwork[varID][levelID][j].ptr[pack[i2]] * cov[i][j];                      
                      o[varID][levelID][i].ptr[pack[i2]] = sum;
                    }                  
                  // NORMALIZING
                  sum=0;                                    
                  for ( i2=0; i2<npack; i2++ )                    
                    sum += w[pack[i2]] * 
                    o[varID][levelID][i].ptr[pack[i2]] * 
                    o[varID][levelID][i].ptr[pack[i2]];                                
                  if ( sum > 0 ) 
                    {
                      sum = sqrt(sum / sum_w);
                      for( i2=0; i2<npack; i2++ )                        
                        o[varID][levelID][i].ptr[pack[i2]] /= sum;                                              
                    }
                  else
                    {
                      for( i2=0; i2<npack; i2++ )
                        o[varID][levelID][i].ptr[pack[i2]] = missval;
                    }
                }                  
            }
        }
    }
  
  streamDefVlist(streamID3, vlistID3);
  streamDefVlist(streamID2, vlistID2);  
  for ( tsID=0; tsID<n; tsID++ )
    {
      taxisDefVdate(taxisID3, 0); 
      taxisDefVtime(taxisID3, 0);
      streamDefTimestep(streamID3, tsID); 
      
      if ( tsID < n_eig ) 
        {
          taxisDefVdate(taxisID2, 0);
          taxisDefVtime(taxisID2, 0);  
          streamDefTimestep(streamID2, tsID); 
        }
      
      for ( varID = 0; varID < nvars; varID++ )
        {
          nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for ( levelID=0; levelID < nlevs;levelID++ )
            {    
              fprintf(stderr, "%4i %4i\n", tsID, n_eig);
              if ( tsID < n_eig )
                {
                  nmiss = 0;              
                  for ( i = 0; i < gridsize; i++ )
                    if ( DBL_IS_EQUAL(o[varID][levelID][tsID].ptr[i], missval) ) nmiss++;              
                  streamDefRecord(streamID2, varID, levelID);              
                  streamWriteRecord(streamID2, o[varID][levelID][tsID].ptr, nmiss);      
                }

              if ( DBL_IS_EQUAL(o2[varID][levelID][tsID].ptr[i], missval) ) nmiss = 1;
              else nmiss = 0;
              streamDefRecord(streamID3, varID, levelID);
              streamWriteRecord(streamID3, o2[varID][levelID][tsID].ptr,nmiss);

            }     
        }
    }
            
  for ( varID=0;varID<nvars;varID++)
    {
      for(levelID=0;levelID<nlevs;levelID++)
        {
          for(i=0;i<gridsize;i++)
            {
              if ( i < n_eig ) 
                free(o[varID][levelID][i].ptr);
              free(o2[varID][levelID][i].ptr);
            }
          free(o[varID][levelID]);
          free(o2[varID][levelID]);
          free(iwork[varID][levelID]);
        }
      free(o[varID]);
      free(o2[varID]);
      free(fwork[varID]);
      free(iwork[varID]);
      
    }
  
  free(o);
  free(o2);
  free(fwork);
  free(iwork);
  free(in.ptr);
  
  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);
  
  cdoFinish();   
  
  return (0);  
  
}
