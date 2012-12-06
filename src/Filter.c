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
 
      Filter    highpass
      Filter    lowpass
      Filter    bandpass
*/

#if defined ( _USE_FFTW3 ) 
#include <fftw3.h>
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#include "stdlib.h"


#define  NALLOC_INC  1000
#define  PI2         6.2832
#define  HALF        0.5

/* FAST FOURIER TRANSFORMATION (bare) */
static
void fft2(double *real, double *imag, int n, int isign)
{
  int nn, mmax, m, j, istep, i;
  double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi, tmp;   
  
  if ( n < 2 || n&(n-1) ) printf("n must be power of 2\n");
  nn = n << 1;
  j = 1;
  
  /* BIT Reversion of data */
  for ( i=1; i<nn; i+=2 )
    {
      if ( j > i ) 
        {        
          /* swap real part */
          tmp = real[j/2]; 
          real[j/2] = real[i/2];          
          real[i/2] = tmp;                    
          
          /* swap imaginary part */
          tmp = imag[j/2]; 
          imag[j/2] = imag[i/2];
          imag[i/2] = tmp;         
        }
      m = n;
      while ( m >= 2 && j > m )
        {
          j -= m;
          m >>= 1;
        }
      j += m;
    }
  
  /* Danielson-Lanzcos algorithm */
  mmax = 2;
  while ( nn > mmax )
    {
      istep = mmax << 1;
      theta = isign*(PI2/mmax);
      wtemp = sin(HALF*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      for( m = 1; m<mmax; m+=2 )
        {
          for ( i = m; i <= nn; i+=istep)
            {              
              j=i+mmax;         
              tempr = wr*real[j/2]-wi*imag[j/2];
              tempi = wr*imag[j/2]+wi*real[j/2];
              real[j/2] = real[i/2]-tempr;
              imag[j/2] = imag[i/2]-tempi;
              real[i/2] += tempr;
              imag[i/2] += tempi;       
            }
          wr = (wtemp=wr)*wpr-wi*wpi+wr;
          wi = wi*wpr+wtemp*wpi+wi;
        }
      mmax = istep;      
    }   
  if ( isign == -1 )
    for( i =0; i<n; i++)
      {
        real[i]/=n;
        imag[i]/=n;
      }
}


/* include from Tinfo.c */
void getTimeInc(double jdelta, int vdate0, int vdate1, int *incperiod, int *incunit);

static
void create_fmasc(int nts, double fdata, double fmin, double fmax, int *fmasc)
{
  double dimin, dimax;
  int i, imin, imax;
  
  dimin = nts*fmin / fdata;
  dimax = nts*fmax / fdata;
 
  imin = dimin<0 ? 0 : (int)floor(dimin);  
  imax = ceil(dimax)>nts/2 ? nts/2 : (int) ceil(dimax);  
  
  fmasc[imin] = 1;
  for ( i = imin+1; i <= imax; i++ )  
    fmasc[i] = fmasc[nts-i] = 1; 
  
}


#if defined ( _USE_FFTW3 ) 
static
void filter_fftw(int nts, const int *fmasc, 
	    fftw_complex *fft_in, fftw_complex *fft_out, fftw_plan *p_T2S, fftw_plan *p_S2T)
{  
  //  fprintf(stderr,"using fftw filter\n");

  int i;

  fftw_execute(*p_T2S);

  for ( i = 0; i < nts; i++ )
    if ( ! fmasc[i] )   {
      fft_out[i][0] = 0;
      fft_out[i][1] = 0;
    }
  
  fftw_execute(*p_S2T);
  
  return;
}

#else 

static
void filter_intrinsic(int nts, const int *fmasc, double *array1, double *array2)
{  
  int i;
  
  fft2(array1, array2, nts, 1);
  for ( i = 0; i < nts; i++ )
    if ( ! fmasc[i] )  array1[i] = array2[i] = 0;
  fft2(array1, array2, nts, -1);
  
  return;
}

#endif

void *Filter(void *argument)
{
  enum {BAND, HIGH, LOW};
  char *tunits[] = {"second", "minute", "hour", "day", "month", "year"};
  int iunits[] = {31536000, 525600, 8760, 365, 12, 1};
  int operatorID;
  int operfunc;
  int gridsize;
  int nrecs;
  int gridID, varID, levelID, recID;
  int tsID;
  int i;
  int nts,nts2;
  int nalloc = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  dtinfo_t *dtinfo = NULL;
  int tunit;
  int incperiod0, incunit0, incunit, dpy, calendar;
  int year0, month0, day0;
  double missval;
  double *array1, *array2;
  double fdata = 0;
  field_t ***vars = NULL;
  double fmin = 0, fmax = 0;
  int *fmasc;
#if defined ( _USE_FFTW3 ) 
  fftw_plan p_T2S, p_S2T;
  fftw_complex *out_fft;
  fftw_complex *in_fft;
#endif
  
  cdoInitialize(argument);

  cdoOperatorAdd("bandpass",  BAND,  0, NULL);
  cdoOperatorAdd("highpass",  HIGH,  0, NULL);
  cdoOperatorAdd("lowpass" ,  LOW,   0, NULL);

  operatorID = cdoOperatorID();
  operfunc   = cdoOperatorF1(operatorID);
  
  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  tunit = taxisInqTunit(taxisID1); 
  calendar = taxisInqCalendar(taxisID1);  
  dpy = calendar_dpy(calendar); /* should be 365 !!! */
  
  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  
  streamDefVlist(streamID2, vlistID2);
  
  nvars = vlistNvars(vlistID1);
  
  tsID = 0;    
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID >= nalloc )
        {
          nalloc += NALLOC_INC;
	  dtinfo = (dtinfo_t *)  realloc(dtinfo, nalloc*sizeof(dtinfo_t));
          vars   = (field_t ***) realloc(vars,   nalloc*sizeof(field_t **));
        }
                       
      taxisInqDTinfo(taxisID1, &dtinfo[tsID]);
   
      vars[tsID] = (field_t **) malloc(nvars*sizeof(field_t *));
      
      for ( varID = 0; varID < nvars; varID++ )
        {
          gridID   = vlistInqVarGrid(vlistID1, varID);
          missval  = vlistInqVarMissval(vlistID1, varID);
          nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          
          vars[tsID][varID] = (field_t *) malloc(nlevel*sizeof(field_t));
          
          for ( levelID = 0; levelID < nlevel; levelID++ )
            {
              vars[tsID][varID][levelID].grid    = gridID;
              vars[tsID][varID][levelID].missval = missval;
              vars[tsID][varID][levelID].ptr     = NULL;
            }
        }
      
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID1, &varID, &levelID);
          gridID   = vlistInqVarGrid(vlistID1, varID);
          gridsize = gridInqSize(gridID);
          vars[tsID][varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));
          streamReadRecord(streamID1, vars[tsID][varID][levelID].ptr, &nmiss);
          vars[tsID][varID][levelID].nmiss = nmiss;
          if ( nmiss ) cdoAbort("Missing value support for operators in module Filter not added yet");
        }

      /* get and check time increment */                   
      if ( tsID > 0)
        {    
	  juldate_t juldate0, juldate;
	  double jdelta;
	  int incperiod = 0;
	  int year, month, day;

          cdiDecodeDate(dtinfo[tsID].v.date,   &year,  &month,  &day);
	  cdiDecodeDate(dtinfo[tsID-1].v.date, &year0, &month0, &day0);               

          juldate0 = juldate_encode(calendar, dtinfo[tsID-1].v.date, dtinfo[tsID-1].v.time);        
          juldate  = juldate_encode(calendar, dtinfo[tsID].v.date, dtinfo[tsID].v.time);         
          jdelta   = juldate_to_seconds(juldate_sub(juldate, juldate0));
          
          if ( tsID == 1 ) 
            {           
              /*printf("%4i %4.4i-%2.2i-%2.2i\n", tsID, year, month, day);
              printf("    %4.4i-%2.2i-%2.2i\n",     year0,month0,day0);*/
              getTimeInc(jdelta, dtinfo[tsID-1].v.date, dtinfo[tsID].v.date, &incperiod0, &incunit0);
              incperiod = incperiod0; 
              if ( incperiod == 0 ) cdoAbort("Time step must be different from zero\n");
              incunit = incunit0;
              if ( cdoVerbose ) fprintf(stdout, "time step %i %s\n", incperiod, tunits[incunit]);
              fdata = 1.*iunits[incunit]/incperiod;
            }
          else 
            getTimeInc(jdelta, dtinfo[tsID-1].v.date, dtinfo[tsID].v.date, &incperiod, &incunit);  
        

	  if ( incunit0 < 4 && month == 2 && day == 29 && 
	       ( day0 != day || month0 != month || year0 != year ) )
	    {
	      cdoWarning("Filtering of multi-year times series only works properly with 365-day-calendar.");
	      cdoWarning("  Please delete the day %i-02-29 (cdo del29feb)", year);
	    }

          if ( ! ( incperiod == incperiod0 && incunit == incunit0 ) )
            cdoWarning("Time increment in step %i (%d%s) differs from step 0 (%d%s)",
		       tsID, incperiod, tunits[incunit], incperiod0, tunits[incunit0]);        
        }
      tsID++;
    }
  
  nts = tsID;
  /*  round up nts to next power of two for (better) performance 
   ** of fast fourier transformation */
#if defined ( _USE_FFTW3 ) 
  nts2 = nts;

  out_fft = (fftw_complex *) malloc ( nts * sizeof(fftw_complex) );
  in_fft  = (fftw_complex *) malloc ( nts * sizeof(fftw_complex) );

  p_T2S = fftw_plan_dft_1d(nts,in_fft,out_fft,  1, FFTW_ESTIMATE);
  p_S2T = fftw_plan_dft_1d(nts,out_fft,in_fft, -1, FFTW_ESTIMATE);
#else 
  nts2 = nts-1;
  nts2 |= nts2 >> 1;  /* handle  2 bit numbers */
  nts2 |= nts2 >> 2;  /* handle  4 bit numbers */
  nts2 |= nts2 >> 4;  /* handle  8 bit numbers */
  nts2 |= nts2 >> 8;  /* handle 16 bit numbers */
  nts2 |= nts2 >> 16; /* handle 32 bit numbers */
  nts2++;

  array1 = (double *) malloc(nts2*sizeof(double));
  array2 = (double *) malloc(nts2*sizeof(double));
#endif

  fmasc  = (int *) calloc(nts2, sizeof(int));
   
  for ( tsID = 0; tsID < nts; tsID++ ) array2[tsID] = 0;

  switch(operfunc)
  {
    case BAND: 
    {
      operatorInputArg("Lower bound of frequency band:");
      fmin = atof(operatorArgv()[0]);
      
      operatorInputArg("Upper bound of frequency band:");
      fmax = atof(operatorArgv()[1]);
      break;
    }
    case HIGH: 
    {              
      operatorInputArg("Lower bound of frequency pass:");
      fmin = atof(operatorArgv()[0]);
      fmax = fdata;
      break;
    }
    case LOW:  
    {
      fmin = 0;
      operatorInputArg("Upper bound of frequency pass:");
      fmax = atof(operatorArgv()[0]);
      break;
    }      
  }
  
  create_fmasc(nts, fdata, fmin, fmax, fmasc); 

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      missval  = vlistInqVarMissval(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

#if defined ( _USE_FFTW3 ) 
      fprintf(stderr," using fftw lib\n");
#endif
      
      for ( levelID = 0; levelID < nlevel; levelID++ )
        { 
#if defined ( _USE_FFTW3 ) 
          for ( i = 0; i < gridsize-1; i++ )
            {
              for ( tsID = 0; tsID < nts; tsID++ )                              
                {
		  in_fft[tsID][0] = vars[tsID][varID][levelID].ptr[i];
		  // in_fft[tsID][1] = vars[tsID][varID][levelID].ptr[i+1];
		  in_fft[tsID][1] = 0;
		}

	      filter_fftw(nts,fmasc,in_fft,out_fft,&p_T2S,&p_S2T);

              for ( tsID = 0; tsID < nts; tsID++ )
		{
		  vars[tsID][varID][levelID].ptr[i]   = in_fft[tsID][0] / nts;  
		  //		  vars[tsID][varID][levelID].ptr[i+1] = in_fft[tsID][1] / nts;  
		}
	    }
#else 
	  for ( i = 0; i < gridsize; i++ )  
	    {
	      // for some reason, the optimization using the complex transform independent of the 
	      // real one in order to transform two time series at the same time does not work
	      // properly here. 

	      memset( array2,0,nts2*sizeof(double) );
	      for ( tsID = 0; tsID < nts; tsID++ )
		array1[tsID] = vars[tsID][varID][levelID].ptr[i];                                         
	      /* zero padding up to next power of to */
              for ( ; tsID < nts2; tsID++ )                
		array1[tsID] = 0;       

	      filter_intrinsic(nts2,fmasc,array1,array2);
              for ( tsID = 0; tsID < nts; tsID++ )
		vars[tsID][varID][levelID].ptr[i]   = array1[tsID];  
	    }
#endif
	}
    }
  
  if ( array1 ) free(array1);
  if ( array2 ) free(array2);

  for ( tsID = 0; tsID < nts; tsID++ )
    {
      taxisDefDTinfo(taxisID2, dtinfo[tsID]);
      streamDefTimestep(streamID2, tsID);
    
      for ( varID = 0; varID < nvars; varID++ )
        {
          nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for ( levelID = 0; levelID < nlevel; levelID++ )
            {
              if ( vars[tsID][varID][levelID].ptr )
                {
                  nmiss = vars[tsID][varID][levelID].nmiss;
		  //fprintf(stderr, "%d %d %d %g\n", tsID, varID, levelID, vars[tsID][varID][levelID].ptr[0]);
		  streamDefRecord(streamID2, varID, levelID);
                  streamWriteRecord(streamID2, vars[tsID][varID][levelID].ptr, nmiss);
                  free(vars[tsID][varID][levelID].ptr);
                }
            }
          free(vars[tsID][varID]);
        }
      free(vars[tsID]);
    }

  if ( vars   ) free(vars);
  if ( dtinfo ) free(dtinfo);
  
  streamClose(streamID2);
  streamClose(streamID1);
  
  cdoFinish();
  
  return (0);
}
