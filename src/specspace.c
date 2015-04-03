#include <stdio.h>
#include <math.h>
#include <string.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "specspace.h"
#include "error.h"
#include "grid.h"


#define  SQUARE_RADIUS   (-PlanetRadius * PlanetRadius)
#define  C_EARTH_RADIUS  (6371000.0)
double PlanetRadius = C_EARTH_RADIUS;


void geninx(long ntr, double *f, double *g)
{
  long m2,n2;
  long m, n ;

  for ( m = 0; m <= ntr; m++ )
    {
      m2 = m * m;
      for ( n = m; n <= ntr; n++ )
	{
	  n2 = n * n;
	  if ( n )
	    {
	      *g++ = -PlanetRadius / n * sqrt((double)(n2-m2)/(double)(4*n2-1));
	      *f++ = -PlanetRadius * m / (double)(n2+n);
	    }
	  else
	    {
	      *g++ = 0.0;
	      *f++ = 0.0;
	    }
	}
    }
}


void legini_old(int ntr, int nlat, double *poli, double *pold,
		double *pol2, double *pol3, double *coslat, double *rcoslat, int flag)
{
  int waves, dimsp;
  int jgl, jm, jn;
  int jsp;
  int pdim;
  double *gmu, *gwt, *pnm;
  double *hnm, gmusq, *ztemp1, *ztemp2;

  waves =  ntr + 1;
  dimsp = (ntr + 1) * (ntr + 2);
  pdim  = dimsp / 2 * nlat;

  gmu  = (double*) malloc(nlat * sizeof(double));
  gwt  = (double*) malloc(nlat * sizeof(double));

  gaussaw(gmu, gwt, nlat);

#if ! defined(_OPENMP)
  pnm    = (double*) malloc(dimsp * sizeof(double));
  hnm    = (double*) malloc(dimsp * sizeof(double));
  ztemp1 = (double*) malloc((waves<<1) * sizeof(double));
  ztemp2 = (double*) malloc((waves<<1) * sizeof(double));
#endif

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(jm, jn, jsp, gmusq, hnm, pnm, ztemp1, ztemp2)
#endif
  for ( jgl = 0; jgl < nlat; jgl++ )
    {
#if defined(_OPENMP)
      pnm    = (double*) malloc(dimsp * sizeof(double));
      hnm    = (double*) malloc(dimsp * sizeof(double));
      ztemp1 = (double*) malloc((waves<<1) * sizeof(double));
      ztemp2 = (double*) malloc((waves<<1) * sizeof(double));
#endif
      gmusq = 1.0 - gmu[jgl]*gmu[jgl];
      coslat[jgl] =  sqrt(gmusq);
      rcoslat[jgl] = 1.0 / coslat[jgl];

      phcs(pnm, hnm, waves, gmu[jgl], ztemp1, ztemp2);

      jsp = jgl;
      for ( jm = 0; jm < waves; jm++ )
	for ( jn = 0; jn < waves - jm; jn++ )
	  {
	              poli[jsp] = pnm[jm*waves+jn] * 2.0;
	              pold[jsp] = pnm[jm*waves+jn] * gwt[jgl];
	    if (flag) pol2[jsp] = hnm[jm*waves+jn] * gwt[jgl] /
                                  (PlanetRadius    * gmusq);
	    if (flag) pol3[jsp] = pnm[jm*waves+jn] * gwt[jgl] * jm /
                                  (PlanetRadius    * gmusq);
	    jsp += nlat;
	  }
#if defined(_OPENMP)
      free(ztemp2);
      free(ztemp1);
      free(pnm);
      free(hnm);
#endif
    }

#if ! defined(_OPENMP)
  free(ztemp2);
  free(ztemp1);
  free(pnm);
  free(hnm);
#endif
  free(gwt);
  free(gmu);
}


void legini(int ntr, int nlat, double *poli, double *pold, double *rcoslat)
{
  int waves, dimsp, dimpnm;
  int jgl, jm, jn, is;
  int isp, latn, lats;
  double *gmu, *gwt, *pnm, *work;

  waves  =  ntr + 1;
  dimsp  = (ntr + 1)*(ntr + 2);
  dimpnm = (ntr + 1)*(ntr + 4)/2;

  gmu  = (double*) malloc(nlat * sizeof(double));
  gwt  = (double*) malloc(nlat * sizeof(double));
  pnm  = (double*) malloc(dimpnm * sizeof(double));
  work = (double*) malloc(3*waves * sizeof(double));

  gaussaw(gmu, gwt, nlat);
  for ( jgl = 0; jgl < nlat; jgl++ ) gwt[jgl] *= 0.5;

  for ( jgl = 0; jgl < nlat; jgl++ )
    rcoslat[jgl] = 1.0 / sqrt(1.0 - gmu[jgl]*gmu[jgl]);

  for ( jgl = 0; jgl < nlat/2; jgl++ )
    {
      jspleg1(pnm, gmu[jgl], ntr, work);

      latn = jgl;
      isp = 0;
      for ( jm = 0; jm < waves; jm++ )
	{
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
	  for ( jn = 0; jn < waves - jm; jn++ )
	    {
	      is = (jn+1)%2 * 2 - 1;
	      lats = latn - jgl + nlat - jgl - 1;
	      poli[latn] = pnm[isp];
	      pold[latn] = pnm[isp] * gwt[jgl];
	      poli[lats] = pnm[isp] * is;
	      pold[lats] = pnm[isp] * gwt[jgl] * is;
	      latn += nlat;
	      isp++;
	    }
	  isp++;
	}
    }

  free(work);
  free(pnm);
  free(gwt);
  free(gmu);
}


void grid2spec(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int ntr, nlat, nlon, nfc;
  int nlev = 1;
  int waves;
  double *fpwork;
    
  ntr  = gridInqTrunc(gridIDout);
  nlon = gridInqXsize(gridIDin);
  nlat = gridInqYsize(gridIDin);

  waves = ntr + 1;
  nfc   = waves * 2;

  fpwork = (double*) malloc(nlat*nfc*nlev*sizeof(double));

  gp2fc(sptrans->trig, sptrans->ifax, arrayIn, fpwork, nlat, nlon, nlev, nfc);
  fc2sp(fpwork, arrayOut, sptrans->pold, nlev, nlat, nfc, ntr);

  free(fpwork);
}
	   
   
void spec2grid(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int ntr, nlat, nlon, nfc;
  int nlev = 1;
  int waves;
  double *fpwork;
    
  ntr  = gridInqTrunc(gridIDin);
  nlon = gridInqXsize(gridIDout);
  nlat = gridInqYsize(gridIDout);

  waves = ntr + 1;
  nfc   = waves * 2;

  fpwork = (double*) malloc(nlat*nfc*nlev*sizeof(double));

  sp2fc(arrayIn, fpwork, sptrans->poli, nlev, nlat, nfc, ntr);
  fc2gp(sptrans->trig, sptrans->ifax, fpwork, arrayOut, nlat, nlon, nlev, nfc);

  free(fpwork);
}


void four2spec(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int ntr, nlat, nfc;
  int nlev = 1;
  int waves;
    
  ntr  = gridInqTrunc(gridIDout);
  nlat = sptrans->nlat;

  waves = ntr + 1;
  nfc   = waves * 2;

  fc2sp(arrayIn, arrayOut, sptrans->pold, nlev, nlat, nfc, ntr);
}


void spec2four(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int ntr, nlat, nfc;
  int nlev = 1;
  int waves;
    
  ntr  = gridInqTrunc(gridIDin);
  nfc  = gridInqSize(gridIDout);
  nlat = nfc2nlat(nfc, ntr);

  waves = ntr + 1;
  nfc   = waves * 2;

  sp2fc(arrayIn, arrayOut, sptrans->poli, nlev, nlat, nfc, ntr);
}


void four2grid(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int ntr, nlat, nlon, nfc;
  int nlev = 1;
  int waves;

  ntr  = gridInqTrunc(gridIDin);
  nlon = gridInqXsize(gridIDout);
  nlat = gridInqYsize(gridIDout);

  waves = ntr + 1;
  nfc   = waves * 2;

  fc2gp(sptrans->trig, sptrans->ifax, arrayIn, arrayOut, nlat, nlon, nlev, nfc);
}


void grid2four(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int nlat, nlon, nfc, ntr;
  int nlev = 1;
  int waves;

  ntr  = gridInqTrunc(gridIDout);
  nlon = gridInqXsize(gridIDin);
  nlat = gridInqYsize(gridIDin);

  waves = ntr + 1;
  nfc   = waves * 2;

  gp2fc(sptrans->trig, sptrans->ifax, arrayIn, arrayOut, nlat, nlon, nlev, nfc);
}


void spec2spec(int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int ntrIn, ntrOut;

  ntrIn  = gridInqTrunc(gridIDin);
  ntrOut = gridInqTrunc(gridIDout);

  sp2sp(arrayIn, ntrIn, arrayOut, ntrOut);
}


void speccut(int gridIDin, double *arrayIn, double *arrayOut, int *waves)
{
  int ntr;

  ntr = gridInqTrunc(gridIDin);

  spcut(arrayIn, arrayOut, ntr, waves);
}


SPTRANS *sptrans_new(int nlon, int nlat, int ntr, int flag)
{
  SPTRANS *sptrans;
  int nsp;

  sptrans = (SPTRANS*) malloc(sizeof(SPTRANS));

  sptrans->nlon = nlon;
  sptrans->nlat = nlat;
  sptrans->ntr  = ntr;

  nsp = (ntr + 1)*(ntr + 2);
  sptrans->poldim = nsp / 2 * nlat;

  sptrans->trig = (double*) malloc(nlon * sizeof(double));
  fft_set(sptrans->trig, sptrans->ifax, nlon);

  sptrans->poli = (double*) malloc(sptrans->poldim * sizeof(double));
  sptrans->pold = (double*) malloc(sptrans->poldim * sizeof(double));
  if ( flag )
    {
      sptrans->pol2 = (double*) malloc(sptrans->poldim * sizeof(double));
      sptrans->pol3 = (double*) malloc(sptrans->poldim * sizeof(double));
    }
  else
    {
      sptrans->pol2 = NULL;
      sptrans->pol3 = NULL;
    }

  sptrans->coslat  = (double*) malloc(nlat * sizeof(double));
  sptrans->rcoslat = (double*) malloc(nlat * sizeof(double));

  if ( flag )
    legini_old(ntr, nlat, sptrans->poli, sptrans->pold,
	       sptrans->pol2, sptrans->pol3, sptrans->coslat, sptrans->rcoslat, flag);
  else
    legini(ntr, nlat, sptrans->poli, sptrans->pold, sptrans->rcoslat);

  return (sptrans);
}


void sptrans_delete(SPTRANS *sptrans)
{
  if ( sptrans )
    {
      if ( sptrans->trig ) { free(sptrans->trig);  sptrans->trig = NULL; }
      if ( sptrans->poli ) { free(sptrans->poli);  sptrans->poli = NULL; }
      if ( sptrans->pold ) { free(sptrans->pold);  sptrans->pold = NULL; }
      if ( sptrans->pol2 ) { free(sptrans->pol2);  sptrans->pol2 = NULL; }
      if ( sptrans->pol3 ) { free(sptrans->pol3);  sptrans->pol3 = NULL; }
      if ( sptrans->coslat  ) { free(sptrans->coslat);   sptrans->coslat = NULL; }
      if ( sptrans->rcoslat ) { free(sptrans->rcoslat);  sptrans->rcoslat = NULL; }

      free(sptrans); sptrans = NULL;
    }
}


DVTRANS *dvtrans_new(int ntr)
{
  DVTRANS *dvtrans;
  int dimsp;

  dvtrans = (DVTRANS*) malloc(sizeof(DVTRANS));

  dvtrans->ntr = ntr;

  dimsp = (ntr + 1)*(ntr + 2);
  dvtrans->fdim = dimsp / 2;

  dvtrans->f1 = (double*) malloc(dvtrans->fdim * sizeof(double));
  dvtrans->f2 = (double*) malloc(dvtrans->fdim * sizeof(double));

  geninx(ntr, dvtrans->f1, dvtrans->f2);

  return (dvtrans);
}


void dvtrans_delete(DVTRANS *dvtrans)
{
  if ( dvtrans )
    {
      if ( dvtrans->f1 ) { free(dvtrans->f1);  dvtrans->f1 = NULL; }
      if ( dvtrans->f2 ) { free(dvtrans->f2);  dvtrans->f2 = NULL; }

      free(dvtrans); dvtrans = NULL;
    }
}


void uv2dv(double *fu, double *fv, double *sd, double *sv,
           double *pol2, double *pol3, int klev, int nlat, int nt)
{
  int lev, jmm, jfc, lat, nfc, nsp2;
  double dir, dii, vor, voi;
  double *ufr, *ufi, *vfr, *vfi;
  double *ful, *fvl, *sdl, *svl;
  double *po2, *po3;

  nsp2 = (nt+1)*(nt+2);
  nfc  = (nt+1)*2;

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(jmm, jfc, lat, po2, po3, ful, fvl, sdl, svl, ufr, ufi, vfr, vfi, dir, dii, vor, voi)
#endif
  for ( lev = 0; lev < klev; lev++ )
    {
      po2 = pol2;
      po3 = pol3;
      ful = fu + lev*nfc*nlat;
      fvl = fv + lev*nfc*nlat;
      sdl = sd + lev*nsp2;
      svl = sv + lev*nsp2;
      for ( jmm = 0; jmm <= nt; jmm++ )
	{
	  for  ( jfc = jmm; jfc <= nt; jfc++ )
	    {
	      ufr = ful;
	      ufi = ful + nlat;
	      vfr = fvl;
	      vfi = fvl + nlat;
	      dir = 0.0;
	      dii = 0.0;
	      vor = 0.0;
	      voi = 0.0;
	      for ( lat = 0; lat < nlat; lat++ )
		{
		  dir += vfr[lat] * po2[lat] - ufi[lat] * po3[lat];
		  dii += vfi[lat] * po2[lat] + ufr[lat] * po3[lat];
		  vor -= ufr[lat] * po2[lat] + vfi[lat] * po3[lat];
		  voi -= ufi[lat] * po2[lat] - vfr[lat] * po3[lat];
		}
	      *sdl++ = dir;
	      *sdl++ = dii;
	      *svl++ = vor;
	      *svl++ = voi;
	      po2 += nlat;
	      po3 += nlat;
	    }
	  ful += 2 * nlat;
	  fvl += 2 * nlat;
	}
    }
}


void dv2uv(double *d, double *o, double *u, double *v, double *f, double *g,
	   int nt, int nsp, int nlev)
{
  /* d(nsp,nlev), o(nsp,nlev)     ! divergence, vorticity        */
  /* u(nsp,nlev), v(nsp,nlev)     ! zonal wind, meridional wind  */
  /* f(nsp/2)   , g(nsp/2)        ! factor tables                */

  int l, m, n;
  int i;

  for ( l = 0; l < nlev; l++ )
    {
      i = 0;

      for ( m = 0; m < nt-1; m++ )
	{
	  /*********/
	  /* n = m */
	  /*********/

	  if ( m == 0 )
	    {
	      *u++ = -g[i+1] * o[2*(i+1)  ];
	      *u++ = -g[i+1] * o[2*(i+1)+1];
	      *v++ =  g[i+1] * d[2*(i+1)  ];
	      *v++ =  g[i+1] * d[2*(i+1)+1];
	    }
	  else
	    {
	      *u++ = -f[i] * d[2*i+1] - g[i+1] * o[2*(i+1)  ];
	      *u++ =  f[i] * d[2*i  ] - g[i+1] * o[2*(i+1)+1];
	      *v++ = -f[i] * o[2*i+1] + g[i+1] * d[2*(i+1)  ];
	      *v++ =  f[i] * o[2*i  ] + g[i+1] * d[2*(i+1)+1];
	    }
	  ++i;

	  /****************/
	  /* m < n < nt-1 */
	  /****************/

	  for ( n = m+1; n < nt-1; n++ )
	    {
	      *u++ =  g[i] * o[2*(i-1)  ] - f[i] * d[2*i+1] - g[i+1] * o[2*(i+1)  ];
	      *u++ =  g[i] * o[2*(i-1)+1] + f[i] * d[2*i  ] - g[i+1] * o[2*(i+1)+1];
	      *v++ = -g[i] * d[2*(i-1)  ] - f[i] * o[2*i+1] + g[i+1] * d[2*(i+1)  ];
	      *v++ = -g[i] * d[2*(i-1)+1] + f[i] * o[2*i  ] + g[i+1] * d[2*(i+1)+1];
	      ++i;
	    }

	  /************/
	  /* n = nt-1 */
	  /************/

	  *u++ =  g[i] * o[2*(i-1)  ] - f[i] * d[2*i+1];
	  *u++ =  g[i] * o[2*(i-1)+1] + f[i] * d[2*i  ];
	  *v++ = -g[i] * d[2*(i-1)  ] - f[i] * o[2*i+1];
	  *v++ = -g[i] * d[2*(i-1)+1] + f[i] * o[2*i  ];
	  ++i;

	  /**********/
	  /* n = nt */
	  /**********/

	  *u++ =  g[i] * o[2*(i-1)  ];
	  *u++ =  g[i] * o[2*(i-1)+1];
	  *v++ = -g[i] * d[2*(i-1)  ];
	  *v++ = -g[i] * d[2*(i-1)+1];
	  ++i;
	}

      /***************************/
      /* m = nt-1  and  n = nt-1 */
      /***************************/

      *u++ = -f[i] * d[2*i+1];
      *u++ =  f[i] * d[2*i  ];
      *v++ = -f[i] * o[2*i+1];
      *v++ =  f[i] * o[2*i  ];
      ++i;

      /*************************/
      /* m = nt-1  and  n = nt */
      /*************************/

      *u++ =  g[i] * o[2*(i-1)  ];
      *u++ =  g[i] * o[2*(i-1)+1];
      *v++ = -g[i] * d[2*(i-1)  ];
      *v++ = -g[i] * d[2*(i-1)+1];
      ++i;

      /***********************/
      /* m = nt  and  n = nt */
      /***********************/

      *u++ = 0.0;
      *u++ = 0.0;
      *v++ = 0.0;
      *v++ = 0.0;

      d += nsp;
      o += nsp;
    }
}


void dv2ps(const double * restrict div, double * restrict pot, long nlev, long ntr)
{
  long l, m, n;
  double fact;

  for ( l = 0; l <  nlev; l++ )
    {
      /* m == 0 */
      *pot++ = 0.0;
      *pot++ = 0.0;
      div += 2;

      for ( n = 1; n <= ntr; n++ )
	{
	  fact = SQUARE_RADIUS / (n*n + n);
	  *pot++ = *div++ * fact;
	  *pot++ = *div++ * fact;
	}

      /* m >= 0 */
      for ( m = 1; m <= ntr; m++ )
	for ( n = m; n <= ntr; n++ )
	  {
	    fact = SQUARE_RADIUS / (n*n + n);
	    *pot++ = *div++ * fact;
	    *pot++ = *div++ * fact;
	  }
    }
}


void scaluv(double *fu, double *rclat, int nlat, int lot)
{
  int l,lat;

  for (l = 0; l < lot; l++)
    for (lat = 0; lat < nlat; lat++)
      {
        *fu *= rclat[lat];
        fu++;
      }
}


void trans_uv2dv(SPTRANS *sptrans, int nlev,
		 int gridID1, double *gu, double *gv,
		 int gridID2, double *sd, double *svo)
{
  int ntr, nlat, nlon, nfc;
  int waves;
  double *fpwork1, *fpwork2;

  if ( gridInqType(gridID1) != GRID_GAUSSIAN )
    Error("unexpected grid1 type: %s instead of Gaussian", gridNamePtr(gridInqType(gridID1)));

  if ( gridInqType(gridID2) != GRID_SPECTRAL )
    Error("unexpected grid2 type: %s instead of spectral", gridNamePtr(gridInqType(gridID2)));
    
  ntr  = gridInqTrunc(gridID2);
  nlon = gridInqXsize(gridID1);
  nlat = gridInqYsize(gridID1);

  waves = ntr + 1;
  nfc   = waves * 2;

  fpwork1 = (double*) malloc(nlat*nfc*nlev*sizeof(double));
  fpwork2 = (double*) malloc(nlat*nfc*nlev*sizeof(double));

  gp2fc(sptrans->trig, sptrans->ifax, gu, fpwork1, nlat, nlon, nlev, nfc);
  gp2fc(sptrans->trig, sptrans->ifax, gv, fpwork2, nlat, nlon, nlev, nfc);

  scaluv(fpwork1, sptrans->coslat, nlat, nfc*nlev);
  scaluv(fpwork2, sptrans->coslat, nlat, nfc*nlev);

  uv2dv(fpwork1, fpwork2, sd, svo, sptrans->pol2, sptrans->pol3, nlev, nlat, ntr);

  free(fpwork1);
  free(fpwork2);
}


void trans_dv2uv(SPTRANS *sptrans, DVTRANS *dvtrans, int nlev,
		 int gridID1, double *sd, double *svo,
		 int gridID2, double *gu, double *gv)
{
  int ntr, nlat, nlon, nfc;
  int waves;
  int dimsp;
  double *fpwork;
  double *su, *sv;

  if ( gridInqType(gridID1) != GRID_SPECTRAL )
    Warning("unexpected grid1 type: %s", gridNamePtr(gridInqType(gridID1)));

  if ( gridInqType(gridID2) != GRID_GAUSSIAN )
    Warning("unexpected grid2 type: %s", gridNamePtr(gridInqType(gridID2)));

  ntr  = gridInqTrunc(gridID1);
  nlon = gridInqXsize(gridID2);
  nlat = gridInqYsize(gridID2);

  waves = ntr + 1;
  nfc   = waves * 2;

  dimsp = (ntr + 1)*(ntr + 2);

  su = gu;
  sv = gv;

  dv2uv(sd, svo, su, sv, dvtrans->f1, dvtrans->f2, ntr, dimsp, nlev);

  fpwork = (double*) malloc(nlat*nfc*nlev*sizeof(double));

  sp2fc(su, fpwork, sptrans->poli, nlev, nlat, nfc, ntr);
  scaluv(fpwork, sptrans->rcoslat, nlat, nfc*nlev);
  fc2gp(sptrans->trig, sptrans->ifax, fpwork, gu, nlat, nlon, nlev, nfc);

  sp2fc(sv, fpwork, sptrans->poli, nlev, nlat, nfc, ntr);
  scaluv(fpwork, sptrans->rcoslat, nlat, nfc*nlev);
  fc2gp(sptrans->trig, sptrans->ifax, fpwork, gv, nlat, nlon, nlev, nfc);

  free(fpwork);
}

