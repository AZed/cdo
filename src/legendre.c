#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h> /* malloc zhlp3 */


void jspleg1(double *pleg, double plat, int ktrunc, double *work)
{
  /*
     jspleg1 - Routine to calculate legendre functions

     Purpose
     --------

     This routine calculates the legendre functions for one latitude.
     (but not their derivatives)


     Interface
     ----------

     jspleg1( pleg, plat, ktrunc)


     Input parameters
     ----------------

     plat      - Latitude in radians
     ktrunc    - Spectral truncation


     Output parameters
     -----------------

     pleg      - Array of legendre functions for one latitude.
                 The array must be at least (KTRUNC+1)*(KTRUNC+4)/2 
                 words long.

     Method
     ------

     Recurrence relation with explicit relations for P(m,m) and 
     P(m,m+1)


     AUTHOR
     ------

     J.D.Chambers         ECMWF        9 November 1993


     Modifications
     -------------

     None

  */
  int itout1, i1m, ilm, jm, jcn, im2;
  double zsin, zcos, zf1m, zre1, zf2m, zn, ze1, ze2;
  double *zhlp1, *zhlp2, *zhlp3;


  /* Initialization */

  itout1 = ktrunc+1;
  /*  zsin   = sin(plat); */
  zsin   = plat;
  zcos   = sqrt(1.-zsin*zsin);

  zhlp1 = work;
  zhlp2 = work + itout1;
  zhlp3 = work + itout1 + itout1;

  /*  Step 1.        M = 0, N = 0 and N = 1 */

  ilm     = 1;
  pleg[0] = 1.0;
  zf1m    = sqrt(3.0);
  pleg[1] = zf1m*zsin;

  /*  Step 2.       Sum for M = 0 to T (T = truncation) */

  for ( jm = 1; jm < itout1; jm++ )
    {
      zhlp1[jm] = sqrt(2.*jm+3.);
      zhlp2[jm] = 1./sqrt(2.*jm);
    }

  zhlp1[0] = sqrt(3.);

  for ( jm = 0; jm < itout1; jm++ )
    {
      i1m  = jm - 1;
      zre1 = zhlp1[jm];
      ze1  = 1./zre1;

      /*   Step 3.       M > 0 only */

      if ( i1m != -1 )
	{
          zf2m = zf1m*zcos*zhlp2[jm];
          zf1m = zf2m*zre1;

	  /*  Step 4.       N = M and N = M+1 */

          ilm       = ilm+1;
          pleg[ilm] = zf2m;
          ilm       = ilm+1;
          pleg[ilm] = zf1m*zsin;

	  /* When output truncation is reached, return to calling program */

          if ( jm == (itout1-1) ) break;
	}

      /*  Step 5.       Sum for N = M+2 to T+1 */

      im2 = i1m+2;

      for ( jcn = im2; jcn < itout1; jcn++ )
	{
          zn         = jcn + 1;
	  zhlp3[jcn] = sqrt((4.*zn*zn-1.)/(zn*zn-jm*jm));
	}

      for ( jcn = im2; jcn < itout1; jcn++ )
	{
          ze2        = zhlp3[jcn];
          ilm        = ilm+1;
          pleg[ilm]  = ze2*(zsin*pleg[ilm-1]-ze1*pleg[ilm-2]);
          ze1        = 1./ze2;
	}
    }
}


/* ============================================= */
/* phcs - Compute values of Legendre polynomials */
/*        and their meridional derivatives       */
/* ============================================= */

void phcs(double *pnm, double *hnm, int waves, double pmu,
	  double *ztemp1, double *ztemp2)
{
  int twowaves;

  int jk, jn, jm;

  double jnmjk;
  double zcos2;
  double lat;
  double zan;
  double zsinpar;
  double zcospar;
  double zsqp;
  double zcosfak;
  double zsinfak;
  double zq;
  double zwm2;
  double zw;
  double zwq;
  double zq2m1;
  double zwm2q2;
  double z2q2;
  double zcnm;
  double zdnm;
  double zenm;

  twowaves  = waves << 1;

  zcos2     = sqrt(1.0 - pmu*pmu);
  lat       = acos(pmu);
  zan       = 1.0;

  ztemp1[0] = 0.5;

  for ( jn = 1; jn < twowaves; jn++ )
    {
      zsqp    = 1.0 / sqrt((double)(jn + jn*jn));
      zan    *= sqrt(1.0 - 1.0/(4*jn*jn));

      zcospar = cos(lat * jn);
      zsinpar = sin(lat * jn) * jn * zsqp;
      zcosfak = 1.0;

      for ( jk = 2; jk < jn; jk += 2 )
	{
	  jnmjk = jn - jk;
	  zcosfak *= (jk-1.0) * (jn+jnmjk+2.0) / (jk * (jn+jnmjk+1.0));
	  zsinfak  = zcosfak * (jnmjk) * zsqp;
	  zcospar += zcosfak * cos(lat * jnmjk);
	  zsinpar += zsinfak * sin(lat * jnmjk);
	}

      /*  code for jk == jn */

      if ((jn & 1) == 0)
	{
	  zcosfak *= (double)((jn-1) * (jn+2)) / (double)(jn * (jn+1));
	  zcospar += zcosfak * 0.5;
	}
      ztemp1[jn  ] = zan * zcospar;
      ztemp2[jn-1] = zan * zsinpar;
    }

  memcpy(pnm, ztemp1, waves*sizeof(double));
  pnm += waves;
  memcpy(pnm, ztemp2, waves*sizeof(double));
  pnm += waves;

  hnm[0] = 0.0;
  for (jn = 1; jn < waves; jn++)
    hnm[jn] = jn * (pmu * ztemp1[jn] - sqrt((jn+jn+1.0) / (jn+jn-1.0)) * ztemp1[jn-1]);

  hnm += waves;

  hnm[0] = pmu * ztemp2[0];

  for (jn = 1; jn < waves; jn++)
    hnm[jn] = (jn+1) * pmu * ztemp2[jn]
            - sqrt(((jn+jn+3.0)*((jn+1)*(jn+1)-1.0)) / (jn+jn+1.0)) * ztemp2[jn-1];
	    
  hnm += waves;

  for (jm = 2; jm < waves; jm++)
    {
      pnm[0] = sqrt(1.0 + 1.0 / (jm+jm)) * zcos2 * ztemp2[0];
      hnm[0] = jm * pmu * pnm[0];
#if defined (CRAY)
#pragma _CRI novector
#endif
#if defined (__uxp__)
#pragma loop scalar
#endif
      for (jn = 1; jn < twowaves-jm; jn++)
	{
          zq      = jm + jm + jn - 1;
          zwm2    = zq + jn;
          zw      = zwm2 + 2;
          zwq     = zw*zq;
          zq2m1   = zq*zq - 1.;
          zwm2q2  = zwm2*zq2m1;
          z2q2    = zq2m1*2;
          zcnm    = sqrt((zwq*(zq-2.))/(zwm2q2-z2q2));
          zdnm    = sqrt((zwq*(jn+1.))/zwm2q2);
          zenm    = sqrt(zw * jn /((zq+1.0) * zwm2));
          pnm[jn] = zcnm * ztemp1[jn] - pmu
                  * (zdnm * ztemp1[jn+1] - zenm * pnm[jn-1]);
          hnm[jn] = (jm + jn) * pmu * pnm[jn]
                  - sqrt(zw * jn * (zq+1) / zwm2) * pnm[jn-1];
	}
      memcpy(ztemp1, ztemp2, twowaves*sizeof(double));
      memcpy(ztemp2, pnm   , twowaves*sizeof(double));
      pnm += waves;
      hnm += waves;
    }
}


/* to slow for nec, 2.0 instead of 2.3 GFlops ( vector length too small ) */
void sp2fctest(double *sa, double *fa, double *poli, int nlev, int nlat, int nfc, int nt)
{
  int lev, jm, jn, latn, lats, nsp2, is;
  double sar, sai;
  double *far, *fai, *pol;
  double *sal, *fal;
  double pval;

  nsp2 = (nt+1)*(nt+2);

  for ( lev = 0; lev < nlev; lev++ )
    {
      pol = poli;
      fal = fa + lev*nfc*nlat;
      sal = sa + lev*nsp2;
      memset(fal, 0, nfc*nlat*sizeof(double));

      for ( jm = 0; jm <= nt; jm++ )
	{
	  for ( jn = 0; jn <= nt - jm; jn++ )
	    {
	      is = (jn+1)%2 * 2 - 1;
	      sar = *sal++;
	      sai = *sal++;
	      far = fal;
	      fai = fal + nlat;
#if defined (SX)
#pragma vdir nodep
#endif
	      for ( latn = 0; latn < nlat/2; latn++ )
		{
		  lats = nlat - latn - 1;
		  pval = pol[latn];
		  far[latn] += pval * sar;
		  fai[latn] += pval * sai;
		  far[lats] += pval * sar * is;
		  fai[lats] += pval * sai * is;
		}
	      pol += nlat;
	    }
	  fal += 2 * nlat;
	}
    }
}


void sp2fc(const double *sa, double *fa, const double *poli, long nlev, long nlat, long nfc, long nt)
{
  long lev, jmm, jfc, lat, nsp2;
  double sar, sai;
  double *fal;
  double * restrict far, * restrict fai;
  const double * restrict pol;
  const double * restrict sal;

  nsp2 = (nt+1)*(nt+2);

#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(jmm, jfc, lat, pol, sar, sai, sal, far, fai, fal)
#endif
  for ( lev = 0; lev < nlev; lev++ )
    {
      pol = poli;
      fal = fa + lev*nfc*nlat;
      sal = sa + lev*nsp2;
      memset(fal, 0, nfc*nlat*sizeof(double));

      for ( jmm = 0; jmm <= nt; jmm++ )
	{
	  for ( jfc = jmm; jfc <= nt; jfc++ )
	    {
	      sar = *sal++;
	      sai = *sal++;
	      far = fal;
	      fai = fal + nlat;
	      for ( lat = 0; lat < nlat; lat++ )
		{
		  far[lat] += pol[lat] * sar;
		  fai[lat] += pol[lat] * sai;
		}
	      pol += nlat;
	    }
	  fal += 2 * nlat;
	}
    }
}


void fc2sp(double *fa, double *sa, double *poli, int nlev, int nlat, int nfc, int nt)
{
  int lev, jmm, jfc, lat, nsp2;
  double sar, sai, *far, *fai, *pol;
  double *sal, *fal;

  nsp2 = (nt+1)*(nt+2);

#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(jmm, jfc, lat, pol, sar, sai, sal, far, fai, fal)
#endif
  for ( lev = 0; lev < nlev; lev++ )
    {
      pol = poli;
      fal = fa + lev*nfc*nlat;
      sal = sa + lev*nsp2;
      for ( jmm = 0; jmm <= nt; jmm++ )
	{
	  for ( jfc = jmm; jfc <= nt; jfc++ )
	    {
	      far = fal;
	      fai = fal + nlat;
	      sar = 0.0;
	      sai = 0.0;
	      for ( lat = 0; lat < nlat; lat++ )
		{
		  sar += pol[lat] * far[lat];
		  sai += pol[lat] * fai[lat];
		}
	      *sal++ = sar;
	      *sal++ = sai;
	      pol += nlat;
	    }
	  fal += 2 * nlat;
	}
    }
}

/* ======================================== */
/* Convert Spectral Array to new truncation */
/* ======================================== */

void sp2sp(double *arrayIn, int truncIn, double *arrayOut, int truncOut)
{
  int n, m;

  if ( truncOut <= truncIn )
    {
      for ( n = 0; n <= truncOut; n++ )
	{
	  for ( m = n; m <= truncOut; m++ )
	    {
	      *arrayOut++ = *arrayIn++ ;
	      *arrayOut++ = *arrayIn++ ;
	    }
	  arrayIn += 2 * (truncIn-truncOut);
	}
    }
  else
    {
      for ( n = 0; n <= truncIn; n++ )
	{
	  for ( m = n; m <= truncIn; m++ )
	    {
	      *arrayOut++ = *arrayIn++ ;
	      *arrayOut++ = *arrayIn++ ;
	    }
	  for ( m = truncIn+1; m <= truncOut; ++m )
	    {
	      *arrayOut++ = 0.0;
	      *arrayOut++ = 0.0;
	    }
	}
      for ( n = truncIn+1; n <= truncOut; ++n )
	for ( m = n; m <= truncOut; ++m )
	  {
	    *arrayOut++ = 0.0;
	    *arrayOut++ = 0.0;
	  }
    }
}

/* ======================================== */
/* Cut spectral wave numbers                */
/* ======================================== */

void spcut(double *arrayIn, double *arrayOut, int trunc, int *waves)
{
  int n, m;

  for ( n = 0; n <= trunc; n++ )
    {
      for ( m = n; m <= trunc; m++ )
	{
	  if ( waves[m] )
	    {
	      *arrayOut++ = *arrayIn++;
	      *arrayOut++ = *arrayIn++;
	    }
	  else
	    {
	      *arrayOut++ = 0.0;
	      *arrayOut++ = 0.0;
	      arrayIn++;
	      arrayIn++;
	    }
	}
    }
}
