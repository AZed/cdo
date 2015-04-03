#include <stdio.h>
#include <stdlib.h> /* exit */
#include <string.h>
#include <math.h>

#ifndef _VINTERP_H
#  include "vinterp.h"
#endif

#define  SCALEHEIGHT     (-7000.)
#define  SCALESLP        (101325.0)

#define  C_EARTH_GRAV    (9.80665)
#define  C_RKBOL         (1.380658e-23)     /* Boltzmann constant in J/K   */
#define  C_RNAVO         (6.0221367e+23)    /* Avogadro constant in 1/mol  */
#define  C_RMD           (28.9644)          /* molecular weight of dry air */
#define  C_R             (C_RKBOL * C_RNAVO)
#define  C_EARTH_RD      (1000. * C_R / C_RMD)

double Grav          = C_EARTH_GRAV;
double RD            = C_EARTH_RD;

int Mars = 0;


void h2p(double *phlev, double *hlev, int nphlev)
{
  int  k;
  double exp_arg;
  double height;

  for ( k = 0; k < nphlev; k++ )
    {
      height  = hlev[k];
      /*
	unitsel == 1 : hlev[k] is given in meters
	unitsel == 2 : hlev[k] is given in kilometers
	h2p needs meters (MKSC-standard)
      */

      exp_arg = height / SCALEHEIGHT;

      phlev[k] = SCALESLP * exp(exp_arg);
    }

}  /* h2p */


void p2h(double *hlev, double *plev, int nphlev)
{
  int  k;

  for ( k = 0; k < nphlev; k++ )
    {
      hlev[k] = log(plev[k]/SCALESLP)*SCALEHEIGHT;
    }

}  /* p2h */


void presh(double *fullp, double *halfp, double *vct, double *ps, int nhlev, int ngp)
{
  int i, lh;
  double zp, ze;
  double *halfpres = halfp;

  if ( ps == NULL )
    {
      fprintf(stderr, "ps undefined!\n");
      exit(EXIT_FAILURE);
    }

  for ( lh = 0; lh < nhlev; lh++ )
    {
      zp = vct[lh];
      ze = vct[lh+nhlev+1];

      for ( i = 0; i < ngp; i++ ) halfpres[i] = zp + ze * ps[i];

      halfpres += ngp;
    }

  memcpy(halfpres, ps, ngp*sizeof(double));

  halfpres = halfp;

  for ( i = 0; i < ngp*nhlev; i++ ) fullp[i] = 0.5 * (halfpres[i] + halfpres[i+ngp]);

} /* presh */


void genind(int *nx, double *plev, double *fullp, int ngp, int nplev, int nhlev)
{
  int  i, lp, lh;
  int *nxl;
  double pres, *fullpres;

  memset(nx, 0, ngp*nplev*sizeof(int));

#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(i, lh, pres, nxl)
#endif
  for ( lp = 0; lp < nplev; lp++ )
    {
      pres = plev[lp];
      nxl  = nx + lp*ngp;
      fullpres = fullp;
      for ( lh = 0; lh < nhlev; lh++ )
	for ( i = 0; i < ngp ; i++ )
	   {
	     if ( pres > fullpres[lh*ngp+i] ) nxl[i] = lh;
	   }
    }

}  /* genind */


void genindmiss(int *nx, double *plev, int ngp, int nplev, double *ps_prog, int *pnmiss)
{
  int i, lp;
  int *nxl;
  double pres;

#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(i, pres, nxl)
#endif
  for ( lp = 0; lp < nplev; lp++ )
    {
      pnmiss[lp] = 0;
      pres = plev[lp];
      nxl  = nx + lp*ngp;
      for ( i = 0; i < ngp; i++ )
	{
	  if ( pres > ps_prog[i] )
	    {
	      nxl[i] = -1;
	      pnmiss[lp]++;
	    }
	}
    }

}  /* genindmiss */


void extra_P(double *slp, double *halfp, double *fullp, double *geop, double *temp, int ngp)
{
  double alpha, tstar, tmsl, zprt, zprtal;
  double zrg;
  double zlapse = 0.0065;
  int j;

  zrg = 1.0 / Grav;

  for ( j = 0; j < ngp; ++j )
    {
      if ( geop[j] < 0.0001 && geop[j] > -0.0001 ) slp[j] = halfp[j];
      else
	{
	  alpha = RD * zlapse * zrg;
	  tstar = (1.0 + alpha * (halfp[j]/fullp[j] - 1.0)) * temp[j];

	  if ( tstar < 255.0 ) tstar = 0.5 * (255.0 + tstar);

	  tmsl = tstar + zlapse * zrg * geop[j];
	  if ( tmsl > 290.5 && tstar > 290.5 )
	    {
	      tstar = 0.5 * (290.5 + tstar);
	      tmsl  = tstar;
	    }

	  if ( tmsl-tstar < 0.000001 && tstar-tmsl < 0.000001 )
	    alpha = 0.0;
	  else if ( geop[j] > 0.0001 || geop[j] < -0.0001 )
	    alpha = RD * (tmsl-tstar) / geop[j];

	  zprt   = geop[j] / (RD * tstar);
	  zprtal = zprt * alpha;
	  slp[j] = halfp[j] * exp(zprt*(1.0-zprtal*(0.5-zprtal/3.0)));
	}
    }

}  /* extrap */


static 
double extra_T(double pres, double halfp, double fullp, double geop, double temp)
{
  double tstar, ztsz, z1, ztmsl, zalph, peval, zhts, zalp;
  double zrg;
  double zlapse = 0.0065;

  zrg   = 1.0 / Grav;
  tstar = (1.0 + zlapse * RD * zrg * (halfp/fullp - 1.0)) * temp;
  ztsz  = tstar;
  z1    = tstar + zlapse * zrg * geop;

  if ( tstar < 255.0 ) tstar = 0.5 * (255.0 + tstar);

  ztmsl = tstar + zlapse * zrg * geop;

  if ( ztmsl > 290.5 && tstar > 290.5 )
    {
      tstar = 0.5 * (290.5 + tstar);
      ztmsl = tstar;
    }

  if ( ztmsl > 290.5 && tstar <= 290.5 ) ztmsl=290.5;

  zalph = RD*zlapse*zrg;

  if ( ztmsl-tstar < 0.000001 && tstar-ztmsl < 0.000001 ) zalph=0.0;

  if ( (ztmsl-tstar > 0.000001 || tstar-ztmsl > 0.000001 ) &&
       (geop > 0.0001 || geop < -0.0001) )
    zalph = RD*(ztmsl-tstar)/geop;

  if ( pres <= halfp )
    peval = ((halfp-pres)*temp+ (pres-fullp)*tstar)/ (halfp-fullp);
  else
    {
      ztmsl = z1;
      tstar = ztsz;
      zhts  = geop * zrg;

      if ( zhts > 2000. && z1 > 298. )
	{
	  ztmsl = 298.;
	  if ( zhts < 2500. ) ztmsl = 0.002*((2500.-zhts)*z1+(zhts-2000.)*ztmsl);
	}

      if ( (ztmsl-tstar) < 0.000001 )
	zalph = 0.;
      else if (geop > 0.0001 || geop < -0.0001)
	zalph = RD*(ztmsl-tstar)/geop;
      else
	zalph = RD*zlapse*zrg;

      zalp  = zalph*log(pres/halfp);
      peval = tstar*(1.0+zalp*(1.0+zalp*(0.5+0.16666666667*zalp)));
    }

  return peval;

}  /* extra_T */


static 
double extra_Z(double pres, double halfp, double fullp, double geop, double temp)
{
  double alpha, tstar, tmsl, zalp, zalpal;
  double zrg;
  double zlapse = 0.0065;

  zrg   = 1.0 / Grav;
  alpha = RD * zlapse * zrg;
  tstar = (1.0 + alpha * (halfp/fullp - 1.0)) * temp;

  if ( tstar < 255.0 ) tstar = 0.5 * (255.0 + tstar);

  tmsl = tstar + zlapse * zrg * geop;

  if ( tmsl > 290.5 && tstar > 290.5 )
    {
      tstar = 0.5 * (290.5 + tstar);
      tmsl  = tstar;
    }

  if ( tmsl > 290.5 && tstar <= 290.5 ) tmsl = 290.5;

  if ( tmsl-tstar < 0.000001 && tstar-tmsl < 0.000001 )
    alpha = 0.0;
  else if ( geop > 0.0001 || geop < -0.0001 )
    alpha = RD * (tmsl-tstar) / geop;

  zalp   = log(pres/halfp);
  zalpal = zalp * alpha;

  return ((geop - RD*tstar*zalp*(1.0 + zalpal*(0.5 + zalpal/6.0)))*zrg);

}  /* extra_Z */


void interp_X(double *gt, double *pt, double *hyb_press, int *nx, double *plev, int nplev,
	     int ngp, int nhlev, double missval)
{
  int lp, i;
  int nl, nh;
  int *nxl;
  double *ptl;
  double pres;

#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(i, pres, nl, nh, nxl, ptl)
#endif
  for ( lp = 0; lp < nplev; lp++ )
    {
      pres = plev[lp];
      nxl  = nx + lp*ngp;
      ptl  = pt + lp*ngp;
      for ( i = 0; i < ngp; i++ )
	{
	  if ( nxl[i] == -1 )
	    ptl[i] = missval;
	  else
	    {
	      nl = nxl[i] * ngp + i;
	      nh = nl + ngp;
	      if ( nh >= ngp*nhlev )
		ptl[i] = gt[nl];
	      else
		ptl[i] =  gt[nl] + (pres-hyb_press[nl])
		       * (gt[nh] - gt[nl]) / (hyb_press[nh] - hyb_press[nl]);
	    }
	}
    }
}  /* interp_X */


void interp_T(double *geop, double *gt, double *pt, double *fullp, double *halfp,
              int *nx, double *plev, int nplev, int ngp, int nhlev, double missval)
{
  int lp, i;
  int nl, nh;
  int *nxl;
  double *ptl;
  double pres;

#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(i, pres, nl, nh, nxl, ptl)
#endif
  for ( lp = 0; lp < nplev; lp++ )
    {
      pres = plev[lp];
      nxl  = nx + lp*ngp;
      ptl  = pt + lp*ngp;
#if defined (CRAY)
#pragma _CRI inline extra_T
#endif
      for ( i = 0; i < ngp; i++ )
	{
	  nl = nxl[i];
	  if ( nl < 0 )
	    ptl[i] = missval;
	  else
	    {
	      if ( nl > nhlev-2 )
		{
		  if ( Mars )
		    ptl[i] = gt[(nhlev-1)*ngp+i];
		  else
#if defined (SX)
#pragma cdir inline
#endif
		    ptl[i] = extra_T(pres, halfp[nhlev*ngp+i],
				     fullp[(nhlev-1)*ngp+i], geop[i],
				     gt[(nhlev-1)*ngp+i]);
		}
	      else
		{
		  nh = nl + 1;
		  ptl[i] =  gt[nl*ngp+i] + (pres-fullp[nl*ngp+i])
                         * (gt[nh*ngp+i] - gt[nl*ngp+i])
                         / (fullp[nh*ngp+i] - fullp[nl*ngp+i]);
		}
	    }
	}
    }
}  /* interp_T */


void interp_Z(double *geop, double *gz, double *pz, double *fullp, double *halfp,
	      int *nx, double *gt, double *plev, int nplev, int ngp, int nhlev, double missval)
{
  int lp, i;
  int nl,nh;
  int *nxl;
  double *pzl;
  double pres;

#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(i, pres, nl, nh, nxl, pzl)
#endif
  for ( lp = 0; lp < nplev; lp++ )
    {
      pres = plev[lp];
      nxl  = nx + lp*ngp;
      pzl  = pz + lp*ngp;
#if defined (CRAY)
#pragma _CRI inline extra_Z
#endif
      for ( i = 0; i < ngp; i++ )
	{
	  nl = nxl[i];
	  if ( nl < 0 )
	    pzl[i] = missval;
	  else
	    {
	      if ( pres > halfp[(nl+1)*ngp+i] ) nl++;

	      if ( nl > nhlev-1 )
		{
		  if ( Mars )
		    pzl[i] = gt[(nhlev-1)*ngp+i];
		  else
#if defined (SX)
#pragma cdir inline
#endif
		    pzl[i] = extra_Z(pres, halfp[nhlev*ngp+i],
				     fullp[(nhlev-1)*ngp+i], geop[i],
				     gt[(nhlev-1)*ngp+i]);
		}
	      else
		{
		  nh = nl + 1;
		  pzl[i] =  gz[nl*ngp+i] + (pres-halfp[nl*ngp+i])
		         * (gz[nh*ngp+i] - gz[nl*ngp+i])
                         / (halfp[nh*ngp+i] - halfp[nl*ngp+i]);
		}
	    }
	}
    }
}  /* interp_Z */
