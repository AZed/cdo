#ifndef _VINTERP_H
#define _VINTERP_H

#if defined(__cplusplus)
extern "C" {
#endif

void h2p(double *phlev, double *hlev, int nphlev);

void presh(double *pf, double *php, double *vct, double *ps, int nhlev, int ngp);

void genind(int *nx, double *plev, double *fullp, int ngp, int nplev, int nhlev);
void genindmiss(int *nx, double *plev, int ngp, int nplev, double *ps_prog, int *pnmiss);

void extra_P(double *slp, double *halfp, double *fullp, double *geop, double *temp, int ngp);

void interp_T(double *geop, double *gt, double *pt, double *fullp, double *halfp,
              int *nx, double *plev, int nplev, int ngp, int nhlev, double missval);
void interp_Z(double *geop, double *gz, double *pz, double *fullp, double *halfp,
	      int *nx, double *gt, double *plev, int nplev, int ngp, int nhlev, double missval);
void interp_X(double *gt, double *pt, double *hyb_press,
	      int *nx, double *plev, int nplev, int ngp, int nhlev, double missval);

#if defined(__cplusplus)
}
#endif

#endif  /* _VINTERP_H */
