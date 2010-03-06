#ifndef _GRID_H
#define _GRID_H

void gridToDegree(const char *units, const char *string, int gridsize, double *array);
int gridToZonal(int gridID);
int gridToMeridional(int gridID);
int gridToCell(int gridID);
int gridToCurvilinear(int gridID);

/* GME grid */
struct cart {
  double x[3];
};

struct geo {
  double lon;
  double lat;
};

void correct_sinxvals(int xsize, int ysize, double *xvals);

double areas(struct cart *dv1, struct cart *dv2, struct cart *dv3);
struct cart gc2cc(struct geo *position);
void factorni(int kni, int *kni2, int *kni3);
void gme_grid_restore(double *p, int ni, int nd);
void gme_grid(int gridsize, double *rlon, double *rlat,
	      double *blon, double *blat, int *imask,
              int ni, int nd, int ni2, int ni3);

/* Rotated grid */
double lamrot_to_lam(double phis, double rlas, double polphi, double pollam, double polgam);
double phirot_to_phi(double phis, double rlas, double polphi, double polgam);
double rl_to_rls(double phi, double rla, double polphi, double pollam);
double ph_to_phs(double phi, double rla, double polphi, double pollam);
void usvs_to_uv(double us, double vs, double phi, double rla,
		double polphi, double pollam, double *u, double *v);


// Projection codes for proj_info structure:
#define PROJ_LATLON  0
#define PROJ_MERC    1
#define PROJ_LC      3
#define PROJ_PS      5


typedef struct {
  int       code;     // Integer code for projection type
  double    lat1;     // SW latitude (1,1) in degrees (-90->90N)
  double    lon1;     // SW longitude (1,1) in degrees (-180->180E)
  double    dx;       // Grid spacing in meters at truelats, used
                      // only for ps, lc, and merc projections
  double    dlat;     // Lat increment for lat/lon grids
  double    dlon;     // Lon increment for lat/lon grids
  double    stdlon;   // Longitude parallel to y-axis (-180->180E)
  double    truelat1; // First true latitude (all projections)
  double    truelat2; // Second true lat (LC only)
  double    hemi;     // 1 for NH, -1 for SH
  double    cone;     // Cone factor for LC projections
  double    polei;    // Computed i-location of pole point
  double    polej;    // Computed j-location of pole point
  double    rsw;      // Computed radius to SW corner
  double    rebydx;   // Earth radius divided by dx
  int       init;     // Flag to indicate if this struct is ready for use
} proj_info_t;

/* Lambert Conformal grid (new version) */
void map_set(int proj_code, double lat1, double lon1, double dx, double stdlon,
	     double truelat1, double truelat2, proj_info_t *proj);
void ijll_lc(double i, double j, proj_info_t proj, double *lat, double *lon);

/* Lambert Conformal grid (old version) */
/*
int W3FB12(double xi, double xj, double alat1, double elon1, double dx,
	   double elonv, double alatan, double *alat, double *elon);
*/
#endif  /* _GRID_H */