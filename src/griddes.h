#ifndef _GRIDDES_H
#define _GRIDDES_H

typedef struct {
  double *xvals;
  double *yvals;
  double *xbounds;
  double *ybounds;
  double *area;
  double  xfirst, yfirst;
  double  xlast, ylast;
  double  xinc, yinc;
  double  xpole, ypole, angle;    /* rotated north pole             */
  double  originLon;              /* lambert                        */
  double  originLat;
  double  lonParY;
  double  lat1;
  double  lat2;
  int     projflag;
  int     scanflag;
  int     def_originLon;
  int     def_originLat;
  int     def_lonParY;
  int     def_lat1;
  int     def_lat2;
  double  a;
  double  lon_0;
  double  lat_0;
  double  lat_1;
  double  lat_2;
  int     def_lon_0;
  int     def_lat_0;
  int     def_lat_1;
  int     def_lat_2;
  int     prec;
  int     isRotated;              /* TRUE for rotated grids         */
  int     type;
  int     ntr;
  int    *rowlon;
  int     nvertex;
  int     size;
  int     xsize;
  int     ysize;
  int     def_xfirst;
  int     def_yfirst;
  int     def_xlast;
  int     def_ylast;
  int     def_xinc;
  int     def_yinc;
  int     nd, ni, ni2, ni3;
  char    xname[128];
  char    xlongname[128];
  char    xunits[128];
  char    yname[128];
  char    ylongname[128];
  char    yunits[128];
}
GRID;

void gridInit(GRID *grid);
int gridDefine(GRID grid);

int gridFromNCfile(const char *gridfile);
int gridFromH5file(const char *gridfile);

#endif  /* _GRIDDES_H */
