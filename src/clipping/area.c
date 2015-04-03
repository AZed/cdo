/**
 * @file area.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 * URL: https://redmine.dkrz.de/doc/YAC/html/index.html
 *
 * This file is part of YAC.
 *
 * YAC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with YAC.  If not, see <http://www.gnu.org/licenses/gpl.txt>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "area.h"
#include "clipping.h"
#include "geometry.h"
#include "utils.h"
#include "ensure_array_size.h"

// area tolerance (100 m*m = 0.0001 km*km)
const double area_tol = 1e-4;

static double scalar_product(double a[], double b[]);

static void cross_product(double a[], double b[], double cross[]);

static double inner_angle ( double plat, double plon, double qlon, double qlat );

static double partial_area ( double a_lon, double a_lat,
                             double b_lon, double b_lat,
                             double c_lon, double c_lat );

/* ----------------------------------- */

double cell_approx_area ( struct grid_cell cell ) {

  /* adopted from Robert.G. Chamberlain and William.H. Duquette */

  int m, M;
  double SUM = 0.0;
  double area;

  M = cell.num_corners;

  for ( m = 0; m < M; m++ )
    SUM += ( cell.coordinates_x[(m+2)%M] - cell.coordinates_x[m] ) *
             sin(cell.coordinates_y[(m+1)%M]);

  area = -EarthRadius2 * SUM;

  return area;
}

/* ----------------------------------- */

double triangle_area ( struct grid_cell cell ) {

  /* taken from the ICON code, mo_base_geometry.f90
     provided by Luis Kornblueh, MPI-M. */

  double s01, s12, s20;
  double ca1, ca2, ca3;
  double a1, a2, a3;

  double triangle[3][3];
  double u01[3], u12[3], u20[3];

  if ( cell.num_corners != 3 ) {
    printf ("Only for triangles!\n");
    return -1;
  }

  /* Convert into cartesian coordinates */

  for (int m = 0; m < 3; m++ )
    LLtoXYZ(cell.coordinates_x[m], cell.coordinates_y[m], triangle[m]);

  /* First, compute cross products Uij = Vi x Vj. */

  cross_product(triangle[0], triangle[1], u01);
  cross_product(triangle[1], triangle[2], u12);
  cross_product(triangle[2], triangle[0], u20);

  /*  Normalize Uij to unit vectors. */

  s01 = scalar_product(u01, u01);
  s12 = scalar_product(u12, u12);
  s20 = scalar_product(u20, u20);

  /* Test for a degenerated triangle associated with collinear vertices. */

  if ( s01 == 0.0 &&
       s12 == 0.0 &&
       s20 == 0.0 )
    return 0.0;

  s01 = sqrt(s01);
  s12 = sqrt(s12);
  s20 = sqrt(s20);

  for (int m = 0; m < 3; m++ ) {
    u01[m] = u01[m]/s01;
    u12[m] = u12[m]/s12;
    u20[m] = u20[m]/s20;
  }

  /*  Compute interior angles Ai as the dihedral angles between planes:

      CA1 = cos(A1) = -<U01,U20>
      CA2 = cos(A2) = -<U12,U01>
      CA3 = cos(A3) = -<U20,U12>

  */

  ca1 = -u01[0]*u20[0]-u01[1]*u20[1]-u01[2]*u20[2];
  ca2 = -u12[0]*u01[0]-u12[1]*u01[1]-u12[2]*u01[2];
  ca3 = -u20[0]*u12[0]-u20[1]*u12[1]-u20[2]*u12[2];

  if ( ca1 < -1.0 ) ca1 = -1.0;
  if ( ca1 >  1.0 ) ca1 =  1.0;
  if ( ca2 < -1.0 ) ca2 = -1.0;
  if ( ca2 >  1.0 ) ca2 =  1.0;
  if ( ca3 < -1.0 ) ca3 = -1.0;
  if ( ca3 >  1.0 ) ca3 =  1.0;

  a1 = acos(ca1);
  a2 = acos(ca2);
  a3 = acos(ca3);

  /*  Compute areas = a1 + a2 + a3 - pi.

      here for a unit sphere: */

  return MAX(( a1+a2+a3-M_PI ) * EarthRadius * EarthRadius, 0.0);
}

/* ----------------------------------- */

double cell_area ( struct grid_cell cell ) {

  /* generalised version based on the ICON code, mo_base_geometry.f90
     provided by Luis Kornblueh, MPI-M. */

  int const M = cell.num_corners; // number of vertices

  double area;
  double s[cell.num_corners];
  double ca[cell.num_corners];
  double a[cell.num_corners];

  double p[cell.num_corners][3];
  double u[cell.num_corners][3];

  /* Convert into cartesian coordinates */

  for (int m = 0; m < M; m++ )
    LLtoXYZ( cell.coordinates_x[m], cell.coordinates_y[m], p[m]);

  /* First, compute cross products Uij = Vi x Vj. */

  for (int m = 0; m < M; m++ )
    cross_product (p[m], p[(m+1)%M], u[m]);

  /*  Normalize Uij to unit vectors. */

  area = 0.0;

  for (int m = 0; m < M; m++ ) {
    s[m] = scalar_product(u[m], u[m]);
    area += s[m];
  }

  /* Test for a degenerated cells associated with collinear vertices. */

  if ( area != 0.0 ) {

    for (int m = 0; m < M; m++ )
      s[m] = sqrt(s[m]);

    for (int m = 0; m < M; m++ )
      for (int i = 0; i < 3; i++ )
        u[m][i] = u[m][i]/s[m];

    /*  Compute interior angles Ai as the dihedral angles between planes
        by using the definition of the scalar product

                    ab = |a| |b| cos (phi)

        As a and b are already normalised this reduces to

                    ab = cos (phi)

        There is no explanation so far for the - in the loop below.
        But otherwise we don't get the correct results for triangles
        and cells. Must have something to do with the theorem.

     */

    for (int m = 0; m < M; m++ ) {
      ca[m] = - scalar_product(u[m], u[(m+1)%M]);
      if ( ca[m] < -1.0 ) ca[m] = -1.0;
      if ( ca[m] >  1.0 ) ca[m] =  1.0;
      a[m] = acos(ca[m]);
    }

    /*  Compute areas = a1 + a2 + a3 - (M-2) * pi.

        here for a unit sphere: */

    area = - (double) (M-2) * M_PI;

    for (int m = 0; m < M; m++ )
      area += a[m];
  }

  return MAX(area * EarthRadius * EarthRadius, 0.0);
}

/* ----------------------------------- */

double girards_area ( struct grid_cell cell  ) {

  /* Bevis and Cambareri, 1987

     this algorithm provides wrong results if one
     of the inner angles becomes less than 1e-18 rad

     (R. Redler, M. Hanke 2013)
  */

  int m;
  const double tol = 1e-18;
  double area = 0.0;

  int M = cell.num_corners;
  if (M < 3) return area;  // a degenerate cell

  double * theta = malloc ( M * sizeof(theta[0]) );

  for ( m = 0; m < M; m++ ) {
     theta[m] = partial_area(cell.coordinates_x[(m+1)%M], cell.coordinates_y[(m+1)%M],
                             cell.coordinates_x[(m+2)%M], cell.coordinates_y[(m+2)%M],
                             cell.coordinates_x[m%M], cell.coordinates_y[m%M]);
     if ( theta[m] < tol ) {
       area = cell3d_area ( cell );
       return area;
     }
  }

  /*  Sum up partial areas */

  area = - (double) (M-2) * M_PI;

  for ( m = 0; m < M; m++ )
    area += theta[m];

  /* Area on Sphere with radius EarthRadius */

  area *= EarthRadius * EarthRadius;

  free ( theta );

  if ( area <= 0.0 ) area = cell3d_area ( cell );

  return area;
}

/* ----------------------------------- */

double cell3d_area( struct grid_cell cell ) {

   /* http://geomalgorithms.com/a01-_area.html */

   int const M = cell.num_corners;
   double area = 0.0;

   if (M < 3) return 0.0;  // a degenerated cell

   double an, ax, ay, az; // abs value of normal and its coords

   int  coord;           // coord to ignore: 1=x[1], 2=x[2], 3=x[3]

   double V[cell.num_corners][3];
   double Norm[3];
   double edge[2][3];

   /* transform vertices into cartesian coordinates on the earth surface
      must be done already here to avoid round-off errors later */
   for (int i = 0; i < M; i++ )
      LLtoXYZ( cell.coordinates_x[i], cell.coordinates_y[i], V[i]);

   /* compute normal vector */
   edge[0][0] = V[0][0] - V[1][0];
   edge[0][1] = V[0][1] - V[1][1];
   edge[0][2] = V[0][2] - V[1][2];
   edge[1][0] = V[0][0] - V[2][0];
   edge[1][1] = V[0][1] - V[2][1];
   edge[1][2] = V[0][2] - V[2][2];

   cross_product(edge[0], edge[1], Norm);

   /* select largest abs coordinate to ignore for projection */
   ax = (Norm[0]>0 ? Norm[0] : -Norm[0]);    // abs x-coord
   ay = (Norm[1]>0 ? Norm[1] : -Norm[1]);    // abs y-coord
   az = (Norm[2]>0 ? Norm[2] : -Norm[2]);    // abs z-coord

   coord = 3;                    // ignore z-coord x[2]
   if (ax > ay) {
       if (ax > az) coord = 1;   // ignore x-coord x[0]
   }
   else if (ay > az) coord = 2;  // ignore y-coord x[1]

   /* compute area of the 2D projection

      fabs is used to remain independent from
      cyclic or anticyclic ordering of vertices */

   for (int i = 1; i<M; i++) {
        switch (coord) {
          case 1:
            area += fabs(V[i%M][1] * (V[(i+1)%M][2] - V[(i-1)%M][2]));
            continue;
          case 2:
            area += fabs(V[i%M][0] * (V[(i+1)%M][2] - V[(i-1)%M][2]));
            continue;
          case 3:
            area += fabs(V[i%M][0] * (V[(i+1)%M][1] - V[(i-1)%M][1]));
            continue;
        }
    }
    switch (coord) {    // wrap-around term
      case 1:
        area += fabs(V[0][1] * (V[1][2] - V[M-1][2]));
        break;
      case 2:
        area += fabs(V[0][0] * (V[1][2] - V[M-1][2]));
        break;
      case 3:
        area += fabs(V[0][0] * (V[1][1] - V[M-1][1]));
        break;
    }

    /* scale to get area before projection */
    an = sqrt( ax*ax + ay*ay + az*az); // length of normal vector
    switch (coord) {
      case 1:
        if ( ax > 0 )
                area *= (an / (2*ax));
        else
          area = 0.0;
              break;
      case 2:
        if ( ay > 0 )
          area *= (an / (2*ay));
        else
          area = 0.0;
              break;
      case 3:
        if ( az > 0 )
          area *= (an / (2*az));
        else
          area = 0.0;
        break;
    }

    return area * EarthRadius * EarthRadius;
}

/** area of a spherical triangle based on L'Huilier's Theorem
  *
  * source code is taken from code by Robert Oehmke of Earth System Modeling
  * Framework (www.earthsystemmodeling.org)
  *
  * it has been extended by a more accurate computation of vector angles
  *
  * the license statement for this routine is as follows:
  * Earth System Modeling Framework
  * Copyright 2002-2013, University Corporation for Atmospheric Research,
  * Massachusetts Institute of Technology, Geophysical Fluid Dynamics
  * Laboratory, University of Michigan, National Centers for Environmental
  * Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
  * NASA Goddard Space Flight Center.
  * Licensed under the University of Illinois-NCSA License.
  */
static double
tri_area(double u[3], double v[3], double w[3]) {

  double a = get_vector_angle(u,v);
  double b = get_vector_angle(u,w);
  double c = get_vector_angle(w,v);

  double s=0.5*(a+b+c);

  double t = tan ( s / 2.0 ) * tan ( ( s - a ) / 2.0 ) *
             tan ( ( s - b ) / 2.0 ) * tan ( ( s - c ) / 2.0 );

  return fabs ( 4.0 * atan ( sqrt (fabs ( t ) ) ) );;
}

double pole_area ( struct grid_cell cell ) {

  double area = 0.0;

  int M = cell.num_corners;

  if (M < 2) return 0.0;

  int closer_to_south_pole = cell.coordinates_y[0] < 0;

  double pole_vec[3] = {0, 0, (closer_to_south_pole)?-1:1};

  // it would also be possible to use the equator instead
  // of the poles as the baseline
  // this could be used as special case for cell close
  // the equator (were the other method is probably most
  // inaccurate)

  for (int i = 0; i < M; ++i) {

    // if one of the points it at the pole
    if (fabs(fabs(cell.coordinates_y[i]) - M_PI_2) < 1e-12) continue;
    if (fabs(fabs(cell.coordinates_y[(i+1)%M]) - M_PI_2) < 1e-12) {
      ++i; // we can skip the next edge
      continue;
    }

    if (cell.edge_type[i] == GREAT_CIRCLE || cell.edge_type[i] == LON_CIRCLE) {

      double a[3];
      double b[3];

      LLtoXYZ(cell.coordinates_x[i], cell.coordinates_y[i], a);
      LLtoXYZ(cell.coordinates_x[(i+1)%M], cell.coordinates_y[(i+1)%M], b);

      double edge_direction = a[0]*b[1]-a[1]*b[0]; // 3. component of cross product

      // if the edge is nearly on a circle of longitude
      if (fabs(edge_direction) < 1e-12) continue;

      double tmp_area = tri_area(a, b, pole_vec);

      // or the other way round
      if (edge_direction > 0) area -= tmp_area;
      else                    area += tmp_area;

    } else if (cell.edge_type[i] == LAT_CIRCLE) {

      // the area of a sphere cap is:
      // A = 2 * PI * r * h (where r == 1 and h == 1 - cos(d_lat))
      // scaled with the longitude angle this is:
      // A' = (d_lon / (2 * PI)) * A
      // => A' = d_lon * (1 - cos(d_lat))

      double d_lon = get_angle(cell.coordinates_x[i], cell.coordinates_x[(i+1)%M]);
      double d_lat = M_PI_2;

      if (closer_to_south_pole)
        d_lat += cell.coordinates_y[i];
      else
        d_lat -= cell.coordinates_y[i];

      double h = 1 - cos(d_lat);

      area += d_lon * h;

    } else {

      abort_message("ERROR: unsupported edge type\n", __FILE__, __LINE__);
    }
  }
  return fabs(area * EarthRadius * EarthRadius);
}

static double
lat_edge_correction(double a[3], double b[3], double lon_a, double lon_b) {

  double const tol = 1e-8;

  if (fabs(a[2] - b[2]) > tol)
    abort_message("ERROR: latitude of both corners is not identical\n",
                  __FILE__, __LINE__);

  double h = fabs(a[2]);

  // if we are at the equator or at a pole
  if (h < tol || fabs(1.0 - h) < tol)
    return 0.0;

  double lat_area = fabs((1.0 - h) * get_angle(lon_a, lon_b));

  double pole[3] = {0, 0, (a[2] >= 0.0)?1:-1};
  double gc_area = tri_area(a, b, pole);

  double correction = lat_area - gc_area;

  if (correction < 0) return 0;
  else return correction;
}

 /*
  * source code is originally based on code by Robert Oehmke of Earth System Modeling
  * Framework (www.earthsystemmodeling.org)
  *
  * it has be extended to support YAC data structures and two types of
  * grid cell edges (great circle and circle of latitude)
  */
double huiliers_area(struct grid_cell cell) {

  double tmp_points[cell.num_corners][3];
  double const tol = 1e-8;

  if (cell.num_corners < 2) return 0;

  // convert lon-lat to xyz
  for (int i = 0; i < cell.num_corners; i++)
    LLtoXYZ(cell.coordinates_x[i], cell.coordinates_y[i], tmp_points[i]);

  // sum areas around cell
  double sum = 0.0;

  for (int i = 2; i < cell.num_corners; i++)
    sum += tri_area(tmp_points[0], tmp_points[i-1], tmp_points[i]);

  // check for edges of latitude
  unsigned num_lat_circle_edges = 0;
  for (unsigned i = 0; i < cell.num_corners; ++i)
    if (cell.edge_type[i] == LAT_CIRCLE)
      num_lat_circle_edges++;

  if (num_lat_circle_edges > 2)
    abort_message("ERROR: invalid cell (has more than two edges that are "
                  "latitude circles)\n", __FILE__, __LINE__);

  if (num_lat_circle_edges > 0) {

    // compute minimum and maximum height of cell
    double min = cell.coordinates_y[0], max = cell.coordinates_y[0];

    for (int i = 1; i < cell.num_corners; ++i)
      if (cell.coordinates_y[i] < min) min = cell.coordinates_y[i];
      else if (cell.coordinates_y[i] > max) max = cell.coordinates_y[i];

    double min_factor = 1.0;

    if (max <= 0.0) {

      double tmp = max;
      max = min;
      min = tmp;

    } else if (signbit(min) != signbit(max))
      min_factor = -1.0;

    for (int i = 0; i < cell.num_corners; ++i) {

      if (cell.edge_type[i] == LAT_CIRCLE) {

        double correction =
          lat_edge_correction(tmp_points[i],
                              tmp_points[(i+1)%cell.num_corners],
                              cell.coordinates_x[i],
                              cell.coordinates_x[(i+1)%cell.num_corners]);

        if (fabs(cell.coordinates_y[i] - min) < tol)
          sum += min_factor * correction;
        else if (fabs(cell.coordinates_y[i] - max) < tol)
          sum -= correction;
        else
          abort_message("ERROR: internal error...should have not occured\n",
                        __FILE__, __LINE__);
      }
    }
  }

  // return area
  return sum * EarthRadius * EarthRadius;
}

/* ----------------------------------- */

double partial_area ( double a_lon, double a_lat,
                      double b_lon, double b_lat,
                      double c_lon, double c_lat ) {

  double theta;
  double angle_f;
  double angle_b;

  angle_f = inner_angle ( a_lat, a_lon, b_lat, b_lon );
  angle_b = inner_angle ( a_lat, a_lon, c_lat, c_lon );

  theta = angle_b - angle_f;

  if ( theta < 0.0 ) theta = theta + M_PI + M_PI;

  return theta;
}

/* ----------------------------------- */

double inner_angle ( double plat, double plon, double qlat, double qlon ) {

  double t = sin((qlon-plon))*cos(qlat);

  double b = sin(qlat)*cos(plat)
           - cos(qlat)*sin(plat) * cos((qlon-plon));

  return atan2(b,t);
}

/* ----------------------------------- */

static double scalar_product(double a[], double b[]) {
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

/* ----------------------------------- */

static void cross_product(double a[], double b[], double cross[]) {
  cross[0] = a[1]*b[2] - a[2]*b[1];
  cross[1] = a[2]*b[0] - a[0]*b[2];
  cross[2] = a[0]*b[1] - a[1]*b[0];
}
