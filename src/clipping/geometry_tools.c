/**
 * @file geometry_tools.c
 * @brief Set of functions to work with coordinates
 *
 * Note: Not all functions are documented by Doxygen. See the source code
 * and \ref geometry.h for further details.
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

#include <math.h>

#include "geometry.h"

double get_vector_angle(double a_vector[3], double b_vector[3]) {

   double dot_product = a_vector[0]*b_vector[0]+a_vector[1]*b_vector[1]+a_vector[2]*b_vector[2];

   double angle;

   // the acos most accurate in the range [-0.5;0.5]
   if (fabs(dot_product) <= 0.5) // the range in which the acos is most accurate
      angle = acos(dot_product);
   else {

#define CROSS_PRODUCT3D(out,a,b) \
  out[0]=a[1]*b[2]-a[2]*b[1]; \
  out[1]=a[2]*b[0]-a[0]*b[2]; \
  out[2]=a[0]*b[1]-a[1]*b[0];
#define NORM(a) sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])

      double temp_vector[3];

      CROSS_PRODUCT3D(temp_vector, a_vector, b_vector);

      double asin_tmp = asin(NORM(temp_vector));

      if (dot_product < 0.0) // if the angle is bigger than (PI / 2)
         angle = M_PI - asin_tmp;
      else
         angle = asin_tmp;

#undef NORM
#undef CROSS_PRODUCT3D
   }

   if (angle < 0.0) return 0;
   if (angle > M_PI) return M_PI;
   return angle;
}

  /* Taken from http://www.geoclub.de/viewtopic.php?f=54&t=29689

     Further information:
     http://en.wikipedia.org/wiki/List_of_common_coordinate_transformations */

void LLtoXYZ(double lon, double lat, double p_out[]) {

   double cos_lat = cos(lat);
   p_out[0] = cos_lat * cos(lon);
   p_out[1] = cos_lat * sin(lon);
   p_out[2] = sin(lat);
}

void LLtoXYZ_deg(double lon, double lat, double p_out[]) {

   lon *= rad;
   lat *= rad;

   LLtoXYZ(lon, lat, p_out);
}

/** \brief compute longitude angle between two points in radian 
**/
double get_angle (double a_lon, double b_lon) {

   while (a_lon - b_lon >   M_PI) b_lon += 2 * M_PI;
   while (a_lon - b_lon < - M_PI) b_lon -= 2 * M_PI;

   return a_lon - b_lon;
}

void XYZtoLL (double p_in[], double * lon, double * lat) {

/*  convert from cartesian to spherical coordinates
    taken from:
    http://www.geoclub.de/viewtopic.php?f=54&t=29689

    Further information:
    http://en.wikipedia.org/wiki/List_of_common_coordinate_transformations */

   *lon = atan2(p_in[1] , p_in[0]);
   *lat = M_PI_2 - acos(p_in[2]);
}
