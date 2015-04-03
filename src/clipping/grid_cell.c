/**
 * @file grid_cell.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 *             Thomas Jahns <jahns@dkrz.de>
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
#include <string.h>
#include <stdio.h>

#include "grid_cell.h"
#include "utils.h"
#include "ensure_array_size.h"
#include "geometry.h"

void init_grid_cell(struct grid_cell * cell) {

   cell->coordinates_x = NULL;
   cell->coordinates_y = NULL;
   cell->coordinates_xyz = NULL;
   cell->edge_type = NULL;
   cell->num_corners = 0;
   cell->array_size = 0;
}

void copy_grid_cell(struct grid_cell in_cell, struct grid_cell * out_cell) {

   if (in_cell.num_corners > out_cell->array_size) {

      free(out_cell->coordinates_x);
      free(out_cell->coordinates_y);
      free(out_cell->coordinates_xyz);
      free(out_cell->edge_type);
      out_cell->coordinates_x = malloc(in_cell.num_corners *
                                       sizeof(*(out_cell->coordinates_x)));
      out_cell->coordinates_y = malloc(in_cell.num_corners *
                                       sizeof(*(out_cell->coordinates_y)));
      out_cell->coordinates_xyz = malloc(3 * in_cell.num_corners *
                                         sizeof(*(out_cell->coordinates_xyz)));
      out_cell->edge_type = malloc(in_cell.num_corners *
                                   sizeof(*(out_cell->edge_type)));
      out_cell->array_size = in_cell.num_corners;
   }

   memcpy(out_cell->coordinates_x, in_cell.coordinates_x,
          in_cell.num_corners * sizeof(*(out_cell->coordinates_x)));
   memcpy(out_cell->coordinates_y, in_cell.coordinates_y,
          in_cell.num_corners * sizeof(*(out_cell->coordinates_y)));
   memcpy(out_cell->coordinates_xyz, in_cell.coordinates_xyz,
          3 * in_cell.num_corners * sizeof(*(out_cell->coordinates_xyz)));
   memcpy(out_cell->edge_type, in_cell.edge_type,
          in_cell.num_corners * sizeof(*(out_cell->edge_type)));
   out_cell->num_corners = in_cell.num_corners;
}

void free_grid_cell(struct grid_cell * cell) {

   if (cell->coordinates_x != NULL) free(cell->coordinates_x);
   if (cell->coordinates_y != NULL) free(cell->coordinates_y);
   if (cell->coordinates_xyz != NULL) free(cell->coordinates_xyz);
   if (cell->edge_type != NULL) free(cell->edge_type);

   init_grid_cell(cell);
}

void pack_grid_cell(struct grid_cell cell, double ** dble_buf,
                    unsigned dble_buf_offset, unsigned * dble_buf_data_size,
                    unsigned * dble_buf_size, unsigned ** uint_buf,
                    unsigned uint_buf_offset, unsigned * uint_buf_data_size,
                    unsigned * uint_buf_size) {

   unsigned required_dble_buf_size, required_uint_buf_size;

   required_dble_buf_size = 2 * cell.num_corners;
   required_uint_buf_size = cell.num_corners + 1;

   ENSURE_ARRAY_SIZE(*dble_buf, *dble_buf_size, dble_buf_offset+required_dble_buf_size);
   ENSURE_ARRAY_SIZE(*uint_buf, *uint_buf_size, uint_buf_offset+required_uint_buf_size);

   memcpy((*dble_buf)+dble_buf_offset, cell.coordinates_x, cell.num_corners * sizeof(double));
   memcpy((*dble_buf)+dble_buf_offset+cell.num_corners, cell.coordinates_y,
          cell.num_corners * sizeof(double));

   *dble_buf_data_size = required_dble_buf_size;

   (*uint_buf)[uint_buf_offset] = cell.num_corners;

   unsigned i;
   for (i = 1; i <= cell.num_corners; ++i)
      (*uint_buf)[uint_buf_offset+i] = cell.edge_type[i-1];

   *uint_buf_data_size = required_uint_buf_size;
}

void unpack_grid_cell(struct grid_cell * cell, double * dble_buf,
                      unsigned * dble_buf_data_size, unsigned * uint_buf,
                      unsigned * uint_buf_data_size) {

   unsigned num_corners;

   num_corners = uint_buf[0];

   *dble_buf_data_size = 2 * num_corners;
   *uint_buf_data_size = num_corners + 1;

   if (num_corners > cell->array_size) {
      cell->coordinates_x = realloc(cell->coordinates_x,
                                    num_corners *
                                    sizeof(cell->coordinates_x[0]));
      cell->coordinates_y = realloc(cell->coordinates_y,
                                    num_corners *
                                    sizeof(cell->coordinates_y[0]));
      cell->coordinates_xyz = realloc(cell->coordinates_xyz,
                                      3 * num_corners *
                                      sizeof(cell->coordinates_xyz[0]));
      cell->edge_type = realloc(cell->edge_type,
                                num_corners * sizeof(cell->edge_type[0]));
      cell->array_size = num_corners;
   }

   cell->num_corners = num_corners;
   memcpy(cell->coordinates_x, dble_buf, num_corners * sizeof(double));
   memcpy(cell->coordinates_y, dble_buf+num_corners, num_corners * sizeof(double));

   unsigned i;
   for (i = 1; i <= num_corners; ++i) {
     cell->edge_type[i-1] = (enum edge_type)uint_buf[i];
     LLtoXYZ(cell->coordinates_x[i], cell->coordinates_y[i],
             cell->coordinates_xyz + 3*i);
   }
}

void print_grid_cell(FILE * stream, struct grid_cell cell, char * name) {

  char * out = NULL;
  unsigned out_array_size = 0;
  unsigned out_size = 0;

  if (name != NULL) {

    out_size = strlen(name) + 1 + 1 + 1;
    ENSURE_ARRAY_SIZE(out, out_array_size, out_size);

    strcpy(out, name);
    strcat(out, ":\n");
  }

  for (unsigned i = 0; i < cell.num_corners; ++i) {

    char buffer[1024];

    sprintf(buffer, "%d x %.16f y %.16f %s\n", i, cell.coordinates_x[i],
           cell.coordinates_y[i],
           (cell.edge_type[i] == LAT_CIRCLE)?("LAT_CIRCLE"):
           ((cell.edge_type[i] == LON_CIRCLE)?("LON_CIRCLE"):
                                              ("GREAT_CIRCLE")));

    out_size += strlen(buffer);

    ENSURE_ARRAY_SIZE(out, out_array_size, out_size);

    strcat(out, buffer);
  }

  if (out != NULL)
    fputs(out, stream);

  free(out);
}
