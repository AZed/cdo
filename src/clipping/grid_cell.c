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

#include "grid_cell.h"
#include "utils.h"
#include "ensure_array_size.h"

void init_grid_cell(struct grid_cell * cell) {

   cell->coordinates_x = NULL;
   cell->coordinates_y = NULL;
   cell->edge_type = NULL;
   cell->num_corners = 0;
}

void copy_grid_cell(struct grid_cell in_cell, struct grid_cell * out_cell) {

   int flag;

   flag = ((out_cell->coordinates_x == NULL) << 0) +
          ((out_cell->coordinates_y == NULL) << 1) +
          ((out_cell->edge_type == NULL)     << 2);

   if (flag != 0 && flag != 7)
      abort_message("ERROR: inconsistent grid_cell data structure\n",
                    __FILE__, __LINE__);

   if ((out_cell->coordinates_x == NULL) ||
       (in_cell.num_corners > out_cell->num_corners)) {

      free(out_cell->coordinates_x);
      free(out_cell->coordinates_y);
      free(out_cell->edge_type);
      out_cell->coordinates_x = malloc(in_cell.num_corners *
                                       sizeof(*(out_cell->coordinates_x)));
      out_cell->coordinates_y = malloc(in_cell.num_corners *
                                       sizeof(*(out_cell->coordinates_y)));
      out_cell->edge_type = malloc(in_cell.num_corners *
                                   sizeof(*(out_cell->edge_type)));
   }

   memcpy(out_cell->coordinates_x, in_cell.coordinates_x,
          in_cell.num_corners * sizeof(*(out_cell->coordinates_x)));
   memcpy(out_cell->coordinates_y, in_cell.coordinates_y,
          in_cell.num_corners * sizeof(*(out_cell->coordinates_y)));
   memcpy(out_cell->edge_type, in_cell.edge_type,
          in_cell.num_corners * sizeof(*(out_cell->edge_type)));
   out_cell->num_corners = in_cell.num_corners;
}

void free_grid_cell(struct grid_cell * cell) {

   if (cell->coordinates_x != NULL) free(cell->coordinates_x);
   if (cell->coordinates_y != NULL) free(cell->coordinates_y);
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

   if (num_corners != cell->num_corners) {
      cell->coordinates_x = realloc (cell->coordinates_x, num_corners * sizeof(cell->coordinates_x[0]));
      cell->coordinates_y = realloc (cell->coordinates_y, num_corners * sizeof(cell->coordinates_y[0]));
      cell->edge_type = realloc (cell->edge_type, num_corners * sizeof(cell->edge_type[0]));
   }

   cell->num_corners = num_corners;
   memcpy(cell->coordinates_x, dble_buf, num_corners * sizeof(double));
   memcpy(cell->coordinates_y, dble_buf+num_corners, num_corners * sizeof(double));

   unsigned i;
   for (i = 1; i <= num_corners; ++i)
     cell->edge_type[i-1] = (enum edge_type)uint_buf[i];
}
