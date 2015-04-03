/**
 * @file grid_cell.h
 * @brief Structs and interfaces to handle grid cells
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

#ifndef GRID_CELL_H
#define GRID_CELL_H

enum edge_type {
   GREAT_CIRCLE = 0, //!< great circle
   LAT_CIRCLE   = 1, //!< latitude circle
   LON_CIRCLE   = 2, //!< longitude circle
};

struct grid_cell {
   double * coordinates_x, * coordinates_y;
   enum edge_type * edge_type;
   unsigned num_corners;
};

/**
 * initiates a grid_cell object
 * before the first being used a grid_cell object has to be initialised
 * @param[in] cell object to be initialised
 * @see free_grid_cell
 * @see get_grid_cell
 */
void init_grid_cell(struct grid_cell * cell);

/**
 * copies a given grid cell
 * @param[in]  in_cell  cell to be copied
 * @param[out] out_cell copied cell
 * @remarks out_cell needs to be a cell that has previously been
 *          initialised or a cell that already contains valid data
 */
void copy_grid_cell(struct grid_cell in_cell, struct grid_cell * out_cell);

/**
 * frees all memory associated with a grid_cell object and reinitialised
 * the cell
 * @param[in,out] cell
 */
void free_grid_cell(struct grid_cell * cell);

void pack_grid_cell(struct grid_cell cell, double ** dble_buf,
                    unsigned dble_buf_offset, unsigned * dble_buf_data_size,
                    unsigned * dble_buf_size, unsigned ** uint_buf,
                    unsigned uint_buf_offset, unsigned * uint_buf_data_size,
                    unsigned * uint_buf_size);

void unpack_grid_cell(struct grid_cell * cell, double * dble_buf,
                      unsigned * dble_buf_data_size, unsigned * uint_buf,
                      unsigned * uint_buf_data_size);

#endif // GRID_CELL_H
