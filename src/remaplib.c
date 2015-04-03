/*
  This is a C library of the Fortran SCRIP version 1.4

  ===>>> Please send bug reports to Uwe.Schulzweida@zmaw.de <<<===

  Spherical Coordinate Remapping and Interpolation Package (SCRIP)
  ================================================================

  SCRIP is a software package which computes addresses and weights for
  remapping and interpolating fields between grids in spherical coordinates.
  It was written originally for remapping fields to other grids in a coupled
  climate model, but is sufficiently general that it can be used in other 
  applications as well. The package should work for any grid on the surface
  of a sphere. SCRIP currently supports four remapping options:

  Conservative remapping
  ----------------------
  First- and second-order conservative remapping as described in
  Jones (1999, Monthly Weather Review, 127, 2204-2210).

  Bilinear interpolation
  ----------------------
  Slightly generalized to use a local bilinear approximation
  (only logically-rectangular grids).

  Bicubic interpolation
  ----------------------
  Similarly generalized (only logically-rectangular grids).

  Distance-weighted averaging
  ---------------------------
  Distance-weighted average of a user-specified number of nearest neighbor values.

  Documentation
  =============

  http://climate.lanl.gov/Software/SCRIP/SCRIPusers.pdf

*/
/*
  2013-11-08 Uwe Schulzweida: split remapgrid class to src_grid and tgt_grid
  2012-01-16 Uwe Schulzweida: alloc grid2_bound_box only for conservative remapping
  2011-01-07 Uwe Schulzweida: Changed remap weights from 2D to 1D array
  2009-05-25 Uwe Schulzweida: Changed restrict data type from double to int
  2009-01-11 Uwe Schulzweida: OpenMP parallelization
 */

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <limits.h>
#include <time.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"
#include "util.h"  /* progressStatus */


#define IS_REG2D_GRID(gridID)  (!gridIsRotated(gridID) && (gridInqType(gridID) == GRID_LONLAT || gridInqType(gridID) == GRID_GAUSSIAN))

/* used for store_link_fast */

#define BLK_SIZE 4096
#define BLK_NUM(x) (x/grid_store->blk_size)
#define BLK_IDX(x) (x%grid_store->blk_size)
struct grid_layer
{
  int *grid2_link;
  struct grid_layer *next;
};

typedef struct grid_layer grid_layer_t;

typedef struct
{
  int blk_size;
  int max_size;
  int nblocks;
  int *blksize;
  int *nlayers;
  grid_layer_t **layers;
} grid_store_t;


/* constants */

/* #define  BABY_STEP  0.001 */ /* original value */
#define  BABY_STEP  0.001

#define  ZERO     0.0
#define  ONE      1.0
#define  TWO      2.0
#define  THREE    3.0
#define  HALF     0.5
#define  QUART    0.25
#define  BIGNUM   1.e+20
#define  TINY     1.e-14


#define  REMAP_GRID_TYPE_REG2D     1
#define  REMAP_GRID_TYPE_CURVE2D   2
#define  REMAP_GRID_TYPE_UNSTRUCT  3


static int  remap_store_link_fast = TRUE;
static int  remap_write_remap     = FALSE;
static int  remap_num_srch_bins   = 180;
#define  DEFAULT_MAX_ITER  100
static long remap_max_iter        = DEFAULT_MAX_ITER;  /* Max iteration count for i, j iteration */

void remap_set_int(int remapvar, int value)
{
  if      ( remapvar == REMAP_STORE_LINK_FAST ) remap_store_link_fast = value;
  else if ( remapvar == REMAP_WRITE_REMAP     ) remap_write_remap     = value;
  else if ( remapvar == REMAP_MAX_ITER        ) remap_max_iter        = value;
  else if ( remapvar == REMAP_NUM_SRCH_BINS   ) remap_num_srch_bins   = value;
  else      cdoAbort("Unsupported remap variable (%d)!", remapvar);
}

/* static double north_thresh =  1.45;  */ /* threshold for coord transformation */
/* static double south_thresh = -2.00;  */ /* threshold for coord transformation */
static double north_thresh =  2.00;  /* threshold for coord transformation */
static double south_thresh = -2.00;  /* threshold for coord transformation */

double intlin(double x, double y1, double x1, double y2, double x2);

extern int timer_remap, timer_remap_con, timer_remap_con_l1, timer_remap_con_l2;
extern int timer_remap_bil, timer_remap_nn;


void remapGridFree(remapgrid_t *grid)
{
  if ( grid->vgpm ) free(grid->vgpm);
  if ( grid->mask ) free(grid->mask);

  if ( grid->reg2d_center_lat ) free(grid->reg2d_center_lat);
  if ( grid->reg2d_center_lon ) free(grid->reg2d_center_lon);
  if ( grid->reg2d_corner_lat ) free(grid->reg2d_corner_lat);
  if ( grid->reg2d_corner_lon ) free(grid->reg2d_corner_lon);

  if ( grid->cell_center_lat ) free(grid->cell_center_lat);
  if ( grid->cell_center_lon ) free(grid->cell_center_lon);
  if ( grid->cell_corner_lat ) free(grid->cell_corner_lat);
  if ( grid->cell_corner_lon ) free(grid->cell_corner_lon);

  if ( grid->cell_area ) free(grid->cell_area);
  if ( grid->cell_frac ) free(grid->cell_frac);

  if ( grid->cell_bound_box ) free(grid->cell_bound_box);

  if ( grid->bin_addr ) free(grid->bin_addr);
  if ( grid->bin_lats ) free(grid->bin_lats);

} /* remapGridFree */

/*****************************************************************************/

void remapVarsFree(remapvars_t *rv)
{
  long i, num_blks;

  if ( rv->pinit == TRUE )
    {
      rv->pinit = FALSE;

      rv->sort_add = FALSE;

      free(rv->src_grid_add);
      free(rv->tgt_grid_add);
      free(rv->wts);

      if ( rv->links.option == TRUE )
	{
	  rv->links.option = FALSE;

	  if ( rv->links.num_blks )
	    {
	      free(rv->links.num_links);
	      num_blks = rv->links.num_blks;
	      for ( i = 0; i < num_blks; ++i )
		{
		  free(rv->links.src_add[i]);
		  free(rv->links.dst_add[i]);
		  free(rv->links.w_index[i]);
		}
	      free(rv->links.src_add);
	      free(rv->links.dst_add);
	      free(rv->links.w_index);
	    }
	}
    }
  else
    fprintf(stderr, "%s Warning: vars not initialized!\n", __func__);

} /* remapVarsFree */

/*****************************************************************************/

void remapgrid_init(remapgrid_t *grid)
{
  grid->remap_grid_type  = -1;
  grid->num_srch_bins    = remap_num_srch_bins; // only for source grid ?

  grid->num_cell_corners = 0;
  grid->luse_cell_corners  = FALSE;
  grid->lneed_cell_corners = FALSE;

  grid->nvgp             = 0;
  grid->vgpm             = NULL;

  grid->mask             = NULL;

  grid->reg2d_center_lon = NULL;
  grid->reg2d_center_lat = NULL;
  grid->reg2d_corner_lon = NULL;
  grid->reg2d_corner_lat = NULL;

  grid->cell_center_lon  = NULL;
  grid->cell_center_lat  = NULL;
  grid->cell_corner_lon  = NULL;
  grid->cell_corner_lat  = NULL;

  grid->cell_area        = NULL;
  grid->cell_frac        = NULL;

  grid->cell_bound_box   = NULL;

  grid->bin_addr         = NULL;
  grid->bin_lats         = NULL;
}

/*****************************************************************************/

void remapGridRealloc(int map_type, remapgrid_t *grid)
{
  long nalloc;

  if ( grid->nvgp )
    grid->vgpm   = realloc(grid->vgpm, grid->nvgp*sizeof(int));

  grid->mask     = realloc(grid->mask, grid->size*sizeof(int));

  if ( remap_write_remap == TRUE || grid->remap_grid_type != REMAP_GRID_TYPE_REG2D )
    {
      grid->cell_center_lon = realloc(grid->cell_center_lon, grid->size*sizeof(double));
      grid->cell_center_lat = realloc(grid->cell_center_lat, grid->size*sizeof(double));
    }

  if ( map_type == MAP_TYPE_CONSERV || map_type == MAP_TYPE_CONSPHERE )
    {
      grid->cell_area = realloc(grid->cell_area, grid->size*sizeof(double));
      memset(grid->cell_area, 0, grid->size*sizeof(double));
    }

  grid->cell_frac = realloc(grid->cell_frac, grid->size*sizeof(double));
  memset(grid->cell_frac, 0, grid->size*sizeof(double));

  if ( grid->lneed_cell_corners )
    {
      if ( grid->num_cell_corners == 0 )
	{
	  cdoAbort("Grid cell corner missing!");
	}
      else
	{
	  nalloc = grid->num_cell_corners*grid->size;

	  grid->cell_corner_lon = realloc(grid->cell_corner_lon, nalloc*sizeof(double));
	  memset(grid->cell_corner_lon, 0, nalloc*sizeof(double));

	  grid->cell_corner_lat = realloc(grid->cell_corner_lat, nalloc*sizeof(double));  
	  memset(grid->cell_corner_lat, 0, nalloc*sizeof(double));
	}
    }
}

/*****************************************************************************/
static
void boundbox_from_corners(long size, long nc, const double *restrict corner_lon,
			   const double *restrict corner_lat, restr_t *restrict bound_box)
{
  long i4, inc, i, j;
  restr_t clon, clat;

#if defined(_OPENMP)
#pragma omp parallel for default(none)        \
  shared(bound_box, corner_lat, corner_lon, nc, size)	\
  private(i4, inc, i, j, clon, clat)
#endif
  for ( i = 0; i < size; ++i )
    {
      i4 = i<<2; // *4
      inc = i*nc;
      clat = RESTR_SCALE(corner_lat[inc]);
      clon = RESTR_SCALE(corner_lon[inc]);
      bound_box[i4  ] = clat;
      bound_box[i4+1] = clat;
      bound_box[i4+2] = clon;
      bound_box[i4+3] = clon;
      for ( j = 1; j < nc; ++j )
	{
	  clat = RESTR_SCALE(corner_lat[inc+j]);
	  clon = RESTR_SCALE(corner_lon[inc+j]);
	  if ( clat < bound_box[i4  ] ) bound_box[i4  ] = clat;
	  if ( clat > bound_box[i4+1] ) bound_box[i4+1] = clat;
	  if ( clon < bound_box[i4+2] ) bound_box[i4+2] = clon;
	  if ( clon > bound_box[i4+3] ) bound_box[i4+3] = clon;
	}
    }
}

static
void boundbox_from_center(int lonIsCyclic, long size, long nx, long ny, const double *restrict center_lon,
			  const double *restrict center_lat, restr_t *restrict bound_box)
{
  long n4, i, j, k, n, ip1, jp1;
  long n_add, e_add, ne_add;
  restr_t tmp_lats[4], tmp_lons[4];  /* temps for computing bounding boxes */

#if defined(_OPENMP)
#pragma omp parallel for default(none)        \
  shared(lonIsCyclic, size, nx, ny, center_lon, center_lat, bound_box)	\
  private(n4, i, j, k, n, ip1, jp1, n_add, e_add, ne_add, tmp_lats, tmp_lons)
#endif
  for ( n = 0; n < size; n++ )
    {
      n4 = n<<2;

      /* Find N,S and NE points to this grid point */
      
      j = n/nx;
      i = n - j*nx;

      if ( i < (nx-1) )
	ip1 = i + 1;
      else
	{
	  /* 2009-01-09 Uwe Schulzweida: bug fix */
	  if ( lonIsCyclic )
	    ip1 = 0;
	  else
	    ip1 = i;
	}

      if ( j < (ny-1) )
	jp1 = j + 1;
      else
	{
	  /* 2008-12-17 Uwe Schulzweida: latitute cyclic ??? (bug fix) */
	  jp1 = j;
	}

      n_add  = jp1*nx + i;
      e_add  = j  *nx + ip1;
      ne_add = jp1*nx + ip1;

      /* Find N,S and NE lat/lon coords and check bounding box */

      tmp_lats[0] = RESTR_SCALE(center_lat[n]);
      tmp_lats[1] = RESTR_SCALE(center_lat[e_add]);
      tmp_lats[2] = RESTR_SCALE(center_lat[ne_add]);
      tmp_lats[3] = RESTR_SCALE(center_lat[n_add]);

      tmp_lons[0] = RESTR_SCALE(center_lon[n]);
      tmp_lons[1] = RESTR_SCALE(center_lon[e_add]);
      tmp_lons[2] = RESTR_SCALE(center_lon[ne_add]);
      tmp_lons[3] = RESTR_SCALE(center_lon[n_add]);

      bound_box[n4  ] = tmp_lats[0];
      bound_box[n4+1] = tmp_lats[0];
      bound_box[n4+2] = tmp_lons[0];
      bound_box[n4+3] = tmp_lons[0];

      for ( k = 1; k < 4; k++ )
	{
	  if ( tmp_lats[k] < bound_box[n4  ] ) bound_box[n4  ] = tmp_lats[k];
	  if ( tmp_lats[k] > bound_box[n4+1] ) bound_box[n4+1] = tmp_lats[k];
	  if ( tmp_lons[k] < bound_box[n4+2] ) bound_box[n4+2] = tmp_lons[k];
	  if ( tmp_lons[k] > bound_box[n4+3] ) bound_box[n4+3] = tmp_lons[k];
	}
    }
}

static
void check_lon_range2(long nc, long nlons, double *corners, double *centers)
{
  long n, k;

  assert(corners != NULL);
  /*
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nlons, lons)
#endif
  */
  for ( n = 0; n < nlons; ++n )
    {
      bool lneg = false, lpos = false;

      for ( k = 0; k < nc; ++k )
	{
 	  if      ( !lneg && corners[n*nc+k] > -PI && corners[n*nc+k] < 0. ) lneg = true;
	  else if ( !lpos && corners[n*nc+k] <  PI && corners[n*nc+k] > 0. ) lpos = true;
	}

      if ( lneg && lpos )
	{
	  if ( centers[n] > PI )
	    for ( k = 0; k < nc; ++k ) corners[n*nc+k] += PI2;
	}
      else
	{
	  for ( k = 0; k < nc; ++k )
	    {
	      if ( corners[n*nc+k] > PI2  ) corners[n*nc+k] -= PI2;
	      if ( corners[n*nc+k] < ZERO ) corners[n*nc+k] += PI2;
	    }
	}
    }
}

static
void check_lon_range(long nlons, double *lons)
{
  long n;

  assert(lons != NULL);

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nlons, lons)
#endif
  for ( n = 0; n < nlons; ++n )
    {
      if ( lons[n] > PI2  ) lons[n] -= PI2;
      if ( lons[n] < ZERO ) lons[n] += PI2;
    }
}

static
void check_lat_range(long nlats, double *lats)
{
  long n;

  assert(lats != NULL);

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nlats, lats)
#endif
  for ( n = 0; n < nlats; ++n )
    {
      if ( lats[n] >  PIH ) lats[n] =  PIH;
      if ( lats[n] < -PIH ) lats[n] = -PIH;
    }
}

static
void check_lon_boundbox_range(long nlons, restr_t *bound_box)
{
  long n, n4;

  assert(bound_box != NULL);

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nlons, bound_box) private(n4)
#endif
  for ( n = 0; n < nlons; ++n )
    {
      n4 = n<<2;
      if ( RESTR_ABS(bound_box[n4+3] - bound_box[n4+2]) > RESTR_SCALE(PI) )
	{
	  bound_box[n4+2] = RESTR_SCALE(0.);
	  bound_box[n4+3] = RESTR_SCALE(PI2);
	}
    }
}

static
void check_lat_boundbox_range(long nlats, restr_t *restrict bound_box, double *restrict lats)
{
  long n, n4;

  assert(bound_box != NULL);

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nlats, bound_box, lats) private(n4)
#endif
  for ( n = 0; n < nlats; ++n )
    {
      n4 = n<<2;
      if ( RESTR_SCALE(lats[n]) < bound_box[n4  ] ) bound_box[n4  ] = RESTR_SCALE(-PIH);
      if ( RESTR_SCALE(lats[n]) > bound_box[n4+1] ) bound_box[n4+1] = RESTR_SCALE( PIH);
    }
}

static
int expand_lonlat_grid(int gridID)
{
  char units[CDI_MAX_NAME];
  int gridIDnew;
  long nx, ny, nxp4, nyp4;
  double *xvals, *yvals;

  nx = gridInqXsize(gridID);
  ny = gridInqYsize(gridID);
  nxp4 = nx+4;
  nyp4 = ny+4;

  xvals = malloc(nxp4*sizeof(double));
  yvals = malloc(nyp4*sizeof(double));
  gridInqXvals(gridID, xvals+2);
  gridInqYvals(gridID, yvals+2);

  gridIDnew = gridCreate(GRID_LONLAT, nxp4*nyp4);
  gridDefXsize(gridIDnew, nxp4);
  gridDefYsize(gridIDnew, nyp4);
	      
  gridInqXunits(gridID,    units);
  gridDefXunits(gridIDnew, units);
  gridInqYunits(gridID,    units);
  gridDefYunits(gridIDnew, units);

  xvals[0] = xvals[2] - 2*gridInqXinc(gridID);
  xvals[1] = xvals[2] - gridInqXinc(gridID);
  xvals[nxp4-2] = xvals[nx+1] + gridInqXinc(gridID);
  xvals[nxp4-1] = xvals[nx+1] + 2*gridInqXinc(gridID);

  yvals[0] = yvals[2] - 2*gridInqYinc(gridID);
  yvals[1] = yvals[2] - gridInqYinc(gridID);
  yvals[nyp4-2] = yvals[ny+1] + gridInqYinc(gridID);
  yvals[nyp4-1] = yvals[ny+1] + 2*gridInqYinc(gridID);

  gridDefXvals(gridIDnew, xvals);
  gridDefYvals(gridIDnew, yvals);

  free(xvals);
  free(yvals);

  if ( gridIsRotated(gridID) )
    {
      gridDefXpole(gridIDnew, gridInqXpole(gridID));
      gridDefYpole(gridIDnew, gridInqYpole(gridID));
    }

  return(gridIDnew);
}

static
int expand_curvilinear_grid(int gridID)
{
  char units[CDI_MAX_NAME];
  int gridIDnew;
  long gridsize, gridsize_new;
  long nx, ny, nxp4, nyp4;
  long i, j;
  double *xvals, *yvals;

  gridsize = gridInqSize(gridID);
  nx = gridInqXsize(gridID);
  ny = gridInqYsize(gridID);
  nxp4 = nx+4;
  nyp4 = ny+4;
  gridsize_new = gridsize + 4*(nx+2) + 4*(ny+2);

  xvals = malloc(gridsize_new*sizeof(double));
  yvals = malloc(gridsize_new*sizeof(double));
  gridInqXvals(gridID, xvals);
  gridInqYvals(gridID, yvals);

  gridIDnew = gridCreate(GRID_CURVILINEAR, nxp4*nyp4);
  gridDefXsize(gridIDnew, nxp4);
  gridDefYsize(gridIDnew, nyp4);

  gridInqXunits(gridID,   units);
  gridDefXunits(gridIDnew, units);
  gridInqYunits(gridID,   units);
  gridDefYunits(gridIDnew, units);

  for ( j = ny-1; j >= 0; j-- )
    for ( i = nx-1; i >= 0; i-- )
      xvals[(j+2)*(nx+4)+i+2] = xvals[j*nx+i];

  for ( j = ny-1; j >= 0; j-- )
    for ( i = nx-1; i >= 0; i-- )
      yvals[(j+2)*(nx+4)+i+2] = yvals[j*nx+i];

  for ( j = 2; j < nyp4-2; j++ )
    {
      xvals[j*nxp4  ] = intlin(3.0, xvals[j*nxp4+3], 0.0, xvals[j*nxp4+2], 1.0);
      xvals[j*nxp4+1] = intlin(2.0, xvals[j*nxp4+3], 0.0, xvals[j*nxp4+2], 1.0); 
      yvals[j*nxp4  ] = intlin(3.0, yvals[j*nxp4+3], 0.0, yvals[j*nxp4+2], 1.0); 
      yvals[j*nxp4+1] = intlin(2.0, yvals[j*nxp4+3], 0.0, yvals[j*nxp4+2], 1.0); 

      xvals[j*nxp4+nxp4-2] = intlin(2.0, xvals[j*nxp4+nxp4-4], 0.0, xvals[j*nxp4+nxp4-3], 1.0); 
      xvals[j*nxp4+nxp4-1] = intlin(3.0, xvals[j*nxp4+nxp4-4], 0.0, xvals[j*nxp4+nxp4-3], 1.0); 
      yvals[j*nxp4+nxp4-2] = intlin(2.0, yvals[j*nxp4+nxp4-4], 0.0, yvals[j*nxp4+nxp4-3], 1.0); 
      yvals[j*nxp4+nxp4-1] = intlin(3.0, yvals[j*nxp4+nxp4-4], 0.0, yvals[j*nxp4+nxp4-3], 1.0); 
    }

  for ( i = 0; i < nxp4; i++ )
    {
      xvals[0*nxp4+i] = intlin(3.0, xvals[3*nxp4+i], 0.0, xvals[2*nxp4+i], 1.0);
      xvals[1*nxp4+i] = intlin(2.0, xvals[3*nxp4+i], 0.0, xvals[2*nxp4+i], 1.0);
      yvals[0*nxp4+i] = intlin(3.0, yvals[3*nxp4+i], 0.0, yvals[2*nxp4+i], 1.0);
      yvals[1*nxp4+i] = intlin(2.0, yvals[3*nxp4+i], 0.0, yvals[2*nxp4+i], 1.0);

      xvals[(nyp4-2)*nxp4+i] = intlin(2.0, xvals[(nyp4-4)*nxp4+i], 0.0, xvals[(nyp4-3)*nxp4+i], 1.0);
      xvals[(nyp4-1)*nxp4+i] = intlin(3.0, xvals[(nyp4-4)*nxp4+i], 0.0, xvals[(nyp4-3)*nxp4+i], 1.0);
      yvals[(nyp4-2)*nxp4+i] = intlin(2.0, yvals[(nyp4-4)*nxp4+i], 0.0, yvals[(nyp4-3)*nxp4+i], 1.0);
      yvals[(nyp4-1)*nxp4+i] = intlin(3.0, yvals[(nyp4-4)*nxp4+i], 0.0, yvals[(nyp4-3)*nxp4+i], 1.0);
    }
  /*
    {
    FILE *fp;
    fp = fopen("xvals.asc", "w");
    for ( i = 0; i < gridsize_new; i++ ) fprintf(fp, "%g\n", xvals[i]);
    fclose(fp);
    fp = fopen("yvals.asc", "w");
    for ( i = 0; i < gridsize_new; i++ ) fprintf(fp, "%g\n", yvals[i]);
    fclose(fp);
    }
  */
  gridDefXvals(gridIDnew, xvals);
  gridDefYvals(gridIDnew, yvals);
  
  free(xvals);
  free(yvals);

  return(gridIDnew);
}

/*****************************************************************************/

static
void grid_check_lat_borders_rad(int n, double *ybounds)
{
  if ( ybounds[0] > ybounds[n-1] )
    {
      if ( RAD2DEG*ybounds[0]   >  88 ) ybounds[0]   =  PIH;
      if ( RAD2DEG*ybounds[n-1] < -88 ) ybounds[n-1] = -PIH;
    }
  else
    {
      if ( RAD2DEG*ybounds[0]   < -88 ) ybounds[0]   = -PIH;
      if ( RAD2DEG*ybounds[n-1] >  88 ) ybounds[n-1] =  PIH;
    }
}

static
void remap_define_reg2d(int gridID, remapgrid_t *grid)
{
  char unitstr[CDI_MAX_NAME];
  long nx, nxm, ny;
  long nxp1, nyp1;

  nx = grid->dims[0];
  ny = grid->dims[1];

  nxp1 = nx + 1;
  nyp1 = ny + 1;

  nxm = nx;
  if ( grid->is_cyclic ) nxm++;

  if ( grid->size != nx*ny ) cdoAbort("Internal error, wrong dimensions!");

  grid->reg2d_center_lon = realloc(grid->reg2d_center_lon, nxm*sizeof(double));
  grid->reg2d_center_lat = realloc(grid->reg2d_center_lat,  ny*sizeof(double));
 
  gridInqXvals(gridID, grid->reg2d_center_lon);
  gridInqYvals(gridID, grid->reg2d_center_lat);

  /* Convert lat/lon units if required */

  gridInqYunits(gridID, unitstr);

  grid_to_radian(unitstr, nx, grid->reg2d_center_lon, "grid reg2d center lon"); 
  grid_to_radian(unitstr, ny, grid->reg2d_center_lat, "grid reg2d center lat"); 

  if ( grid->is_cyclic ) grid->reg2d_center_lon[nx] = grid->reg2d_center_lon[0] + PI2;

  grid->reg2d_corner_lon = malloc(nxp1*sizeof(double));
  grid->reg2d_corner_lat = malloc(nyp1*sizeof(double));

  grid_gen_corners(nx, grid->reg2d_center_lon, grid->reg2d_corner_lon);
  grid_gen_corners(ny, grid->reg2d_center_lat, grid->reg2d_corner_lat);
  grid_check_lat_borders_rad(ny+1, grid->reg2d_corner_lat);

  //for ( long i = 0; i < nxp1; ++i ) printf("lon %ld %g\n", i, grid->reg2d_corner_lon[i]);
  //for ( long i = 0; i < nyp1; ++i ) printf("lat %ld %g\n", i, grid->reg2d_corner_lat[i]);

}

static
void remap_define_grid(int map_type, int gridID, remapgrid_t *grid)
{
  char xunitstr[CDI_MAX_NAME];
  char yunitstr[CDI_MAX_NAME];
  long gridsize;
  long i;
  int lgrid_destroy = FALSE;
  int lgrid_gen_bounds = FALSE;
  int gridID_gme = -1;

  if ( gridInqType(grid->gridID) != GRID_UNSTRUCTURED && gridInqType(grid->gridID) != GRID_CURVILINEAR )
    {
      if ( gridInqType(grid->gridID) == GRID_GME )
	{
	  gridID_gme = gridToUnstructured(grid->gridID, 1);
	  grid->nvgp = gridInqSize(gridID_gme);
	  gridID = gridDuplicate(gridID_gme);
	  gridCompress(gridID);
	  grid->luse_cell_corners = TRUE;
	}
      else if ( remap_write_remap == TRUE || grid->remap_grid_type != REMAP_GRID_TYPE_REG2D )
	{
	  lgrid_destroy = TRUE;
	  gridID = gridToCurvilinear(grid->gridID, 1);
	  lgrid_gen_bounds = TRUE;
	}
    }

  gridsize = grid->size = gridInqSize(gridID);

  grid->dims[0] = gridInqXsize(gridID);
  grid->dims[1] = gridInqYsize(gridID);

  grid->is_cyclic = gridIsCircular(gridID);

  if ( gridInqType(gridID) == GRID_UNSTRUCTURED )
    grid->rank = 1;
  else
    grid->rank = 2;

  if ( gridInqType(gridID) == GRID_UNSTRUCTURED )
    grid->num_cell_corners = gridInqNvertex(gridID);
  else
    grid->num_cell_corners = 4;

  remapGridRealloc(map_type, grid);

  /* Initialize logical mask */

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(gridsize, grid)
#endif
  for ( i = 0; i < gridsize; ++i ) grid->mask[i] = TRUE;


  if ( remap_write_remap == FALSE && grid->remap_grid_type == REMAP_GRID_TYPE_REG2D ) return;


  gridInqXvals(gridID, grid->cell_center_lon);
  gridInqYvals(gridID, grid->cell_center_lat);

  gridInqXunits(gridID, xunitstr);
  gridInqYunits(gridID, yunitstr);


  if ( grid->lneed_cell_corners )
    {
      if ( gridInqYbounds(gridID, NULL) && gridInqXbounds(gridID, NULL) )
	{
	  gridInqXbounds(gridID, grid->cell_corner_lon);
	  gridInqYbounds(gridID, grid->cell_corner_lat);
	}
      else if ( lgrid_gen_bounds )
	{
	  grid_cell_center_to_bounds_X2D(xunitstr, grid->dims[0], grid->dims[1], grid->cell_center_lon, grid->cell_corner_lon, 0);
	  grid_cell_center_to_bounds_Y2D(yunitstr, grid->dims[0], grid->dims[1], grid->cell_center_lat, grid->cell_corner_lat);
	}
      else
	{
	  cdoAbort("Grid corner missing!");
	}
    }


  if ( gridInqType(grid->gridID) == GRID_GME ) gridInqMaskGME(gridID_gme, grid->vgpm);

  /* Convert lat/lon units if required */

  grid_to_radian(xunitstr, grid->size, grid->cell_center_lon, "grid center lon"); 
  grid_to_radian(yunitstr, grid->size, grid->cell_center_lat, "grid center lat"); 
  /* Note: using units from cell center instead from bounds */
  if ( grid->num_cell_corners && grid->lneed_cell_corners )
    {
      grid_to_radian(xunitstr, grid->num_cell_corners*grid->size, grid->cell_corner_lon, "grid corner lon"); 
      grid_to_radian(yunitstr, grid->num_cell_corners*grid->size, grid->cell_corner_lat, "grid corner lat"); 
    }

  if ( lgrid_destroy ) gridDestroy(gridID);

  /* Convert longitudes to 0,2pi interval */

  check_lon_range(grid->size, grid->cell_center_lon);

  if ( grid->num_cell_corners && grid->lneed_cell_corners )
    {
      //check_lon_range2(grid->num_cell_corners, grid->size, grid->cell_corner_lon, grid->cell_center_lon);
	check_lon_range(grid->num_cell_corners*grid->size, grid->cell_corner_lon);
    }
  /*  Make sure input latitude range is within the machine values for +/- pi/2 */

  check_lat_range(grid->size, grid->cell_center_lat);

  if ( grid->num_cell_corners && grid->lneed_cell_corners )
    check_lat_range(grid->num_cell_corners*grid->size, grid->cell_corner_lat);
}

/*  Compute bounding boxes for restricting future grid searches */
static
void cell_bounding_boxes(remapgrid_t *grid, int remap_grid_basis)
{
  if ( remap_grid_basis == REMAP_GRID_BASIS_SRC || grid->luse_cell_corners )
    grid->cell_bound_box = realloc(grid->cell_bound_box, 4*grid->size*sizeof(restr_t));

  if ( grid->luse_cell_corners )
    {
      if ( grid->lneed_cell_corners )
	{
	  if ( cdoVerbose ) cdoPrint("Grid: boundbox_from_corners");

	  boundbox_from_corners(grid->size, grid->num_cell_corners, 
				grid->cell_corner_lon, grid->cell_corner_lat, grid->cell_bound_box);
	}
      else /* full grid search */
	{
	  long gridsize;
	  long i, i4;
	  
	  gridsize = grid->size;
  
	  if ( cdoVerbose ) cdoPrint("Grid: bounds missing -> full grid search!");

	  for ( i = 0; i < gridsize; ++i )
	    {
	      i4 = i<<2;
	      grid->cell_bound_box[i4  ] = RESTR_SCALE(-PIH);
	      grid->cell_bound_box[i4+1] = RESTR_SCALE( PIH);
	      grid->cell_bound_box[i4+2] = RESTR_SCALE(0.);
	      grid->cell_bound_box[i4+3] = RESTR_SCALE(PI2);
	    }
	}
    }
  else if ( remap_grid_basis == REMAP_GRID_BASIS_SRC )
    {
      long nx, ny;

      if ( grid->rank != 2 ) cdoAbort("Internal problem, grid rank = %d!", grid->rank);

      nx = grid->dims[0];
      ny = grid->dims[1];

      if ( cdoVerbose ) cdoPrint("Grid: boundbox_from_center");

      boundbox_from_center(grid->is_cyclic, grid->size, nx, ny, 
			   grid->cell_center_lon, grid->cell_center_lat, grid->cell_bound_box);
    }

  if ( remap_grid_basis == REMAP_GRID_BASIS_SRC || grid->lneed_cell_corners )
    check_lon_boundbox_range(grid->size, grid->cell_bound_box);

  /* Try to check for cells that overlap poles */

  if ( remap_grid_basis == REMAP_GRID_BASIS_SRC || grid->lneed_cell_corners )
    check_lat_boundbox_range(grid->size, grid->cell_bound_box, grid->cell_center_lat);
}


void remap_grids_init(int map_type, int lextrapolate, int gridID1, remapgrid_t *src_grid, int gridID2, remapgrid_t *tgt_grid)
{
  int reg2d_src_gridID = gridID1;
  int reg2d_tgt_gridID = gridID2;

  /* Initialize remapgrid structure */
  remapgrid_init(src_grid);
  remapgrid_init(tgt_grid);

  if ( map_type == MAP_TYPE_BILINEAR || map_type == MAP_TYPE_BICUBIC ||
       map_type == MAP_TYPE_DISTWGT  || map_type == MAP_TYPE_CONSPHERE )
    {
      if ( IS_REG2D_GRID(gridID1) ) src_grid->remap_grid_type = REMAP_GRID_TYPE_REG2D;
      // src_grid->remap_grid_type = 0;
    }

  if ( map_type == MAP_TYPE_CONSPHERE && src_grid->remap_grid_type == REMAP_GRID_TYPE_REG2D )
    {
      if ( IS_REG2D_GRID(gridID2) ) tgt_grid->remap_grid_type = REMAP_GRID_TYPE_REG2D;
      else src_grid->remap_grid_type = -1;
    }

  if ( lextrapolate > 0 )
    src_grid->lextrapolate = TRUE;
  else
    src_grid->lextrapolate = FALSE;

  if ( map_type == MAP_TYPE_CONSERV || map_type == MAP_TYPE_CONSPHERE )
    {
      if ( src_grid->remap_grid_type != REMAP_GRID_TYPE_REG2D )
	{
	  src_grid->luse_cell_corners  = TRUE;
	  src_grid->lneed_cell_corners = TRUE;
	}

      if ( tgt_grid->remap_grid_type != REMAP_GRID_TYPE_REG2D )
	{
	  tgt_grid->luse_cell_corners  = TRUE;
	  tgt_grid->lneed_cell_corners = TRUE;
	}
    }

  src_grid->gridID = gridID1;
  tgt_grid->gridID = gridID2;

  if ( !src_grid->lextrapolate && gridInqSize(src_grid->gridID) > 1 &&
       map_type == MAP_TYPE_DISTWGT &&
       ((gridInqType(gridID1) == GRID_LONLAT && gridIsRotated(gridID1)) ||
	(gridInqType(gridID1) == GRID_LONLAT && src_grid->non_global)) )
    {
      src_grid->gridID = gridID1 = expand_lonlat_grid(gridID1);
    }

  if ( gridInqType(gridID1) == GRID_UNSTRUCTURED )
    {
      if ( gridInqYvals(gridID1, NULL) == 0 || gridInqXvals(gridID1, NULL) == 0 )
	{
	  if ( gridInqNumber(gridID1) > 0 )
	    {
	      src_grid->gridID = gridID1 = referenceToGrid(gridID1);
	      if ( gridID1 == -1 ) cdoAbort("Reference to source grid not found!");
	    }
	}
    }

  if ( gridInqType(gridID2) == GRID_UNSTRUCTURED )
    {
      if ( gridInqYvals(gridID2, NULL) == 0 || gridInqXvals(gridID2, NULL) == 0 )
	{
	  if ( gridInqNumber(gridID2) > 0 )
	    {
	      tgt_grid->gridID = gridID2 = referenceToGrid(gridID2);
	      if ( gridID2 == -1 ) cdoAbort("Reference to target grid not found!");
	    }
	}
    }

  if ( gridInqSize(src_grid->gridID) > 1 && 
       (gridInqType(src_grid->gridID) == GRID_LCC || 
	gridInqType(src_grid->gridID) == GRID_LAEA || 
	gridInqType(src_grid->gridID) == GRID_SINUSOIDAL) )
    {
      src_grid->gridID = gridID1 = gridToCurvilinear(src_grid->gridID, 1);
    }

  if ( !src_grid->lextrapolate && gridInqSize(src_grid->gridID) > 1 &&
       map_type == MAP_TYPE_DISTWGT &&
       (gridInqType(gridID1) == GRID_CURVILINEAR && src_grid->non_global) )
    {
      src_grid->gridID = gridID1 = expand_curvilinear_grid(gridID1);
    }

  if ( map_type == MAP_TYPE_DISTWGT )
    {
      if ( gridInqType(src_grid->gridID) == GRID_UNSTRUCTURED )
	{
	  src_grid->luse_cell_corners  = TRUE;
	  src_grid->lneed_cell_corners = FALSE; /* full grid search */
	}
      /* not used !
      if ( gridInqType(tgt_grid->gridID) == GRID_UNSTRUCTURED )
	{
	  tgt_grid->luse_cell_corners  = TRUE;
	  tgt_grid->lneed_cell_corners = FALSE;
	}
      */
    }

  //if ( src_grid->remap_grid_type != REMAP_GRID_TYPE_REG2D )
    remap_define_grid(map_type, gridID1, src_grid);

  remap_define_grid(map_type, gridID2, tgt_grid);

  if ( src_grid->remap_grid_type == REMAP_GRID_TYPE_REG2D && tgt_grid->remap_grid_type == REMAP_GRID_TYPE_REG2D )
    {
      remap_define_reg2d(reg2d_src_gridID, src_grid);
      remap_define_reg2d(reg2d_tgt_gridID, tgt_grid);
    }
  else if ( src_grid->remap_grid_type == REMAP_GRID_TYPE_REG2D )
    {
      remap_define_reg2d(reg2d_src_gridID, src_grid);
    }
  else
    {
      cell_bounding_boxes(src_grid, REMAP_GRID_BASIS_SRC);
      cell_bounding_boxes(tgt_grid, REMAP_GRID_BASIS_TGT);
      /*
	Set up and assign address ranges to search bins in order to further restrict later searches
      */
      calc_lat_bins(src_grid, tgt_grid, map_type);
    }

}  /* remapGridInit */

/*****************************************************************************/

/*
    This routine initializes some variables and provides an initial
    allocation of arrays (fairly large so frequent resizing unnecessary).
*/
void remap_vars_init(int map_type, long src_grid_size, long tgt_grid_size, remapvars_t *rv)
{
  /* Initialize all pointer */
  if ( rv->pinit == FALSE )
    {
      rv->pinit = TRUE;

      rv->src_grid_add = NULL;
      rv->tgt_grid_add = NULL;
      rv->wts          = NULL;
    }

  /* Determine the number of weights */

#if defined(_OPENMP)
  if ( ompNumThreads > 1 )
    {
      if      ( map_type == MAP_TYPE_CONSERV   ) rv->sort_add = TRUE;
      else if ( map_type == MAP_TYPE_CONSPHERE ) rv->sort_add = TRUE;
      else if ( map_type == MAP_TYPE_BILINEAR  ) rv->sort_add = TRUE;
      else if ( map_type == MAP_TYPE_BICUBIC   ) rv->sort_add = TRUE;
      else if ( map_type == MAP_TYPE_DISTWGT   ) rv->sort_add = TRUE;
      else cdoAbort("Unknown mapping method!");
    }
  else
#endif
    {
      if      ( map_type == MAP_TYPE_CONSERV   ) rv->sort_add = TRUE;
      else if ( map_type == MAP_TYPE_CONSPHERE ) rv->sort_add = TRUE;
      else if ( map_type == MAP_TYPE_BILINEAR  ) rv->sort_add = FALSE;
      else if ( map_type == MAP_TYPE_BICUBIC   ) rv->sort_add = FALSE;
      else if ( map_type == MAP_TYPE_DISTWGT   ) rv->sort_add = TRUE;
      else cdoAbort("Unknown mapping method!");
    }

  if      ( map_type == MAP_TYPE_CONSERV   ) rv->num_wts = 3;
  else if ( map_type == MAP_TYPE_CONSPHERE ) rv->num_wts = 1;
  else if ( map_type == MAP_TYPE_BILINEAR  ) rv->num_wts = 1;
  else if ( map_type == MAP_TYPE_BICUBIC   ) rv->num_wts = 4;
  else if ( map_type == MAP_TYPE_DISTWGT   ) rv->num_wts = 1;
  else cdoAbort("Unknown mapping method!");

   /*
    Initialize num_links and set max_links to four times the largest 
    of the destination grid sizes initially (can be changed later).
    Set a default resize increment to increase the size of link
    arrays if the number of links exceeds the initial size
  */
  rv->num_links = 0;
  rv->max_links = 4 * tgt_grid_size;

  rv->resize_increment = (int) (0.1 * MAX(src_grid_size, tgt_grid_size));

  /*  Allocate address and weight arrays for mapping 1 */

  rv->src_grid_add = realloc(rv->src_grid_add, rv->max_links*sizeof(int));
  rv->tgt_grid_add = realloc(rv->tgt_grid_add, rv->max_links*sizeof(int));

  rv->wts = realloc(rv->wts, rv->num_wts*rv->max_links*sizeof(double));

  rv->links.option    = FALSE;
  rv->links.max_links = 0;
  rv->links.num_blks  = 0;
  rv->links.num_links = NULL;
  rv->links.src_add   = NULL;
  rv->links.dst_add   = NULL;
  rv->links.w_index   = NULL;

} /* remapVarsInit */

/*****************************************************************************/

/*
   This routine resizes remapping arrays by increasing(decreasing) the max_links by increment
*/
void resize_remap_vars(remapvars_t *rv, int increment)
{
  /*
    Input variables:
    int  increment  ! the number of links to add(subtract) to arrays
  */

  /*  Reallocate arrays at new size */

  rv->max_links += increment;

  if ( rv->max_links )
    {
      rv->src_grid_add = realloc(rv->src_grid_add, rv->max_links*sizeof(int));
      rv->tgt_grid_add = realloc(rv->tgt_grid_add, rv->max_links*sizeof(int));

      rv->wts = realloc(rv->wts, rv->num_wts*rv->max_links*sizeof(double));
    }

} /* resize_remap_vars */


/*
  -----------------------------------------------------------------------

  Performs the remapping based on weights computed elsewhere
     
  -----------------------------------------------------------------------
*/
void remap(double *restrict dst_array, double missval, long dst_size, long num_links, double *restrict map_wts, 
	   long num_wts, const int *restrict dst_add, const int *restrict src_add, const double *restrict src_array, 
	   const double *restrict src_grad1, const double *restrict src_grad2, const double *restrict src_grad3,
	   remaplink_t links)
{
  /*
    Input arrays:

    int *dst_add         ! destination address for each link
    int *src_add         ! source      address for each link

    int num_wts          ! num of weights used in remapping

    double *map_wts      ! remapping weights for each link

    double *src_array    ! array with source field to be remapped

    optional:

    double *src_grad1    ! gradient arrays on source grid necessary for
    double *src_grad2    ! higher-order remappings
    double *src_grad3

    output variables:

    double *dst_array    ! array for remapped field on destination grid
  */

  /* Local variables */
  long n;
  int iorder;

  /* Check the order of the interpolation */

  if ( src_grad1 )
    iorder = 2;
  else
    iorder = 1;

  for ( n = 0; n < dst_size; ++n ) dst_array[n] = missval;

  if ( cdoTimer ) timer_start(timer_remap);

#ifdef SX
#pragma cdir nodep
#endif
  for ( n = 0; n < num_links; ++n ) dst_array[dst_add[n]] = ZERO;

  if ( iorder == 1 )   /* First order remapping */
    {
      if ( links.option == TRUE )
	{
	  long j;
	  for ( j = 0; j < links.num_blks; ++j )
	    {
#ifdef SX
#pragma cdir nodep
#endif
	      for ( n = 0; n < links.num_links[j]; ++n )
		{
		  dst_array[links.dst_add[j][n]] += src_array[links.src_add[j][n]]*map_wts[num_wts*links.w_index[j][n]];
		}
	    }
	}
      else
	{
	  for ( n = 0; n < num_links; ++n )
	    {
	      /*
		printf("%5d %5d %5d %g # dst_add src_add n\n", dst_add[n], src_add[n], n, map_wts[num_wts*n]);
	      */
	      dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[num_wts*n];
	    }
	}
    }
  else                 /* Second order remapping */
    {
      if ( num_wts == 3 )
	{
	  for ( n = 0; n < num_links; ++n )
	    {
	      dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[num_wts*n] +
                                       src_grad1[src_add[n]]*map_wts[num_wts*n+1] +
                                       src_grad2[src_add[n]]*map_wts[num_wts*n+2];
	    }
	}
      else if ( num_wts == 4 )
	{
      	  for ( n = 0; n < num_links; ++n )
	    {
              dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[num_wts*n] +
                                       src_grad1[src_add[n]]*map_wts[num_wts*n+1] +
                                       src_grad2[src_add[n]]*map_wts[num_wts*n+2] +
                                       src_grad3[src_add[n]]*map_wts[num_wts*n+3];
	    }
	}
    }

  if ( cdoTimer ) timer_stop(timer_remap);
}

static
long get_max_add(long num_links, long size, const int *restrict add)
{
  long n, i;
  long max_add;
  int *isum;

  isum = malloc(size*sizeof(int));
  memset(isum, 0, size*sizeof(int));

  for ( n = 0; n < num_links; ++n ) isum[add[n]]++;

  max_add = 0;
  for ( i = 0; i < size; ++i ) if ( isum[i] > max_add ) max_add = isum[i];
  free(isum);

  return (max_add);
}

static 
long binary_search_int(const int *array, long len, int value)
{       
  long low = 0, high = len - 1, midpoint = 0;
 
  while ( low <= high )
    {
      midpoint = low + (high - low)/2;      
 
      // check to see if value is equal to item in array
      if ( value == array[midpoint] ) return midpoint;

      if ( value < array[midpoint] )
	high = midpoint - 1;
      else
	low  = midpoint + 1;
    }
 
  // item was not found
  return -1L;
}

/*
  -----------------------------------------------------------------------

  Performs the remapping based on weights computed elsewhere
     
  -----------------------------------------------------------------------
*/
void remap_laf(double *restrict dst_array, double missval, long dst_size, long num_links, double *restrict map_wts,
	       long num_wts, const int *restrict dst_add, const int *restrict src_add, const double *restrict src_array)
{
  /*
    Input arrays:

    int *dst_add         ! destination address for each link
    int *src_add         ! source      address for each link

    int num_wts          ! num of weights used in remapping

    double *map_wts      ! remapping weights for each link

    double *src_array    ! array with source field to be remapped

    output variables:

    double *dst_array    ! array for remapped field on destination grid
  */

  /* Local variables */
  long i, n, k, ncls, imax;
  long max_cls;
  double wts;
  double *src_cls;
  double *src_wts;
#if defined(_OPENMP)
  double **src_cls2;
  double **src_wts2;
  int ompthID;
#endif

  for ( i = 0; i < dst_size; ++i ) dst_array[i] = missval;

  if ( num_links == 0 ) return;

  max_cls = get_max_add(num_links, dst_size, dst_add);

#if defined(_OPENMP)
  src_cls2 = malloc(ompNumThreads*sizeof(double *));
  src_wts2 = malloc(ompNumThreads*sizeof(double *));
  for ( i = 0; i < ompNumThreads; ++i )
    {
      src_cls2[i] = malloc(max_cls*sizeof(double));
      src_wts2[i] = malloc(max_cls*sizeof(double));
    }
#else
  src_cls = malloc(max_cls*sizeof(double));
  src_wts = malloc(max_cls*sizeof(double));
#endif

  for ( n = 0; n < num_links; ++n )
    if ( DBL_IS_EQUAL(dst_array[dst_add[n]], missval) ) dst_array[dst_add[n]] = ZERO;

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(dst_size, src_cls2, src_wts2, num_links, dst_add, src_add, src_array, map_wts, \
	 num_wts, dst_array, max_cls)					\
  private(i, n, k, ompthID, src_cls, src_wts, ncls, imax, wts) \
  schedule(dynamic,1)
#endif
  for ( i = 0; i < dst_size; ++i )
    {
#if defined(_OPENMP)
      ompthID = omp_get_thread_num();
      src_cls = src_cls2[ompthID];
      src_wts = src_wts2[ompthID];
#endif
      memset(src_cls, 0, max_cls*sizeof(double));
      memset(src_wts, 0, max_cls*sizeof(double));
      /*
      ncls = 0;
      for ( n = 0; n < num_links; n++ )
	{
	  if ( i == dst_add[n] )
	    {
	      for ( k = 0; k < ncls; k++ )
		if ( IS_EQUAL(src_array[src_add[n]], src_cls[k]) ) break;
	      
	      if ( k == ncls )
		{
		  src_cls[k] = src_array[src_add[n]];
		  ncls++;
		}
	      
	      src_wts[k] += map_wts[num_wts*n];
	    }
	}
      */
      /* only for sorted dst_add! */
      {
      long min_add = 1, max_add = 0;

      n = binary_search_int(dst_add, num_links, (int)i);

      if ( n >= 0 && n < num_links )
	{
	  min_add = n;
	  
	  for ( n = min_add+1; n < num_links; ++n )
	    if ( i != dst_add[n] ) break;

	  max_add = n;

	  for ( n = min_add; n > 0; --n )
	    if ( i != dst_add[n-1] ) break;

	  min_add = n;
	}

      ncls = 0;
      for ( n = min_add; n < max_add; ++n )
	{
	  for ( k = 0; k < ncls; ++k )
	    if ( IS_EQUAL(src_array[src_add[n]], src_cls[k]) ) break;
	      
	  if ( k == ncls )
	    {
	      src_cls[k] = src_array[src_add[n]];
	      ncls++;
	    }
	      
	  src_wts[k] += map_wts[num_wts*n];
	}
      }
      
      if ( ncls )
	{
	  imax = 0;
	  wts = src_wts[0];
	  for ( k = 1; k < ncls; ++k )
	    {
	      if ( src_wts[k] > wts )
		{
		  wts  = src_wts[k];
		  imax = k;
		}
	    }

	  dst_array[i] = src_cls[imax];
	}
    }

#if defined(_OPENMP)
  for ( i = 0; i < ompNumThreads; ++i )
    {
      free(src_cls2[i]);
      free(src_wts2[i]);
    }

  free(src_cls2);
  free(src_wts2);
#else
  free(src_cls);
  free(src_wts);
#endif
}

/*
  -----------------------------------------------------------------------

  Performs the remapping based on weights computed elsewhere
     
  -----------------------------------------------------------------------
*/
void remap_sum(double *restrict dst_array, double missval, long dst_size, long num_links, double *restrict map_wts,
	       long num_wts, const int *restrict dst_add, const int *restrict src_add, const double *restrict src_array)
{
  /*
    Input arrays:

    int *dst_add         ! destination address for each link
    int *src_add         ! source      address for each link

    int num_wts          ! num of weights used in remapping

    double *map_wts      ! remapping weights for each link

    double *src_array    ! array with source field to be remapped

    output variables:

    double *dst_array    ! array for remapped field on destination grid
  */
  /* Local variables */
  long n;

  for ( n = 0; n < dst_size; ++n ) dst_array[n] = missval;

#ifdef SX
#pragma cdir nodep
#endif
  for ( n = 0; n < num_links; ++n )
    if ( DBL_IS_EQUAL(dst_array[dst_add[n]], missval) ) dst_array[dst_add[n]] = ZERO;

  for ( n = 0; n < num_links; ++n )
    {
      /*
	printf("%5d %5d %5d %g # dst_add src_add n\n", dst_add[n], src_add[n], n, map_wts[num_wts*n]);
      */
      //dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[num_wts*n];
      dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[num_wts*n];
      printf("%ld %d %d %g %g %g\n", n, dst_add[n], src_add[n],
	     src_array[src_add[n]], map_wts[num_wts*n], dst_array[dst_add[n]]);
    }
}


static double  converge = 1.e-10;            /* Convergence criterion */

/* threshold for coord transformation */
void remap_set_threshhold(double threshhold)
{
  north_thresh =  threshhold;
  south_thresh = -threshhold;  

  if ( cdoVerbose ) cdoPrint("threshhold: north=%g  south=%g", north_thresh, south_thresh);
}


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      BILINEAR INTERPOLATION                                             */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

long find_element(double x, long nelem, const double *restrict array);
int rect_grid_search(long *ii, long *jj, double x, double y, long nxm, long nym, const double *restrict xm, const double *restrict ym);

static
int grid_search_reg2d_nn(long nx, long ny, int *restrict nbr_add, double *restrict nbr_dist, double plat, double plon,
			 const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  int search_result = 0;
  long n, srch_add;
  long i;
  long ii, jj;
  long jjskip;
  double coslat, sinlat;
  double dist_min, distance; /* For computing dist-weighted avg */
  double *sincoslon;
  double coslat_dst = cos(plat);
  double sinlat_dst = sin(plat);
  double coslon_dst = cos(plon);
  double sinlon_dst = sin(plon);

  dist_min = BIGNUM;
  for ( n = 0; n < 4; ++n ) nbr_dist[n] = BIGNUM;  

  long jjf = 0, jjl = ny-1;
  if ( plon >= src_center_lon[0] && plon <= src_center_lon[nx-1] )
    {
      if ( src_center_lat[0] < src_center_lat[ny-1] )
	{
	  if ( plat <= src_center_lat[0] )
	    { jjf = 0; jjl = 1; }
	  else
	    { jjf = ny-2; jjl = ny-1; }
	}
      else
	{
	  if ( plat >= src_center_lat[0] )
	    { jjf = 0; jjl = 1; }
	  else
	    { jjf = ny-2; jjl = ny-1; }
	}
    }

  sincoslon = malloc(nx*sizeof(double));

  for ( ii = 0; ii < nx; ++ii )
    sincoslon[ii] = coslon_dst*cos(src_center_lon[ii]) + sinlon_dst*sin(src_center_lon[ii]);

  for ( jj = jjf; jj <= jjl; ++jj )
    {
      coslat = coslat_dst*cos(src_center_lat[jj]);
      sinlat = sinlat_dst*sin(src_center_lat[jj]);

      jjskip = jj > 1 && jj < (ny-2);

      for ( ii = 0; ii < nx; ++ii )
	{
	  if ( jjskip && ii > 1 && ii < (nx-2) ) continue;

	  srch_add = jj*nx + ii;

	  distance = acos(coslat*sincoslon[ii] + sinlat);

	  if ( distance < dist_min )
	    {
	      for ( n = 0; n < 4; ++n )
		{
		  if ( distance < nbr_dist[n] )
		    {
		      for ( i = 3; i > n; --i )
			{
			  nbr_add [i] = nbr_add [i-1];
			  nbr_dist[i] = nbr_dist[i-1];
			}
		      search_result = -1;
		      nbr_add [n] = srch_add;
		      nbr_dist[n] = distance;
		      dist_min = nbr_dist[3];
		      break;
		    }
		}
	    }
	}
    }

  free(sincoslon);

  for ( n = 0; n < 4; ++n ) nbr_dist[n] = ONE/(nbr_dist[n] + TINY);
  distance = 0.0;
  for ( n = 0; n < 4; ++n ) distance += nbr_dist[n];
  for ( n = 0; n < 4; ++n ) nbr_dist[n] /= distance;

  return (search_result);
}

static
int grid_search_reg2d(remapgrid_t *src_grid, int *restrict src_add, double *restrict src_lats, 
		      double *restrict src_lons,  double plat, double plon, const int *restrict src_grid_dims,
		      const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  /*
    Output variables:

    int    src_add[4]              ! address of each corner point enclosing P
    double src_lats[4]             ! latitudes  of the four corner points
    double src_lons[4]             ! longitudes of the four corner points

    Input variables:

    double plat                    ! latitude  of the search point
    double plon                    ! longitude of the search point

    int src_grid_dims[2]           ! size of each src grid dimension

    double src_center_lat[]        ! latitude  of each src grid center 
    double src_center_lon[]        ! longitude of each src grid center
  */
  /*  Local variables */
  int search_result = 0;
  int lfound;
  long n;
  long nx, nxm, ny;
  long ii, iix, jj;

  for ( n = 0; n < 4; ++n ) src_add[n] = 0;

  nx = src_grid_dims[0];
  ny = src_grid_dims[1];

  nxm = nx;
  if ( src_grid->is_cyclic ) nxm++;

  if ( /*plon < 0   &&*/ plon < src_center_lon[0]     ) plon += PI2;
  if ( /*plon > PI2 &&*/ plon > src_center_lon[nxm-1] ) plon -= PI2;

  lfound = rect_grid_search(&ii, &jj, plon, plat, nxm, ny, src_center_lon, src_center_lat);

  if ( lfound )
    {
      iix = ii;
      if ( src_grid->is_cyclic && iix == (nxm-1) ) iix = 0;
      src_add[0] = (jj-1)*nx+(ii-1);
      src_add[1] = (jj-1)*nx+(iix);
      src_add[2] = (jj)*nx+(iix);
      src_add[3] = (jj)*nx+(ii-1);

      src_lons[0] = src_center_lon[ii-1];
      src_lons[1] = src_center_lon[iix];
      /* For consistency, we must make sure all lons are in same 2pi interval */
      if ( src_lons[0] > PI2 ) src_lons[0] -= PI2;
      if ( src_lons[0] < 0   ) src_lons[0] += PI2;
      if ( src_lons[1] > PI2 ) src_lons[1] -= PI2;
      if ( src_lons[1] < 0   ) src_lons[1] += PI2;
      src_lons[2] = src_lons[1];
      src_lons[3] = src_lons[0];

      src_lats[0] = src_center_lat[jj-1];
      src_lats[1] = src_lats[0];
      src_lats[2] = src_center_lat[jj];
      src_lats[3] = src_lats[2];

      search_result = 1;
      
      return (search_result);
    }

  /*
    If no cell found, point is likely either in a box that straddles either pole or is outside 
    the grid. Fall back to a distance-weighted average of the four closest points.
    Go ahead and compute weights here, but store in src_lats and return -add to prevent the 
    parent routine from computing bilinear weights.
  */
  if ( !src_grid->lextrapolate ) return (search_result);

  /*
    printf("Could not find location for %g %g\n", plat*RAD2DEG, plon*RAD2DEG);
    printf("Using nearest-neighbor average for this point\n");
  */
  search_result = grid_search_reg2d_nn(nx, ny, src_add, src_lats, plat, plon, src_center_lat, src_center_lon);

  return (search_result);
}  /* grid_search_reg2d */


static
int grid_search_nn(long min_add, long max_add, int *restrict nbr_add, double *restrict nbr_dist, 
		   double plat, double plon,
		   const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  int search_result = 0;
  long n, srch_add;
  long i;
  double dist_min, distance; /* For computing dist-weighted avg */
  double coslat_dst = cos(plat);
  double sinlat_dst = sin(plat);
  double coslon_dst = cos(plon);
  double sinlon_dst = sin(plon);

  dist_min = BIGNUM;
  for ( n = 0; n < 4; ++n ) nbr_dist[n] = BIGNUM;
  for ( srch_add = min_add; srch_add <= max_add; ++srch_add )
    {
      distance = acos(coslat_dst*cos(src_center_lat[srch_add])*
		     (coslon_dst*cos(src_center_lon[srch_add]) +
                      sinlon_dst*sin(src_center_lon[srch_add]))+
		      sinlat_dst*sin(src_center_lat[srch_add]));

      if ( distance < dist_min )
	{
          for ( n = 0; n < 4; ++n )
	    {
	      if ( distance < nbr_dist[n] )
		{
		  for ( i = 3; i > n; --i )
		    {
		      nbr_add [i] = nbr_add [i-1];
		      nbr_dist[i] = nbr_dist[i-1];
		    }
		  search_result = -1;
		  nbr_add [n] = srch_add;
		  nbr_dist[n] = distance;
		  dist_min = nbr_dist[3];
		  break;
		}
	    }
        }
    }

  for ( n = 0; n < 4; ++n ) nbr_dist[n] = ONE/(nbr_dist[n] + TINY);
  distance = 0.0;
  for ( n = 0; n < 4; ++n ) distance += nbr_dist[n];
  for ( n = 0; n < 4; ++n ) nbr_dist[n] /= distance;

  return (search_result);
}

static
int grid_search(remapgrid_t *src_grid, int *restrict src_add, double *restrict src_lats, 
		double *restrict src_lons,  double plat, double plon, const int *restrict src_grid_dims,
		const double *restrict src_center_lat, const double *restrict src_center_lon,
		const restr_t *restrict src_grid_bound_box, const int *restrict src_bin_add)
{
  /*
    Output variables:

    int    src_add[4]              ! address of each corner point enclosing P
    double src_lats[4]             ! latitudes  of the four corner points
    double src_lons[4]             ! longitudes of the four corner points

    Input variables:

    double plat                    ! latitude  of the search point
    double plon                    ! longitude of the search point

    int src_grid_dims[2]           ! size of each src grid dimension

    double src_center_lat[]        ! latitude  of each src grid center 
    double src_center_lon[]        ! longitude of each src grid center

    restr_t src_grid_bound_box[][4] ! bound box for source grid

    int src_bin_add[][2]           ! latitude bins for restricting
  */
  /*  Local variables */
  long n, n2, next_n, srch_add, srch_add4;    /* dummy indices                    */
  long nx, ny;                                /* dimensions of src grid           */
  long min_add, max_add;                      /* addresses for restricting search */
  long i, j, jp1, ip1, n_add, e_add, ne_add;  /* addresses                        */
  long nbins;
  /* Vectors for cross-product check */
  double vec1_lat, vec1_lon;
  double vec2_lat, vec2_lon, cross_product;
  int scross[4], scross_last = 0;
  int search_result = 0;
  restr_t rlat, rlon;
  restr_t *bin_lats = src_grid->bin_lats;

  nbins = src_grid->num_srch_bins;

  rlat = RESTR_SCALE(plat);
  rlon = RESTR_SCALE(plon);

  /* restrict search first using bins */

  for ( n = 0; n < 4; ++n ) src_add[n] = 0;

  min_add = src_grid->size-1;
  max_add = 0;

  for ( n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      if ( rlat >= bin_lats[n2] && rlat <= bin_lats[n2+1] )
	{
	  if ( src_bin_add[n2  ] < min_add ) min_add = src_bin_add[n2  ];
	  if ( src_bin_add[n2+1] > max_add ) max_add = src_bin_add[n2+1];
	}
    }
 
  /* Now perform a more detailed search */

  nx = src_grid_dims[0];
  ny = src_grid_dims[1];

  /* srch_loop */
  for ( srch_add = min_add; srch_add <= max_add; ++srch_add )
    {
      srch_add4 = srch_add<<2;
      /* First check bounding box */
      if ( rlon >= src_grid_bound_box[srch_add4+2] &&
	   rlon <= src_grid_bound_box[srch_add4+3] &&
	   rlat >= src_grid_bound_box[srch_add4  ] &&
	   rlat <= src_grid_bound_box[srch_add4+1])
	{
	  /* We are within bounding box so get really serious */

          /* Determine neighbor addresses */
          j = srch_add/nx;
          i = srch_add - j*nx;

          if ( i < (nx-1) )
            ip1 = i + 1;
          else
	    {
	      /* 2009-01-09 Uwe Schulzweida: bug fix */
	      if ( src_grid->is_cyclic )
		ip1 = 0;
	      else
		ip1 = i;
	    }

          if ( j < (ny-1) )
            jp1 = j + 1;
          else
	    {
	      /* 2008-12-17 Uwe Schulzweida: latitute cyclic ??? (bug fix) */
	      jp1 = j;
	    }

          n_add  = jp1*nx + i;
          e_add  = j  *nx + ip1;
	  ne_add = jp1*nx + ip1;

          src_lons[0] = src_center_lon[srch_add];
          src_lons[1] = src_center_lon[e_add];
          src_lons[2] = src_center_lon[ne_add];
          src_lons[3] = src_center_lon[n_add];

          src_lats[0] = src_center_lat[srch_add];
          src_lats[1] = src_center_lat[e_add];
          src_lats[2] = src_center_lat[ne_add];
          src_lats[3] = src_center_lat[n_add];

	  /* For consistency, we must make sure all lons are in same 2pi interval */

          vec1_lon = src_lons[0] - plon;
          if      ( vec1_lon >  PI ) src_lons[0] -= PI2;
          else if ( vec1_lon < -PI ) src_lons[0] += PI2;

          for ( n = 1; n < 4; ++n )
	    {
	      vec1_lon = src_lons[n] - src_lons[0];
	      if      ( vec1_lon >  PI ) src_lons[n] -= PI2;
	      else if ( vec1_lon < -PI ) src_lons[n] += PI2;
	    }

          /* corner_loop */
          for ( n = 0; n < 4; ++n )
	    {
	      next_n = (n+1)%4;

	      /*
		Here we take the cross product of the vector making 
		up each box side with the vector formed by the vertex
		and search point.  If all the cross products are 
		positive, the point is contained in the box.
	      */
	      vec1_lat = src_lats[next_n] - src_lats[n];
	      vec1_lon = src_lons[next_n] - src_lons[n];
	      vec2_lat = plat - src_lats[n];
	      vec2_lon = plon - src_lons[n];

	      /* Check for 0,2pi crossings */

	      if      ( vec1_lon >  THREE*PIH ) vec1_lon -= PI2;
	      else if ( vec1_lon < -THREE*PIH ) vec1_lon += PI2;

	      if      ( vec2_lon >  THREE*PIH ) vec2_lon -= PI2;
	      else if ( vec2_lon < -THREE*PIH ) vec2_lon += PI2;

	      cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat;

	      /* If cross product is less than ZERO, this cell doesn't work    */
	      /* 2008-10-16 Uwe Schulzweida: bug fix for cross_product eq zero */

	      scross[n] = cross_product < 0 ? -1 : cross_product > 0 ? 1 : 0;

	      if ( n == 0 ) scross_last = scross[n];

	      if ( (scross[n] < 0 && scross_last > 0) || (scross[n] > 0 && scross_last < 0) ) break;

	      scross_last = scross[n];
	    } /* corner_loop */

	  if ( n >= 4 )
	    {
	      n = 0;
	      if      ( scross[0]>=0 && scross[1]>=0 && scross[2]>=0 && scross[3]>=0 ) n = 4;
	      else if ( scross[0]<=0 && scross[1]<=0 && scross[2]<=0 && scross[3]<=0 ) n = 4;
	    }

	  /* If cross products all same sign, we found the location */
          if ( n >= 4 )
	    {
	      src_add[0] = srch_add;
	      src_add[1] = e_add;
	      src_add[2] = ne_add;
	      src_add[3] = n_add;

	      search_result = 1;

	      return (search_result);
	    }

	  /* Otherwise move on to next cell */

        } /* Bounding box check */
    } /* srch_loop */

  /*
    If no cell found, point is likely either in a box that straddles either pole or is outside 
    the grid. Fall back to a distance-weighted average of the four closest points.
    Go ahead and compute weights here, but store in src_lats and return -add to prevent the 
    parent routine from computing bilinear weights.
  */
  if ( !src_grid->lextrapolate ) return (search_result);

  /*
    printf("Could not find location for %g %g\n", plat*RAD2DEG, plon*RAD2DEG);
    printf("Using nearest-neighbor average for this point\n");
  */
  search_result = grid_search_nn(min_add, max_add, src_add, src_lats, plat, plon, src_center_lat, src_center_lon);

  return (search_result);
}  /* grid_search */

/*
  This routine stores the address and weight for four links associated with one destination
  point in the appropriate address and weight arrays and resizes those arrays if necessary.
*/
void store_link_bilin(remapvars_t *rv, int dst_add, int *restrict src_add, double *restrict weights)
{
  /*
    Input variables:
    int dst_add       ! address on destination grid
    int src_add[4]    ! addresses on source grid
    double weights[4] ! array of remapping weights for these links
  */
  long n;
  long nlink;

  /*
     Increment number of links and check to see if remap arrays need
     to be increased to accomodate the new link. Then store the link.
  */
  nlink = rv->num_links;
  rv->num_links += 4;

  if ( rv->num_links >= rv->max_links ) 
    resize_remap_vars(rv, rv->resize_increment);

  if ( rv->sort_add == FALSE )
    {
      for ( n = 0; n < 3; ++n )
	if ( src_add[n] > src_add[n+1] ) break;

      if ( n < 2 ) rv->sort_add = TRUE;
      else if ( n == 2 ) // swap 3rd and 4th elements
	{
	  {
	    int itmp;
	    itmp = src_add[2];
	    src_add[2] = src_add[3];
	    src_add[3] = itmp;
	  }
	  {
	    double dtmp;
	    dtmp = weights[2];
	    weights[2] = weights[3];
	    weights[3] = dtmp;
	  }
	}
    }

  for ( n = 0; n < 4; ++n )
    {
      rv->src_grid_add[nlink+n] = src_add[n];
      rv->tgt_grid_add[nlink+n] = dst_add;
      rv->wts         [nlink+n] = weights[n];
    }

} /* store_link_bilin */

static
long find_ij_weights(double plon, double plat, double *restrict src_lats, double *restrict src_lons, double *ig, double *jg)
{
  long iter;                     /*  iteration counters   */
  double iguess, jguess;         /*  current guess for bilinear coordinate  */
  double deli, delj;             /*  corrections to iw,jw                   */
  double dth1, dth2, dth3;       /*  some latitude  differences             */
  double dph1, dph2, dph3;       /*  some longitude differences             */
  double dthp, dphp;             /*  difference between point and sw corner */
  double mat1, mat2, mat3, mat4; /*  matrix elements                        */
  double determinant;            /*  matrix determinant                     */

  /* Iterate to find iw,jw for bilinear approximation  */

  dth1 = src_lats[1] - src_lats[0];
  dth2 = src_lats[3] - src_lats[0];
  dth3 = src_lats[2] - src_lats[1] - dth2;

  dph1 = src_lons[1] - src_lons[0];
  dph2 = src_lons[3] - src_lons[0];
  dph3 = src_lons[2] - src_lons[1];

  if ( dph1 >  THREE*PIH ) dph1 -= PI2;
  if ( dph2 >  THREE*PIH ) dph2 -= PI2;
  if ( dph3 >  THREE*PIH ) dph3 -= PI2;
  if ( dph1 < -THREE*PIH ) dph1 += PI2;
  if ( dph2 < -THREE*PIH ) dph2 += PI2;
  if ( dph3 < -THREE*PIH ) dph3 += PI2;

  dph3 = dph3 - dph2;

  iguess = HALF;
  jguess = HALF;

  for ( iter = 0; iter < remap_max_iter; ++iter )
    {
      dthp = plat - src_lats[0] - dth1*iguess - dth2*jguess - dth3*iguess*jguess;
      dphp = plon - src_lons[0];
      
      if ( dphp >  THREE*PIH ) dphp -= PI2;
      if ( dphp < -THREE*PIH ) dphp += PI2;

      dphp = dphp - dph1*iguess - dph2*jguess - dph3*iguess*jguess;

      mat1 = dth1 + dth3*jguess;
      mat2 = dth2 + dth3*iguess;
      mat3 = dph1 + dph3*jguess;
      mat4 = dph2 + dph3*iguess;

      determinant = mat1*mat4 - mat2*mat3;

      deli = (dthp*mat4 - dphp*mat2)/determinant;
      delj = (dphp*mat1 - dthp*mat3)/determinant;

      if ( fabs(deli) < converge && fabs(delj) < converge ) break;

      iguess += deli;
      jguess += delj;
    }

  *ig = iguess;
  *jg = jguess;

  return (iter);
}

/*
  -----------------------------------------------------------------------

  This routine computes the weights for a bilinear interpolation.

  -----------------------------------------------------------------------
*/
void remap_bilin(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /*   Local variables */
  int  lwarn = TRUE;
  int  search_result;
  long tgt_grid_size;
  long dst_add;                  /*  destination addresss */
  long n, icount;
  long iter;                     /*  iteration counters   */

  int src_add[4];                /*  address for the four source points     */

  double src_lats[4];            /*  latitudes  of four bilinear corners    */
  double src_lons[4];            /*  longitudes of four bilinear corners    */
  double wgts[4];                /*  bilinear weights for four corners      */

  double plat, plon;             /*  lat/lon coords of destination point    */
  double findex = 0;
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoTimer ) timer_start(timer_remap_bil);

  progressInit();

  tgt_grid_size = tgt_grid->size;

  /* Compute mappings from source to target grid */

  if ( src_grid->rank != 2 )
    cdoAbort("Can not do bilinear interpolation when source grid rank != 2"); 

  /* Loop over destination grid */

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoTimer, cdoVerbose, remap_grid_type, tgt_grid_size, src_grid, tgt_grid, rv, remap_max_iter, converge, lwarn, findex) \
  private(dst_add, n, icount, iter, src_add, src_lats, src_lons, wgts, plat, plon, search_result)    \
  schedule(dynamic,1)
#endif
  /* grid_loop1 */
  for ( dst_add = 0; dst_add < tgt_grid_size; ++dst_add )
    {
      int lprogress = 1;
#if defined(_OPENMP)
      if ( omp_get_thread_num() != 0 ) lprogress = 0;
#endif
#if defined(_OPENMP)
#pragma omp atomic
#endif
      findex++;
      if ( lprogress ) progressStatus(0, 1, findex/tgt_grid_size);

      if ( ! tgt_grid->mask[dst_add] ) continue;

      plat = tgt_grid->cell_center_lat[dst_add];
      plon = tgt_grid->cell_center_lon[dst_add];

      /* Find nearest square of grid points on source grid  */
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	search_result = grid_search_reg2d(src_grid, src_add, src_lats, src_lons, 
					  plat, plon, src_grid->dims,
					  src_grid->reg2d_center_lat, src_grid->reg2d_center_lon);
      else
	search_result = grid_search(src_grid, src_add, src_lats, src_lons, 
				    plat, plon, src_grid->dims,
				    src_grid->cell_center_lat, src_grid->cell_center_lon,
				    src_grid->cell_bound_box, src_grid->bin_addr);

      /* Check to see if points are land points */
      if ( search_result > 0 )
	{
	  for ( n = 0; n < 4; ++n )
	    if ( ! src_grid->mask[src_add[n]] ) search_result = 0;
	}

      /* If point found, find local iw,jw coordinates for weights  */
      if ( search_result > 0 )
	{
	  double iw, jw;  /*  current guess for bilinear coordinate  */

          tgt_grid->cell_frac[dst_add] = ONE;

	  iter = find_ij_weights(plon, plat, src_lats, src_lons, &iw, &jw);

          if ( iter < remap_max_iter )
	    {
	      /* Successfully found iw,jw - compute weights */

	      wgts[0] = (ONE-iw) * (ONE-jw);
	      wgts[1] =      iw  * (ONE-jw);
	      wgts[2] =      iw  *      jw;
	      wgts[3] = (ONE-iw) *      jw;

#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_bilin(rv, dst_add, src_add, wgts);
	    }
          else
	    {
	      if ( cdoVerbose )
		{
		  cdoPrint("Point coords: %g %g", plat, plon);
		  cdoPrint("Src grid lats: %g %g %g %g", src_lats[0], src_lats[1], src_lats[2], src_lats[3]);
		  cdoPrint("Src grid lons: %g %g %g %g", src_lons[0], src_lons[1], src_lons[2], src_lons[3]);
		  cdoPrint("Src grid addresses: %d %d %d %d", src_add[0], src_add[1], src_add[2], src_add[3]);
		  cdoPrint("Src grid lats: %g %g %g %g",
			   src_grid->cell_center_lat[src_add[0]], src_grid->cell_center_lat[src_add[1]],
			   src_grid->cell_center_lat[src_add[2]], src_grid->cell_center_lat[src_add[3]]);
		  cdoPrint("Src grid lons: %g %g %g %g",
			   src_grid->cell_center_lon[src_add[0]], src_grid->cell_center_lon[src_add[1]],
			   src_grid->cell_center_lon[src_add[2]], src_grid->cell_center_lon[src_add[3]]);
		  cdoPrint("Current iw,jw : %g %g", iw, jw);
		}

	      if ( cdoVerbose || lwarn )
		{
		  lwarn = FALSE;
		  //  cdoWarning("Iteration for iw,jw exceed max iteration count of %d!", remap_max_iter);
		  cdoWarning("Bilinear interpolation failed for some grid points - used a distance-weighted average instead!");
		}

	      search_result = -1;
	    }
	}

      /*
	Search for bilinear failed - use a distance-weighted average instead (this is typically near the pole)
	Distance was stored in src_lats!
      */
      if ( search_result < 0 )
	{
          icount = 0;
          for ( n = 0; n < 4; ++n )
	    {
	      if ( src_grid->mask[src_add[n]] )
		icount++;
	      else
		src_lats[n] = ZERO;
	    }

          if ( icount > 0 )
	    {
	      /* Renormalize weights */
	      double sum_wgts = 0.0; /* sum of weights for normalization */
	      /* 2012-05-08 Uwe Schulzweida: using absolute value of src_lats (bug fix) */
	      for ( n = 0; n < 4; ++n ) sum_wgts += fabs(src_lats[n]);
	      for ( n = 0; n < 4; ++n ) wgts[n] = fabs(src_lats[n])/sum_wgts;

	      tgt_grid->cell_frac[dst_add] = ONE;

#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_bilin(rv, dst_add, src_add, wgts);
	    }
        }
    } /* grid_loop1 */

  if ( cdoTimer ) timer_stop(timer_remap_bil);
} /* remap_bilin */


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      BICUBIC INTERPOLATION                                              */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/*
  This routine stores the address and weight for four links associated with one destination 
  point in the appropriate address and weight arrays and resizes those arrays if necessary.
*/
static
void store_link_bicub(remapvars_t *rv, int dst_add, int *restrict src_add, double weights[4][4])
{
  /*
    Input variables:
    int dst_add          ! address on destination grid
    int src_add[4]       ! addresses on source grid
    double weights[4][4] ! array of remapping weights for these links
  */
  long n, k;
  long nlink;

  /*
     Increment number of links and check to see if remap arrays need
     to be increased to accomodate the new link. Then store the link.
  */
  nlink = rv->num_links;
  rv->num_links += 4;

  if ( rv->num_links >= rv->max_links ) 
    resize_remap_vars(rv, rv->resize_increment);

  if ( rv->sort_add == FALSE )
    {
      for ( n = 0; n < 3; ++n )
	if ( src_add[n] > src_add[n+1] ) break;

      if ( n < 2 ) rv->sort_add = TRUE;
      else if ( n == 2 ) // swap 3rd and 4th elements
	{
	  {
	    int itmp;
	    itmp = src_add[2];
	    src_add[2] = src_add[3];
	    src_add[3] = itmp;
	  }
	  {
	    double dtmp[4];
	    for ( k = 0; k < 4; ++k ) dtmp[k] = weights[k][2];
	    for ( k = 0; k < 4; ++k ) weights[k][2] = weights[k][3];
	    for ( k = 0; k < 4; ++k ) weights[k][3] = dtmp[k];
	  }
	}
    }

  for ( n = 0; n < 4; ++n )
    {
      rv->src_grid_add[nlink+n] = src_add[n];
      rv->tgt_grid_add[nlink+n] = dst_add;
      for ( k = 0; k < 4; ++k )
	rv->wts[4*(nlink+n)+k] = weights[k][n];
    }

} /* store_link_bicub */

/*
  -----------------------------------------------------------------------

  This routine computes the weights for a bicubic interpolation.

  -----------------------------------------------------------------------
*/
void remap_bicub(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /*   Local variables */
  int  lwarn = TRUE;
  int  search_result;
  long tgt_grid_size;
  long n, icount;
  long dst_add;        /*  destination addresss */
  long iter;           /*  iteration counters   */

  int src_add[4];     /* address for the four source points */

  double src_lats[4]; /*  latitudes  of four bilinear corners */
  double src_lons[4]; /*  longitudes of four bilinear corners */
  double wgts[4][4];  /*  bicubic weights for four corners    */

  double plat, plon;             /*  lat/lon coords of destination point    */
  double findex = 0;
  int remap_grid_type = src_grid->remap_grid_type;

  progressInit();

  tgt_grid_size = tgt_grid->size;

  /* Compute mappings from source to target grid */

  if ( src_grid->rank != 2 )
    cdoAbort("Can not do bicubic interpolation when source grid rank != 2"); 

  /* Loop over destination grid */

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoTimer, cdoVerbose, remap_grid_type, tgt_grid_size, src_grid, tgt_grid, rv, remap_max_iter, converge, lwarn, findex) \
  private(dst_add, n, icount, iter, src_add, src_lats, src_lons, wgts, plat, plon, search_result) \
  schedule(dynamic,1)
#endif
  /* grid_loop1 */
  for ( dst_add = 0; dst_add < tgt_grid_size; ++dst_add )
    {
      int lprogress = 1;
#if defined(_OPENMP)
      if ( omp_get_thread_num() != 0 ) lprogress = 0;
#endif
#if defined(_OPENMP)
#pragma omp atomic
#endif
      findex++;
      if ( lprogress ) progressStatus(0, 1, findex/tgt_grid_size);

      if ( ! tgt_grid->mask[dst_add] ) continue;

      plat = tgt_grid->cell_center_lat[dst_add];
      plon = tgt_grid->cell_center_lon[dst_add];

      /* Find nearest square of grid points on source grid  */
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	search_result = grid_search_reg2d(src_grid, src_add, src_lats, src_lons, 
					  plat, plon, src_grid->dims,
					  src_grid->reg2d_center_lat, src_grid->reg2d_center_lon);
      else
	search_result = grid_search(src_grid, src_add, src_lats, src_lons, 
				    plat, plon, src_grid->dims,
				    src_grid->cell_center_lat, src_grid->cell_center_lon,
				    src_grid->cell_bound_box, src_grid->bin_addr);

      /* Check to see if points are land points */
      if ( search_result > 0 )
	{
	  for ( n = 0; n < 4; ++n )
	    if ( ! src_grid->mask[src_add[n]] ) search_result = 0;
	}

      /* If point found, find local iw,jw coordinates for weights  */
      if ( search_result > 0 )
	{
	  double iw, jw;  /*  current guess for bilinear coordinate  */

          tgt_grid->cell_frac[dst_add] = ONE;

	  iter = find_ij_weights(plon, plat, src_lats, src_lons, &iw, &jw);

          if ( iter < remap_max_iter )
	    {
	      /* Successfully found iw,jw - compute weights */

	      wgts[0][0] = (ONE-jw*jw*(THREE-TWO*jw)) * (ONE-iw*iw*(THREE-TWO*iw));
	      wgts[0][1] = (ONE-jw*jw*(THREE-TWO*jw)) *      iw*iw*(THREE-TWO*iw);
	      wgts[0][2] =      jw*jw*(THREE-TWO*jw)  *      iw*iw*(THREE-TWO*iw);
	      wgts[0][3] =      jw*jw*(THREE-TWO*jw)  * (ONE-iw*iw*(THREE-TWO*iw));
	      wgts[1][0] = (ONE-jw*jw*(THREE-TWO*jw)) *      iw*(iw-ONE)*(iw-ONE);
	      wgts[1][1] = (ONE-jw*jw*(THREE-TWO*jw)) *      iw*iw*(iw-ONE);
	      wgts[1][2] =      jw*jw*(THREE-TWO*jw)  *      iw*iw*(iw-ONE);
	      wgts[1][3] =      jw*jw*(THREE-TWO*jw)  *      iw*(iw-ONE)*(iw-ONE);
	      wgts[2][0] =      jw*(jw-ONE)*(jw-ONE)  * (ONE-iw*iw*(THREE-TWO*iw));
	      wgts[2][1] =      jw*(jw-ONE)*(jw-ONE)  *      iw*iw*(THREE-TWO*iw);
	      wgts[2][2] =      jw*jw*(jw-ONE)        *      iw*iw*(THREE-TWO*iw);
	      wgts[2][3] =      jw*jw*(jw-ONE)        * (ONE-iw*iw*(THREE-TWO*iw));
	      wgts[3][0] =      iw*(iw-ONE)*(iw-ONE)  *      jw*(jw-ONE)*(jw-ONE);
              wgts[3][1] =      iw*iw*(iw-ONE)        *      jw*(jw-ONE)*(jw-ONE);
	      wgts[3][2] =      iw*iw*(iw-ONE)        *      jw*jw*(jw-ONE);
	      wgts[3][3] =      iw*(iw-ONE)*(iw-ONE)  *      jw*jw*(jw-ONE);

#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_bicub(rv, dst_add, src_add, wgts);
	    }
          else
	    {
	      if ( cdoVerbose || lwarn )
		{
		  lwarn = FALSE;
		  // cdoWarning("Iteration for iw,jw exceed max iteration count of %d!", remap_max_iter);
		  cdoWarning("Bicubic interpolation failed for some grid points - used a distance-weighted average instead!");
		}

	      search_result = -1;
	    }
	}
	  
      /*
	Search for bicubic failed - use a distance-weighted average instead (this is typically near the pole)
      */
      if ( search_result < 0 )
	{
          icount = 0;
          for ( n = 0; n < 4; ++n )
	    {
	      if ( src_grid->mask[src_add[n]] )
		icount++;
	      else
		src_lats[n] = ZERO;
	    }

          if ( icount > 0 )
	    {
	      /* Renormalize weights */
	      double sum_wgts = 0.0; /* sum of weights for normalization */
	      /* 2012-05-08 Uwe Schulzweida: using absolute value of src_lats (bug fix) */
	      for ( n = 0; n < 4; ++n ) sum_wgts += fabs(src_lats[n]);
	      for ( n = 0; n < 4; ++n ) wgts[0][n] = fabs(src_lats[n])/sum_wgts;
	      for ( n = 0; n < 4; ++n ) wgts[1][n] = ZERO;
	      for ( n = 0; n < 4; ++n ) wgts[2][n] = ZERO;
	      for ( n = 0; n < 4; ++n ) wgts[3][n] = ZERO;

	      tgt_grid->cell_frac[dst_add] = ONE;

#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_bicub(rv, dst_add, src_add, wgts);
	    }
        }
    } /* grid_loop1 */

} /* remap_bicub */


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      INTERPOLATION USING A DISTANCE-WEIGHTED AVERAGE                    */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

static
void get_restrict_add(remapgrid_t *src_grid, double plat, const int *restrict src_bin_add, long *minadd, long *maxadd)
{
  long n, n2;
  long min_add = 0, max_add = 0, nm1, np1;
  long nbins;
  restr_t rlat;
  restr_t *bin_lats = src_grid->bin_lats;

  nbins = src_grid->num_srch_bins;

  rlat = RESTR_SCALE(plat);

  for ( n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      if ( rlat >= bin_lats[n2  ] && rlat <= bin_lats[n2+1] )
	{
	  min_add = src_bin_add[n2  ];
	  max_add = src_bin_add[n2+1];

	  nm1 = MAX(n-1, 0);
	  np1 = MIN(n+1, nbins-1);

	  min_add = MIN(min_add, src_bin_add[2*nm1  ]);
	  max_add = MAX(max_add, src_bin_add[2*nm1+1]);
	  min_add = MIN(min_add, src_bin_add[2*np1  ]);
	  max_add = MAX(max_add, src_bin_add[2*np1+1]);
	}
    }

  *minadd = min_add;
  *maxadd = max_add;
  /*
  if ( cdoVerbose )
    printf("plon %g plat %g min_add %ld max_add %ld diff %ld\n",
	   plon, plat, min_add, max_add, max_add-min_add);
  */
}

/*
   This routine finds the closest num_neighbor points to a search 
   point and computes a distance to each of the neighbors.
*/
static
void grid_search_nbr_reg2d(int num_neighbors, remapgrid_t *src_grid, int *restrict nbr_add, double *restrict nbr_dist, 
			   double plat, double plon, const int *restrict src_grid_dims,
			   double coslat_dst, double coslon_dst, double sinlat_dst, double sinlon_dst,
			   const double *restrict sinlat, const double *restrict coslat,
			   const double *restrict sinlon, const double *restrict coslon,
			   const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  /*
    Output variables:

    int nbr_add[num_neighbors]     ! address of each of the closest points
    double nbr_dist[num_neighbors] ! distance to each of the closest points

    Input variables:

    double plat,         ! latitude  of the search point
    double plon,         ! longitude of the search point
  */
  /*  Local variables */
  int lfound;
  long n, nadd, nchk;
  long nx, nxm, ny;
  long ii, jj;
  long i, j, ix;
  int src_add[25];
  long num_add = 0;
  double distance;   //  Angular distance
  /*
  double coslat_dst = cos(plat);  // cos(lat)  of the search point
  double coslon_dst = cos(plon);  // cos(lon)  of the search point
  double sinlat_dst = sin(plat);  // sin(lat)  of the search point
  double sinlon_dst = sin(plon);  // sin(lon)  of the search point
  */
  nx = src_grid_dims[0];
  ny = src_grid_dims[1];

  nxm = nx;
  if ( src_grid->is_cyclic ) nxm++;

  if ( plon < src_center_lon[0]     ) plon += PI2;
  if ( plon > src_center_lon[nxm-1] ) plon -= PI2;

  lfound = rect_grid_search(&ii, &jj, plon, plat, nxm, ny, src_center_lon, src_center_lat);

  if ( lfound )
    {
      if ( src_grid->is_cyclic && ii == (nxm-1) ) ii = 0;

      for ( j = (jj-2); j <= (jj+2); ++j )
	for ( i = (ii-2); i <= (ii+2); ++i )
	  {
	    ix = i;
	    
	    if ( src_grid->is_cyclic )
	      {
		if ( ix <   0 ) ix += nx;
		if ( ix >= nx ) ix -= nx;
	      }

	    if ( ix >= 0 && ix < nx && j >= 0 && j < ny )
	      src_add[num_add++] = j*nx+ix;
	  }
      /*
      num_add = 0;

      for ( j = (jj-1); j <= jj; ++j )
	for ( i = (ii-1); i <= ii; ++i )
	  {
	    ix = i;
	    if ( src_grid->is_cyclic && ix == (nxm-1) ) ix = 0;

	    src_add[num_add++] = j*nx+ix;
	  }
      */
    }

  /* Initialize distance and address arrays */
  for ( n = 0; n < num_neighbors; ++n )
    {
      nbr_add[n]  = -1;
      nbr_dist[n] = BIGNUM;
    }

  if ( lfound )
    {
      long ix, iy;

      for ( long na = 0; na < num_add; ++na )
	{
	  nadd = src_add[na];

	  iy = nadd/nx;
	  ix = nadd - iy*nx;

	  /* Find distance to this point */
	  distance =  sinlat_dst*sinlat[iy] + coslat_dst*coslat[iy]*
	             (coslon_dst*coslon[ix] + sinlon_dst*sinlon[ix]);
	  /*
	  distance =  sinlat_dst*sinlat[nadd] + coslat_dst*coslat[nadd]*
	             (coslon_dst*coslon[nadd] + sinlon_dst*sinlon[nadd]);
	  */
	  /* 2008-07-30 Uwe Schulzweida: check that distance is inside the range of -1 to 1,
	                                 otherwise the result of acos(distance) is NaN */
	  if ( distance >  1 ) distance =  1;
	  if ( distance < -1 ) distance = -1;
	  distance = acos(distance);

	  /* Uwe Schulzweida: if distance is zero, set to small number */
	  if ( IS_EQUAL(distance, 0) ) distance = TINY;

	  /* Store the address and distance if this is one of the smallest four so far */
	  for ( nchk = 0; nchk < num_neighbors; ++nchk )
	    {
	      if ( distance < nbr_dist[nchk] )
		{
		  for ( n = num_neighbors-1; n > nchk; --n )
		    {
		      nbr_add[n]  = nbr_add[n-1];
		      nbr_dist[n] = nbr_dist[n-1];
		    }
		  nbr_add[nchk]  = nadd;
		  nbr_dist[nchk] = distance;
		  break;
		}
	      else if ( num_neighbors == 1 && distance <= nbr_dist[0] && nadd < nbr_add[0] )
		{
		  nbr_add[0]  = nadd;
		  nbr_dist[0] = distance;
		}
	    }
	}
    }
  else if ( src_grid->lextrapolate )
    {
      int search_result;
      search_result = grid_search_reg2d_nn(nx, ny, nbr_add, nbr_dist, plat, plon, src_center_lat, src_center_lon);
      
      if ( search_result >= 0 )
	for ( n = 0; n < 4; ++n ) nbr_add[n] = -1;
    }
}  /*  grid_search_nbr_reg2d  */

static
void grid_search_nbr(int num_neighbors, remapgrid_t *src_grid, int *restrict nbr_add, double *restrict nbr_dist, 
		     double plat, double plon, const int *restrict src_bin_add,
		     double coslat_dst, double coslon_dst, double sinlat_dst, double sinlon_dst,
		     const double *restrict sinlat, const double *restrict coslat,
		     const double *restrict sinlon, const double *restrict coslon)
{
  /*
    Output variables:

    int nbr_add[num_neighbors]     ! address of each of the closest points
    double nbr_dist[num_neighbors] ! distance to each of the closest points

    Input variables:

    int src_bin_add[][2]  ! search bins for restricting search

    double plat,         ! latitude  of the search point
    double plon,         ! longitude of the search point
  */
  /*  Local variables */
  long n, nadd, nchk;
  long min_add, max_add;
  double distance;     /* Angular distance */
  /* result changed a little on a few points with high resolution grid
  double xcoslat_dst = cos(plat);  // cos(lat)  of the search point
  double xcoslon_dst = cos(plon);  // cos(lon)  of the search point
  double xsinlat_dst = sin(plat);  // sin(lat)  of the search point
  double xsinlon_dst = sin(plon);  // sin(lon)  of the search point
  */
  /* Loop over source grid and find nearest neighbors                         */
  /* restrict the search using search bins expand the bins to catch neighbors */

  get_restrict_add(src_grid, plat, src_bin_add, &min_add, &max_add);

  /* Initialize distance and address arrays */
  for ( n = 0; n < num_neighbors; ++n )
    {
      nbr_add[n]  = -1;
      nbr_dist[n] = BIGNUM;
    }

  for ( nadd = min_add; nadd <= max_add; ++nadd )
    {
      /* Find distance to this point */
      distance =  sinlat_dst*sinlat[nadd] + coslat_dst*coslat[nadd]*
	         (coslon_dst*coslon[nadd] + sinlon_dst*sinlon[nadd]);
      /* 2008-07-30 Uwe Schulzweida: check that distance is inside the range of -1 to 1,
                                     otherwise the result of acos(distance) is NaN */
      if ( distance >  1 ) distance =  1;
      if ( distance < -1 ) distance = -1;
      distance = acos(distance);

      /* Uwe Schulzweida: if distance is zero, set to small number */
      if ( IS_EQUAL(distance, 0) ) distance = TINY;

      /* Store the address and distance if this is one of the smallest four so far */
      for ( nchk = 0; nchk < num_neighbors; ++nchk )
	{
          if ( distance < nbr_dist[nchk] )
	    {
	      for ( n = num_neighbors-1; n > nchk; --n )
		{
		  nbr_add[n]  = nbr_add[n-1];
		  nbr_dist[n] = nbr_dist[n-1];
		}
	      nbr_add[nchk]  = nadd;
	      nbr_dist[nchk] = distance;
	      break;
	    }
        }
    }

}  /*  grid_search_nbr  */

/*
  This routine stores the address and weight for this link in the appropriate 
  address and weight arrays and resizes those arrays if necessary.
*/
static
void store_link_nbr(remapvars_t *rv, int add1, int add2, double weights)
{
  /*
    Input variables:
    int  add1         ! address on source grid
    int  add2         ! address on target grid
    double weights    ! remapping weight for this link
  */
  long nlink;

  /*
     Increment number of links and check to see if remap arrays need
     to be increased to accomodate the new link. Then store the link.
  */
  nlink = rv->num_links;
  rv->num_links++;

  if ( rv->num_links >= rv->max_links ) 
    resize_remap_vars(rv, rv->resize_increment);

  rv->src_grid_add[nlink] = add1;
  rv->tgt_grid_add[nlink] = add2;
  rv->wts[nlink]          = weights;

} /* store_link_nbr */

/*
  -----------------------------------------------------------------------

   This routine computes the inverse-distance weights for a
   nearest-neighbor interpolation.

  -----------------------------------------------------------------------
*/
void remap_distwgt(int num_neighbors, remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /*  Local variables */

  long src_grid_size;
  long tgt_grid_size;
  long n;
  long dst_add;                   /* destination address                     */
  int nbr_mask[num_neighbors];    /* mask at nearest neighbors               */
  int nbr_add[num_neighbors];     /* source address at nearest neighbors     */
  double nbr_dist[num_neighbors]; /* angular distance four nearest neighbors */
  double dist_tot;         /* sum of neighbor distances (for normalizing) */
  double coslat_dst;       /* cos(lat) of destination grid point */
  double coslon_dst;       /* cos(lon) of destination grid point */
  double sinlat_dst;       /* sin(lat) of destination grid point */
  double sinlon_dst;       /* sin(lon) of destination grid point */
  double *coslat, *sinlat; /* cosine, sine of grid lats (for distance)    */
  double *coslon, *sinlon; /* cosine, sine of grid lons (for distance)    */
  double wgtstmp;          /* hold the link weight                        */
  double plat, plon;             /*  lat/lon coords of destination point    */
  double findex = 0;
  int remap_grid_type = src_grid->remap_grid_type;

  progressInit();

  /* Compute mappings from source to target grid */

  src_grid_size = src_grid->size;
  tgt_grid_size = tgt_grid->size;

  /* Compute cos, sin of lat/lon on source grid for distance calculations */

  if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
    {
      long nx = src_grid->dims[0];
      long ny = src_grid->dims[1];

      coslat = malloc(ny*sizeof(double));
      coslon = malloc(nx*sizeof(double));
      sinlat = malloc(ny*sizeof(double));
      sinlon = malloc(nx*sizeof(double));

      for ( n = 0; n < nx; ++n )
	{
	  double rlon = src_grid->reg2d_center_lon[n];
	  if ( rlon > PI2  ) rlon -= PI2;
	  if ( rlon < ZERO ) rlon += PI2;
	  coslon[n] = cos(rlon);
	  sinlon[n] = sin(rlon);
	}
      for ( n = 0; n < ny; ++n )
	{
	  coslat[n] = cos(src_grid->reg2d_center_lat[n]);
	  sinlat[n] = sin(src_grid->reg2d_center_lat[n]);
	}
    }
  else
    {
      coslat = malloc(src_grid_size*sizeof(double));
      coslon = malloc(src_grid_size*sizeof(double));
      sinlat = malloc(src_grid_size*sizeof(double));
      sinlon = malloc(src_grid_size*sizeof(double));

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(src_grid, src_grid_size, coslat, coslon, sinlat, sinlon)
#endif
      for ( n = 0; n < src_grid_size; ++n )
	{
	  coslat[n] = cos(src_grid->cell_center_lat[n]);
	  coslon[n] = cos(src_grid->cell_center_lon[n]);
	  sinlat[n] = sin(src_grid->cell_center_lat[n]);
	  sinlon[n] = sin(src_grid->cell_center_lon[n]);
	}
    }

  /* Loop over destination grid  */
#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoTimer, num_neighbors, remap_grid_type, src_grid, tgt_grid, rv, tgt_grid_size, coslat, coslon, sinlat, sinlon, findex) \
  private(dst_add, n, coslat_dst, coslon_dst, sinlat_dst, sinlon_dst, dist_tot, nbr_add, nbr_dist, nbr_mask, wgtstmp, plat, plon) \
  schedule(dynamic,1)
#endif
  for ( dst_add = 0; dst_add < tgt_grid_size; ++dst_add )
    {
      int lprogress = 1;
#if defined(_OPENMP)
      if ( omp_get_thread_num() != 0 ) lprogress = 0;
#endif
#if defined(_OPENMP)
#pragma omp atomic
#endif
      findex++;
      if ( lprogress ) progressStatus(0, 1, findex/tgt_grid_size);

      if ( ! tgt_grid->mask[dst_add] ) continue;
	
      plat = tgt_grid->cell_center_lat[dst_add];
      plon = tgt_grid->cell_center_lon[dst_add];

      coslat_dst = cos(plat);
      coslon_dst = cos(plon);
      sinlat_dst = sin(plat);
      sinlon_dst = sin(plon);

      /* Find nearest grid points on source grid and distances to each point */
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	grid_search_nbr_reg2d(num_neighbors, src_grid, nbr_add, nbr_dist, 
			      plat, plon, src_grid->dims,
			      coslat_dst, coslon_dst, sinlat_dst, sinlon_dst,
			      sinlat, coslat, sinlon, coslon,
			      src_grid->reg2d_center_lat, src_grid->reg2d_center_lon);
      else
	grid_search_nbr(num_neighbors, src_grid, nbr_add, nbr_dist, 
			plat, plon, src_grid->bin_addr,
			coslat_dst, coslon_dst, sinlat_dst, sinlon_dst,
			sinlat, coslat, sinlon, coslon);

      /* Compute weights based on inverse distance if mask is false, eliminate those points */

      dist_tot = ZERO;
      for ( n = 0; n < num_neighbors; ++n )
	{
	  // printf("dst_add %ld %ld %d %g\n", dst_add, n, nbr_add[n], nbr_dist[n]);
	  nbr_mask[n] = FALSE;

	  /* Uwe Schulzweida: check if nbr_add is valid */
	  if ( nbr_add[n] >= 0 )
	    if ( src_grid->mask[nbr_add[n]] )
	      {
		nbr_dist[n] = ONE/nbr_dist[n];
		dist_tot = dist_tot + nbr_dist[n];
		nbr_mask[n] = TRUE;
	      }
	}

      /* Normalize weights and store the link */

      for ( n = 0; n < num_neighbors; ++n )
	{
          if ( nbr_mask[n] )
	    {
	      wgtstmp = nbr_dist[n]/dist_tot;

	      tgt_grid->cell_frac[dst_add] = ONE;
#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_nbr(rv, nbr_add[n], dst_add, wgtstmp);
	    }
	}
    } /* for ( dst_add = 0; dst_add < tgt_grid_size; ++dst_add ) */

  free(coslat);
  free(coslon);
  free(sinlat);
  free(sinlon);

}  /* remap_distwgt */


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      CONSERVATIVE INTERPOLATION                                         */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/*
    This routine is identical to the intersection routine except
    that a coordinate transformation (using a Lambert azimuthal
    equivalent projection) is performed to treat polar cells more
    accurately.
*/
static
void pole_intersection(long *location, double *intrsct_lat, double *intrsct_lon, int *lcoinc,
		       int *lthresh, double beglat, double beglon, double endlat, double endlon,
		       double *begseg, int lrevers,
		       long num_srch_cells, long srch_corners, const int *restrict srch_add,
		       const double *restrict srch_corner_lat, const double *restrict srch_corner_lon,
		       int *luse_last, double *intrsct_x, double *intrsct_y,
		       int *avoid_pole_count, double *avoid_pole_offset)
{
  /*
    Intent(in): 
    double beglat, beglon,  ! beginning lat/lon endpoints for segment
    double endlat, endlon   ! ending    lat/lon endpoints for segment
    int    lrevers          ! flag true if segment integrated in reverse

    Intent(inout) :
    double begseg[2] ! begin lat/lon of full segment
    int *location    ! address in destination array containing this
                     ! segment -- also may contain last location on entry
    int *lthresh     ! flag segment crossing threshold boundary

    intent(out): 
    *int lcoinc      ! flag segment coincident with grid line
    double *intrsct_lat, *intrsct_lon ! lat/lon coords of next intersect.
  */
  /* Local variables */
  long n, next_n, cell;
  long ioffset;
  int loutside; /* flags points outside grid */

  double pi4, rns;           /*  north/south conversion */
  double x1, x2;             /*  local x variables for segment */
  double y1, y2;             /*  local y variables for segment */
  double begx, begy;         /*  beginning x,y variables for segment */
  double endx, endy;         /*  beginning x,y variables for segment */
  double begsegx, begsegy;   /*  beginning x,y variables for segment */
  double grdx1, grdx2;       /*  local x variables for grid cell */
  double grdy1, grdy2;       /*  local y variables for grid cell */
  double vec1_y, vec1_x;     /*  vectors and cross products used */
  double vec2_y, vec2_x;     /*  during grid search */
  double cross_product, eps; /*  eps=small offset away from intersect */
  double s1, s2, determ;     /*  variables used for linear solve to   */
  double mat1, mat2, mat3, mat4, rhs1, rhs2;  /* find intersection */

  double *srch_corner_x;     /*  x of each corner of srch cells */
  double *srch_corner_y;     /*  y of each corner of srch cells */

  /*printf("pole_intersection: %g %g %g %g\n", beglat, beglon, endlat, endlon);*/

  /* Initialize defaults, flags, etc. */

  if ( ! *lthresh ) *location = -1;
  *lcoinc      = FALSE;
  *intrsct_lat = endlat;
  *intrsct_lon = endlon;

  loutside = FALSE;
  s1 = ZERO;

  /* Convert coordinates */

  srch_corner_x = malloc(srch_corners*num_srch_cells*sizeof(double));
  srch_corner_y = malloc(srch_corners*num_srch_cells*sizeof(double));

  if ( beglat > ZERO )
    {
      pi4 = QUART*PI;
      rns = ONE;
    }
  else
    {
      pi4 = -QUART*PI;
      rns = -ONE;
    }

  if ( *luse_last )
    {
      x1 = *intrsct_x;
      y1 = *intrsct_y;
    }
  else
    {
      x1 = rns*TWO*sin(pi4 - HALF*beglat)*cos(beglon);
      y1 =     TWO*sin(pi4 - HALF*beglat)*sin(beglon);
      *luse_last = TRUE;
    }

  x2 = rns*TWO*sin(pi4 - HALF*endlat)*cos(endlon);
  y2 =     TWO*sin(pi4 - HALF*endlat)*sin(endlon);

  for ( n = 0; n < srch_corners*num_srch_cells; ++n )
    {
      srch_corner_x[n] = rns*TWO*sin(pi4 - HALF*srch_corner_lat[n])*
	                         cos(srch_corner_lon[n]);
      srch_corner_y[n] =     TWO*sin(pi4 - HALF*srch_corner_lat[n])*
	                         sin(srch_corner_lon[n]);
    }

  begx = x1;
  begy = y1;
  endx = x2;
  endy = y2;
  begsegx = rns*TWO*sin(pi4 - HALF*begseg[0])*cos(begseg[1]);
  begsegy =     TWO*sin(pi4 - HALF*begseg[0])*sin(begseg[1]);
  *intrsct_x = endx;
  *intrsct_y = endy;

  /*
     Search for location of this segment in ocean grid using cross
     product method to determine whether a point is enclosed by a cell
  */
  while ( TRUE ) /* srch_loop */
    {
      /* If last segment crossed threshold, use that location */

      if ( *lthresh )
	{
	  for ( cell=0; cell < num_srch_cells; ++cell )
	    if ( srch_add[cell] == *location )
	      {
		eps = TINY;
		goto after_srch_loop;
	      }
	}

      /* Otherwise normal search algorithm */

      for ( cell = 0; cell < num_srch_cells; ++cell ) /* cell_loop  */
	{
	  ioffset = cell*srch_corners;
	  for ( n = 0; n < srch_corners; ++n ) /* corner_loop */
	    {
	      next_n = (n+1)%srch_corners;
	      /*
		Here we take the cross product of the vector making 
		up each cell side with the vector formed by the vertex
		and search point.  If all the cross products are 
		positive, the point is contained in the cell.
	      */
	      vec1_x = srch_corner_x[ioffset+next_n] - srch_corner_x[ioffset+n];
	      vec1_y = srch_corner_y[ioffset+next_n] - srch_corner_y[ioffset+n];
	      vec2_x = x1 - srch_corner_x[ioffset+n];
	      vec2_y = y1 - srch_corner_y[ioffset+n];

	      /* If endpoint coincident with vertex, offset the endpoint */

	      if ( IS_EQUAL(vec2_x, 0) && IS_EQUAL(vec2_y, 0) )
		{
		  x1 += 1.e-10*(x2 - x1);
		  y1 += 1.e-10*(y2 - y1);
		  vec2_x = x1 - srch_corner_x[ioffset+n];
		  vec2_y = y1 - srch_corner_y[ioffset+n];
		}

	      cross_product = vec1_x*vec2_y - vec2_x*vec1_y;

	      /*
		If the cross product for a side is ZERO, the point 
                  lies exactly on the side or the length of a side
                  is ZERO.  If the length is ZERO set det > 0.
                  otherwise, perform another cross 
                  product between the side and the segment itself. 
	        If this cross product is also ZERO, the line is 
	          coincident with the cell boundary - perform the 
                  dot product and only choose the cell if the dot 
                  product is positive (parallel vs anti-parallel).
	      */
	      if ( IS_EQUAL(cross_product, 0) )
		{
		  if ( IS_NOT_EQUAL(vec1_x, 0) || IS_NOT_EQUAL(vec1_y, 0) )
		    {
		      vec2_x = x2 - x1;
		      vec2_y = y2 - y1;
		      cross_product = vec1_x*vec2_y - vec2_x*vec1_y;
		    }
		  else
		    cross_product = ONE;

		  if ( IS_EQUAL(cross_product, 0) )
		    {
		      *lcoinc = TRUE;
		      cross_product = vec1_x*vec2_x + vec1_y*vec2_y;
		      if ( lrevers ) cross_product = -cross_product;
		    }
		}

	      /* If cross product is less than ZERO, this cell doesn't work */

	      if ( cross_product < ZERO ) break; /* corner_loop */
	     
	    } /* corner_loop */

	  /* If cross products all positive, we found the location */

	  if  ( n >= srch_corners )
	    {
	      *location = srch_add[cell];
	      /*
		If the beginning of this segment was outside the
		grid, invert the segment so the intersection found
		will be the first intersection with the grid
	      */
	      if ( loutside )
		{
		  x2 = begx;
		  y2 = begy;
		  *location = -1;
		  eps  = -TINY;
		}
	      else
		eps  = TINY;
            
	      goto after_srch_loop;
	    }

	  /* Otherwise move on to next cell */

	} /* cell_loop */

     /*
       If no cell found, the point lies outside the grid.
       take some baby steps along the segment to see if any
       part of the segment lies inside the grid.  
     */
      loutside = TRUE;
      s1 = s1 + BABY_STEP;
      x1 = begx + s1*(x2 - begx);
      y1 = begy + s1*(y2 - begy);

      /* Reached the end of the segment and still outside the grid return no intersection */

      if ( s1 >= ONE )
	{
          free(srch_corner_y);
          free(srch_corner_x);
          *luse_last = FALSE;
          return;
	}
    } /* srch_loop */

 after_srch_loop:

  /*
    Now that a cell is found, search for the next intersection.
    Loop over sides of the cell to find intersection with side
    must check all sides for coincidences or intersections
  */

  ioffset = cell*srch_corners;

  for ( n = 0; n < srch_corners; ++n ) /* intrsct_loop */
    {
      next_n = (n+1)%srch_corners;

      grdy1 = srch_corner_y[ioffset+n];
      grdy2 = srch_corner_y[ioffset+next_n];
      grdx1 = srch_corner_x[ioffset+n];
      grdx2 = srch_corner_x[ioffset+next_n];

      /* Set up linear system to solve for intersection */

      mat1 = x2 - x1;
      mat2 = grdx1 - grdx2;
      mat3 = y2 - y1;
      mat4 = grdy1 - grdy2;
      rhs1 = grdx1 - x1;
      rhs2 = grdy1 - y1;

      determ = mat1*mat4 - mat2*mat3;

      /*
         If the determinant is ZERO, the segments are either 
           parallel or coincident.  Coincidences were detected 
           above so do nothing.
         If the determinant is non-ZERO, solve for the linear 
           parameters s for the intersection point on each line 
           segment.
         If 0<s1,s2<1 then the segment intersects with this side.
           Return the point of intersection (adding a small
           number so the intersection is off the grid line).
      */
      if ( fabs(determ) > 1.e-30 )
	{
          s1 = (rhs1*mat4 - mat2*rhs2)/determ;
          s2 = (mat1*rhs2 - rhs1*mat3)/determ;

	  /* Uwe Schulzweida: s1 >= ZERO! (bug fix) */
          if ( s2 >= ZERO && s2 <= ONE && s1 >= ZERO && s1 <= ONE )
	    {
	      /*
		Recompute intersection based on full segment
		so intersections are consistent for both sweeps
	      */
	      if ( ! loutside )
		{
		  mat1 = x2 - begsegx;
		  mat3 = y2 - begsegy;
		  rhs1 = grdx1 - begsegx;
		  rhs2 = grdy1 - begsegy;
		}
	      else
		{
		  mat1 = x2 - endx;
		  mat3 = y2 - endy;
		  rhs1 = grdx1 - endx;
		  rhs2 = grdy1 - endy;
		}

	      determ = mat1*mat4 - mat2*mat3;

	      /*
		Sometimes due to roundoff, the previous 
		determinant is non-ZERO, but the lines
		are actually coincident.  If this is the
		case, skip the rest.
	      */
	      if ( IS_NOT_EQUAL(determ, 0) )
	       {
		 s1 = (rhs1*mat4 - mat2*rhs2)/determ;
		 s2 = (mat1*rhs2 - rhs1*mat3)/determ;

		 if ( ! loutside )
		   {
		     *intrsct_x = begsegx + s1*mat1;
		     *intrsct_y = begsegy + s1*mat3;
		   }
		 else 
		   {
		     *intrsct_x = endx + s1*mat1;
		     *intrsct_y = endy + s1*mat3;
		   }

		 /* Convert back to lat/lon coordinates */

		 *intrsct_lon = rns*atan2(*intrsct_y, *intrsct_x);
		 if ( *intrsct_lon < ZERO ) 
		   *intrsct_lon = *intrsct_lon + PI2;
		 
		 if ( fabs(*intrsct_x) > 1.e-10 )
		   *intrsct_lat = (pi4 - asin(rns*HALF*(*intrsct_x)/cos(*intrsct_lon)))*TWO;
		 else if ( fabs(*intrsct_y) > 1.e-10 )
		   *intrsct_lat = (pi4 - asin(HALF*(*intrsct_y)/sin(*intrsct_lon)))*TWO;
		 else
		   *intrsct_lat = TWO*pi4;

		 /* Add offset in transformed space for next pass. */

		 if ( s1 - eps/determ < ONE )
		   {
		     *intrsct_x = *intrsct_x - mat1*(eps/determ);
		     *intrsct_y = *intrsct_y - mat3*(eps/determ);
		   }
		 else
		   {
		     if ( ! loutside)
		       {
			 *intrsct_x = endx;
			 *intrsct_y = endy;
			 *intrsct_lat = endlat;
			 *intrsct_lon = endlon;
		       }
		     else 
		       {
			 *intrsct_x = begsegx;
			 *intrsct_y = begsegy;
			 *intrsct_lat = begseg[0];
			 *intrsct_lon = begseg[1];
		       }
		   }

		 break; /* intrsct_loop */
	       }
	    }
	}
 
      /* No intersection this side, move on to next side */

    } /* intrsct_loop */

  free(srch_corner_y);
  free(srch_corner_x);

  /*
     If segment manages to cross over pole, shift the beginning 
     endpoint in order to avoid hitting pole directly
     (it is ok for endpoint to be pole point)
  */

  if ( fabs(*intrsct_x) < 1.e-10 && fabs(*intrsct_y) < 1.e-10 &&
       (IS_NOT_EQUAL(endx, 0) && IS_NOT_EQUAL(endy, 0)) )
    {
      if ( *avoid_pole_count > 2 )
	{
	  *avoid_pole_count  = 0;
	  *avoid_pole_offset = 10.*(*avoid_pole_offset);
        }

      cross_product = begsegx*(endy-begsegy) - begsegy*(endx-begsegx);
      *intrsct_lat = begseg[0];
      if ( cross_product*(*intrsct_lat) > ZERO )
	{
          *intrsct_lon = beglon    + *avoid_pole_offset;
          begseg[1]    = begseg[1] + *avoid_pole_offset;
	}
      else
	{
          *intrsct_lon = beglon    - *avoid_pole_offset;
          begseg[1]    = begseg[1] - *avoid_pole_offset;
        }

      *avoid_pole_count = *avoid_pole_count + 1;
      *luse_last = FALSE;
    }
  else
    {
      *avoid_pole_count  = 0;
      *avoid_pole_offset = TINY;
    }

  /*
     If the segment crosses a pole threshold, reset the intersection
     to be the threshold latitude and do not reuse x,y intersect
     on next entry.  Only check if did not cross threshold last
     time - sometimes the coordinate transformation can place a
     segment on the other side of the threshold again
  */
  if ( *lthresh )
    {
      if ( *intrsct_lat > north_thresh || *intrsct_lat < south_thresh )
	*lthresh = FALSE;
    }
  else if ( beglat > ZERO && *intrsct_lat < north_thresh )
    {
      mat4 = endlat - begseg[0];
      mat3 = endlon - begseg[1];
      if ( mat3 >  PI ) mat3 = mat3 - PI2;
      if ( mat3 < -PI ) mat3 = mat3 + PI2;
      *intrsct_lat = north_thresh - TINY;
      s1 = (north_thresh - begseg[0])/mat4;
      *intrsct_lon = begseg[1] + s1*mat3;
      *luse_last = FALSE;
      *lthresh = TRUE;
    }
  else if ( beglat < ZERO && *intrsct_lat > south_thresh )
    {
      mat4 = endlat - begseg[0];
      mat3 = endlon - begseg[1];
      if ( mat3 >  PI ) mat3 = mat3 - PI2;
      if ( mat3 < -PI ) mat3 = mat3 + PI2;
      *intrsct_lat = south_thresh + TINY;
      s1 = (south_thresh - begseg[0])/mat4;
      *intrsct_lon = begseg[1] + s1*mat3;
      *luse_last = FALSE;
      *lthresh = TRUE;
    }

  /* If reached end of segment, do not use x,y intersect on next entry */

  if ( IS_EQUAL(*intrsct_lat, endlat) && IS_EQUAL(*intrsct_lon, endlon) ) *luse_last = FALSE;

}  /* pole_intersection */


/*
   This routine finds the next intersection of a destination grid line with 
   the line segment given by beglon, endlon, etc.
   A coincidence flag is returned if the segment is entirely coincident with 
   an ocean grid line.  The cells in which to search for an intersection must 
   have already been restricted in the calling routine.
*/
static
void intersection(long *location, double *intrsct_lat, double *intrsct_lon, int *lcoinc,
		  double beglat, double beglon, double endlat, double endlon, double *begseg,
		  int lbegin, int lrevers,
		  long num_srch_cells, long srch_corners, const int *restrict srch_add,
		  const double *restrict srch_corner_lat, const double *restrict srch_corner_lon,
		  int *last_loc, int *lthresh, double *intrsct_lat_off, double *intrsct_lon_off,
		  int *luse_last, double *intrsct_x, double *intrsct_y,
		  int *avoid_pole_count, double *avoid_pole_offset)
{
  /*
    Intent(in): 
    int lbegin,             ! flag for first integration along this segment
    int lrevers             ! flag whether segment integrated in reverse
    double beglat, beglon,  ! beginning lat/lon endpoints for segment
    double endlat, endlon   ! ending    lat/lon endpoints for segment

    Intent(inout) :: 
    double *begseg          ! begin lat/lon of full segment

    intent(out): 
    int *location           ! address in destination array containing this segment
    int *lcoinc             ! flag segments which are entirely coincident with a grid line
    double *intrsct_lat, *intrsct_lon ! lat/lon coords of next intersect.
  */
  /* Local variables */
  long n, next_n, cell;
  long ioffset;

  int  loutside;             /* flags points outside grid */

  double lon1, lon2;         /* local longitude variables for segment */
  double lat1, lat2;         /* local latitude  variables for segment */
  double grdlon1, grdlon2;   /* local longitude variables for grid cell */
  double grdlat1, grdlat2;   /* local latitude  variables for grid cell */
  double vec1_lat, vec1_lon; /* vectors and cross products used */
  double vec2_lat, vec2_lon; /* during grid search */
  double cross_product; 
  double eps, offset;        /* small offset away from intersect */
  double s1, s2, determ;     /* variables used for linear solve to */
  double mat1 = 0, mat2, mat3 = 0, mat4, rhs1, rhs2;  /* find intersection */

  /* Initialize defaults, flags, etc. */

  *location    = -1;
  *lcoinc      = FALSE;
  *intrsct_lat = endlat;
  *intrsct_lon = endlon;

  if ( num_srch_cells == 0 ) return;

  if ( beglat > north_thresh || beglat < south_thresh )
    {
      if ( *lthresh ) *location = *last_loc;
      pole_intersection(location,
			intrsct_lat, intrsct_lon, lcoinc, lthresh,
			beglat, beglon, endlat, endlon, begseg, lrevers,
			num_srch_cells, srch_corners, srch_add,
			srch_corner_lat, srch_corner_lon,
			luse_last, intrsct_x, intrsct_y,
			avoid_pole_count, avoid_pole_offset);

      if ( *lthresh )
	{
          *last_loc = *location;
          *intrsct_lat_off = *intrsct_lat;
          *intrsct_lon_off = *intrsct_lon;
        }
      return;
    }

  loutside = FALSE;
  if ( lbegin )
    {
      lat1 = beglat;
      lon1 = beglon;
    }
  else
    {
      lat1 = *intrsct_lat_off;
      lon1 = *intrsct_lon_off;
    }

  lat2 = endlat;
  lon2 = endlon;
  if      ( (lon2-lon1) >  THREE*PIH ) lon2 -= PI2;
  else if ( (lon2-lon1) < -THREE*PIH ) lon2 += PI2;

  s1 = ZERO;

  /*
     Search for location of this segment in ocean grid using cross
     product method to determine whether a point is enclosed by a cell
  */
  while ( TRUE ) /* srch_loop */
    {
      /* If last segment crossed threshold, use that location */

      if ( *lthresh )
       {
         for ( cell = 0; cell < num_srch_cells; ++cell )
	   if ( srch_add[cell] == *last_loc )
	     {
               *location = *last_loc;
               eps = TINY;
               goto after_srch_loop;
	     }
       }

      /* Otherwise normal search algorithm */

      for ( cell = 0; cell < num_srch_cells; ++cell ) /* cell_loop  */
	{
	  ioffset = cell*srch_corners;
	  for ( n = 0; n < srch_corners; n++ ) /* corner_loop */
	    {
	      next_n = (n+1)%srch_corners;
	      /*
		Here we take the cross product of the vector making 
		up each cell side with the vector formed by the vertex
		and search point.  If all the cross products are 
		positive, the point is contained in the cell.
	      */
	      vec1_lat = srch_corner_lat[ioffset+next_n] - srch_corner_lat[ioffset+n];
	      vec1_lon = srch_corner_lon[ioffset+next_n] - srch_corner_lon[ioffset+n];
	      vec2_lat = lat1 - srch_corner_lat[ioffset+n];
	      vec2_lon = lon1 - srch_corner_lon[ioffset+n];

	      /* If endpoint coincident with vertex, offset the endpoint */

	      if ( IS_EQUAL(vec2_lat, 0) && IS_EQUAL(vec2_lon, 0) )
		{
		  lat1 += 1.e-10*(lat2-lat1);
		  lon1 += 1.e-10*(lon2-lon1);
		  vec2_lat = lat1 - srch_corner_lat[ioffset+n];
		  vec2_lon = lon1 - srch_corner_lon[ioffset+n];
		}

	      /* Check for 0,2pi crossings */

	      if      ( vec1_lon >  PI ) vec1_lon -= PI2;
	      else if ( vec1_lon < -PI ) vec1_lon += PI2;

	      if      ( vec2_lon >  PI ) vec2_lon -= PI2;
	      else if ( vec2_lon < -PI ) vec2_lon += PI2;

	      cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat;

	      /*
	       If the cross product for a side is ZERO, the point 
                 lies exactly on the side or the side is degenerate
                 (ZERO length).  If degenerate, set the cross 
                 product to a positive number.  Otherwise perform 
                 another cross product between the side and the 
                 segment itself. 
	       If this cross product is also ZERO, the line is 
                 coincident with the cell boundary - perform the 
                 dot product and only choose the cell if the dot 
                 product is positive (parallel vs anti-parallel).
	      */
	      if ( IS_EQUAL(cross_product, 0) )
		{
		  if ( IS_NOT_EQUAL(vec1_lat, 0) || IS_NOT_EQUAL(vec1_lon, 0) )
		    {
		      vec2_lat = lat2 - lat1;
		      vec2_lon = lon2 - lon1;

		      if      ( vec2_lon >  PI ) vec2_lon -= PI2;
		      else if ( vec2_lon < -PI ) vec2_lon += PI2;

		      cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat;
		    }
		  else
		    cross_product = ONE;

		  if ( IS_EQUAL(cross_product, 0) )
		    {
		      *lcoinc = TRUE;
		      cross_product = vec1_lon*vec2_lon + vec1_lat*vec2_lat;
		      if ( lrevers ) cross_product = -cross_product;
		    }
		}

	      /* If cross product is less than ZERO, this cell doesn't work */

	      if ( cross_product < ZERO ) break; /* corner_loop */

	    } /* corner_loop */

	  /* If cross products all positive, we found the location */

	  if ( n >= srch_corners )
	    {
	      *location = srch_add[cell];
	      /*
		If the beginning of this segment was outside the
		grid, invert the segment so the intersection found
		will be the first intersection with the grid
	      */
	      if ( loutside )
		{
		  lat2 = beglat;
		  lon2 = beglon;
		  *location = -1;
		  eps  = -TINY;
		}
	      else
		eps  = TINY;

	      goto after_srch_loop;
	    }

	  /* Otherwise move on to next cell */

	} /* cell_loop */

      /*
	If still no cell found, the point lies outside the grid.
	Take some baby steps along the segment to see if any
	part of the segment lies inside the grid.  
      */
      loutside = TRUE;
      s1 = s1 + BABY_STEP;
      lat1 = beglat + s1*(endlat - beglat);
      lon1 = beglon + s1*(lon2   - beglon);

      /* Reached the end of the segment and still outside the grid return no intersection */

      if ( s1 >= ONE ) return;

    } /* srch_loop */

 after_srch_loop:

  /*
    Now that a cell is found, search for the next intersection.
    Loop over sides of the cell to find intersection with side
    must check all sides for coincidences or intersections
  */

  ioffset = cell*srch_corners;

  for ( n = 0; n < srch_corners; ++n ) /* intrsct_loop */
    {
      next_n = (n+1)%srch_corners;

      grdlon1 = srch_corner_lon[ioffset+n];
      grdlon2 = srch_corner_lon[ioffset+next_n];
      grdlat1 = srch_corner_lat[ioffset+n];
      grdlat2 = srch_corner_lat[ioffset+next_n];

      /* Set up linear system to solve for intersection */

      mat1 = lat2 - lat1;
      mat2 = grdlat1 - grdlat2;
      mat3 = lon2 - lon1;
      mat4 = grdlon1 - grdlon2;
      rhs1 = grdlat1 - lat1;
      rhs2 = grdlon1 - lon1;

      if      ( mat3 >  PI ) mat3 -= PI2;
      else if ( mat3 < -PI ) mat3 += PI2;

      if      ( mat4 >  PI ) mat4 -= PI2;
      else if ( mat4 < -PI ) mat4 += PI2;

      if      ( rhs2 >  PI ) rhs2 -= PI2;
      else if ( rhs2 < -PI ) rhs2 += PI2;

      determ = mat1*mat4 - mat2*mat3;

      /*
         If the determinant is ZERO, the segments are either 
           parallel or coincident.  Coincidences were detected 
           above so do nothing.
         If the determinant is non-ZERO, solve for the linear 
           parameters s for the intersection point on each line 
           segment.
         If 0<s1,s2<1 then the segment intersects with this side.
           Return the point of intersection (adding a small
           number so the intersection is off the grid line).
      */
      if ( fabs(determ) > 1.e-30 )
	{
	  s1 = (rhs1*mat4 - mat2*rhs2)/determ;
	  s2 = (mat1*rhs2 - rhs1*mat3)/determ;

	  if ( s2 >= ZERO && s2 <= ONE && s1 >= ZERO && s1 <= ONE )
	    {
	      /*
		Recompute intersection based on full segment
		so intersections are consistent for both sweeps
	      */
	      if ( ! loutside )
		{
		  mat1 = lat2 - begseg[0];
		  mat3 = lon2 - begseg[1];
		  rhs1 = grdlat1 - begseg[0];
		  rhs2 = grdlon1 - begseg[1];
		}
	      else
		{
		  mat1 = begseg[0] - endlat;
		  mat3 = begseg[1] - endlon;
		  rhs1 = grdlat1 - endlat;
		  rhs2 = grdlon1 - endlon;
		}

	      if      ( mat3 >  PI ) mat3 -= PI2;
	      else if ( mat3 < -PI ) mat3 += PI2;

	      if      ( rhs2 > PI  ) rhs2 -= PI2;
	      else if ( rhs2 < -PI ) rhs2 += PI2;

	      determ = mat1*mat4 - mat2*mat3;

	      /*
		Sometimes due to roundoff, the previous 
		determinant is non-ZERO, but the lines
		are actually coincident.  If this is the
		case, skip the rest.
	      */
	      if ( IS_NOT_EQUAL(determ, 0) )
		{
		  s1 = (rhs1*mat4 - mat2*rhs2)/determ;
		  s2 = (mat1*rhs2 - rhs1*mat3)/determ;

		  offset = s1 + eps/determ;
		  if ( offset > ONE ) offset = ONE;

		  if ( ! loutside )
		    {
		      *intrsct_lat = begseg[0] + mat1*s1;
		      *intrsct_lon = begseg[1] + mat3*s1;
		      *intrsct_lat_off = begseg[0] + mat1*offset;
		      *intrsct_lon_off = begseg[1] + mat3*offset;
		    }
		  else
		    {
		      *intrsct_lat = endlat + mat1*s1;
		      *intrsct_lon = endlon + mat3*s1;
		      *intrsct_lat_off = endlat + mat1*offset;
		      *intrsct_lon_off = endlon + mat3*offset;
		    }
		  break; /* intrsct_loop */
		}
	    }
	}

      /* No intersection this side, move on to next side */

    } /* intrsct_loop */

  /*
     If the segment crosses a pole threshold, reset the intersection
     to be the threshold latitude.  Only check if this was not a
     threshold segment since sometimes coordinate transform can end
     up on other side of threshold again.
  */
  if ( *lthresh )
    {
      if ( *intrsct_lat < north_thresh || *intrsct_lat > south_thresh )
	*lthresh = FALSE;
    }
  else if ( lat1 > ZERO && *intrsct_lat > north_thresh )
    {
      *intrsct_lat = north_thresh + TINY;
      *intrsct_lat_off = north_thresh + eps*mat1;
      s1 = (*intrsct_lat - begseg[0])/mat1;
      *intrsct_lon     = begseg[1] + s1*mat3;
      *intrsct_lon_off = begseg[1] + (s1+eps)*mat3;
      *last_loc = *location;
      *lthresh = TRUE;
    }
  else if ( lat1 < ZERO && *intrsct_lat < south_thresh )
    {
      *intrsct_lat = south_thresh - TINY;
      *intrsct_lat_off = south_thresh + eps*mat1;
      s1 = (*intrsct_lat - begseg[0])/mat1;
      *intrsct_lon     = begseg[1] + s1*mat3;
      *intrsct_lon_off = begseg[1] + (s1+eps)*mat3;
      *last_loc = *location;
      *lthresh = TRUE;
    }

}  /* intersection */


/*
   This routine computes the line integral of the flux function 
   that results in the interpolation weights.  The line is defined
   by the input lat/lon of the endpoints.
*/
static
void line_integral(double *weights, double in_phi1, double in_phi2, 
		   double theta1, double theta2, double grid1_lon, double grid2_lon)
{
  /*
    Intent(in): 
    double in_phi1, in_phi2,     ! Longitude endpoints for the segment
    double theta1, theta2,       ! Latitude  endpoints for the segment
    double grid1_lon,            ! Reference coordinates for each
    double grid2_lon             ! Grid (to ensure correct 0,2pi interv.)

    Intent(out):
    double weights[6]            ! Line integral contribution to weights
  */

  /*  Local variables  */
  double dphi, sinth1, sinth2, costh1, costh2, fac;
  double phi1, phi2;
  double f1, f2, fint;

  /*  Weights for the general case based on a trapezoidal approx to the integrals. */

  sinth1 = sin(theta1);
  sinth2 = sin(theta2);
  costh1 = cos(theta1);
  costh2 = cos(theta2);

  dphi = in_phi1 - in_phi2;
  if      ( dphi >  PI ) dphi -= PI2;
  else if ( dphi < -PI ) dphi += PI2;
      
  dphi = HALF*dphi;

  /*
     The first weight is the area overlap integral. The second and
     fourth are second-order latitude gradient weights.
  */
  weights[0] = dphi*(sinth1 + sinth2);
  weights[1] = dphi*(costh1 + costh2 + (theta1*sinth1 + theta2*sinth2));
  weights[3] = weights[0];
  weights[4] = weights[1];

  /*
     The third and fifth weights are for the second-order phi gradient
     component.  Must be careful of longitude range.
  */
  f1 = HALF*(costh1*sinth1 + theta1);
  f2 = HALF*(costh2*sinth2 + theta2);

  phi1 = in_phi1 - grid1_lon;
  if      ( phi1 >  PI ) phi1 -= PI2;
  else if ( phi1 < -PI ) phi1 += PI2;

  phi2 = in_phi2 - grid1_lon;
  if      ( phi2 >  PI ) phi2 -= PI2;
  else if ( phi2 < -PI ) phi2 += PI2;

  if ( (phi2-phi1) <  PI && (phi2-phi1) > -PI )
    weights[2] = dphi*(phi1*f1 + phi2*f2);
  else
    {
      if ( phi1 > ZERO )
	fac = PI;
      else
	fac = -PI;

      fint = f1 + (f2-f1)*(fac-phi1)/fabs(dphi);
      weights[2] = HALF*phi1*(phi1-fac)*f1 -
	           HALF*phi2*(phi2+fac)*f2 +
	           HALF*fac*(phi1+phi2)*fint;
    }

  phi1 = in_phi1 - grid2_lon;
  if      ( phi1 >  PI ) phi1 -= PI2;
  else if ( phi1 < -PI ) phi1 += PI2;

  phi2 = in_phi2 - grid2_lon;
  if      ( phi2 >  PI ) phi2 -= PI2;
  else if ( phi2 < -PI ) phi2 += PI2;

  if ( (phi2-phi1) <  PI  && (phi2-phi1) > -PI )
    weights[5] = dphi*(phi1*f1 + phi2*f2);
  else
    {
      if ( phi1 > ZERO ) fac =  PI;
      else               fac = -PI;

      fint = f1 + (f2-f1)*(fac-phi1)/fabs(dphi);
      weights[5] = HALF*phi1*(phi1-fac)*f1 -
     	           HALF*phi2*(phi2+fac)*f2 +
	           HALF*fac*(phi1+phi2)*fint;
    }

}  /* line_integral */

static
void grid_store_init(grid_store_t *grid_store, long gridsize)
{
  long iblk;
  long blksize[] = {128, 256, 512, 1024, 2048, 4096, 8192};
  long nblks = sizeof(blksize)/sizeof(long);
  long nblocks;

  for ( iblk = nblks-1; iblk >= 0; --iblk )
    if ( gridsize/blksize[iblk] > 99 ) break;

  if ( iblk < 0 ) iblk = 0;

  /* grid_store->blk_size = BLK_SIZE; */
  grid_store->blk_size = blksize[iblk];
  grid_store->max_size = gridsize;

  grid_store->nblocks = grid_store->max_size/grid_store->blk_size;
  if ( grid_store->max_size%grid_store->blk_size > 0 ) grid_store->nblocks++;

  if ( cdoVerbose )
    fprintf(stdout, "blksize = %d  lastblksize = %d  max_size = %d  nblocks = %d\n", 
	    grid_store->blk_size, grid_store->max_size%grid_store->blk_size, 
	    grid_store->max_size, grid_store->nblocks);

  grid_store->blksize = malloc(grid_store->nblocks*sizeof(int));
  grid_store->nlayers = malloc(grid_store->nblocks*sizeof(int));
  grid_store->layers  = malloc(grid_store->nblocks*sizeof(grid_layer_t *));

  nblocks = grid_store->nblocks;
  for ( iblk = 0; iblk < nblocks; ++iblk )
    {
      grid_store->blksize[iblk] = grid_store->blk_size;
      grid_store->nlayers[iblk] = 0;
      grid_store->layers[iblk]  = NULL;
    }
  if ( grid_store->max_size%grid_store->blk_size > 0 )
    grid_store->blksize[grid_store->nblocks-1] = grid_store->max_size%grid_store->blk_size;
}

static
void grid_store_delete(grid_store_t *grid_store)
{
  grid_layer_t *grid_layer, *grid_layer_f;
  long ilayer;
  long i, j;
  long iblk;

  for ( iblk = 0; iblk < grid_store->nblocks; ++iblk )
    {
      j = 0;
      grid_layer = grid_store->layers[iblk];
      for ( ilayer = 0; ilayer < grid_store->nlayers[iblk]; ++ilayer )
	{
	  if ( cdoVerbose )
	    {
	      for ( i = 0; i < grid_store->blksize[iblk]; ++i )
		if ( grid_layer->grid2_link[i] != -1 ) j++;
	    }
	      
	  grid_layer_f = grid_layer;
	  free(grid_layer->grid2_link);
	  grid_layer = grid_layer->next;
	  free(grid_layer_f);
	}

      if ( cdoVerbose )
	{
	  fprintf(stderr, "block = %ld nlayers = %d  allocated = %d  used = %ld\n",
		  iblk+1, grid_store->nlayers[iblk], 
		  grid_store->nlayers[iblk]*grid_store->blksize[iblk], j);
	}
    }

  free(grid_store->blksize);
  free(grid_store->layers);
  free(grid_store->nlayers);  
}

/*
    This routine stores the address and weight for this link in the appropriate 
    address and weight arrays and resizes those arrays if necessary.
*/
static
void store_link_cnsrv_fast(remapvars_t *rv, long add1, long add2, long num_wts, double *weights, grid_store_t *grid_store)
{
  /*
    Input variables:
    int  add1         ! address on source grid
    int  add2         ! address on target grid
    double weights[]  ! array of remapping weights for this link
  */
  /* Local variables */
  long nlink; /* link index */
  long ilayer, i, iblk, iadd2;
  long nlayer, blksize;
  int lstore_link;
  grid_layer_t *grid_layer, **grid_layer2;

  /*  If all weights are ZERO, do not bother storing the link */

  if ( num_wts == 3 )
    {
      if ( IS_EQUAL(weights[0], 0) && IS_EQUAL(weights[1], 0) && IS_EQUAL(weights[2], 0) ) return;
    }
  else
    {
      if ( IS_EQUAL(weights[0], 0) ) return;
    }
    
  /* If the link already exists, add the weight to the current weight arrays */

  iblk  = BLK_NUM(add2);
  iadd2 = BLK_IDX(add2);

  lstore_link = FALSE;
  grid_layer2 = &grid_store->layers[iblk];
  nlayer = grid_store->nlayers[iblk];
  for ( ilayer = 0; ilayer < nlayer; ++ilayer )
    {
      grid_layer = *grid_layer2;
      nlink = grid_layer->grid2_link[iadd2];
      if ( nlink == -1 )
	{
	  break;
	}
      else if ( add1 == rv->src_grid_add[nlink] )
	{
	  lstore_link = TRUE;
	  break;
	}
      grid_layer2 = &(*grid_layer2)->next;
    }

  if ( lstore_link )
    {
      for ( i = 0; i < num_wts; ++i ) rv->wts[num_wts*nlink+i] += weights[i];	      
      return;
    }

  /*
     If the link does not yet exist, increment number of links and 
     check to see if remap arrays need to be increased to accomodate 
     the new link. Then store the link.
  */
  nlink = rv->num_links;

  if ( ilayer < grid_store->nlayers[iblk] )
    {
      grid_layer->grid2_link[iadd2] = nlink;
    }
  else
    {
      grid_layer = malloc(sizeof(grid_layer_t));
      grid_layer->next = NULL;
      grid_layer->grid2_link = malloc(grid_store->blksize[iblk]*sizeof(int));

      blksize = grid_store->blksize[iblk];
      for ( i = 0; i < blksize; ++i )
	grid_layer->grid2_link[i] = -1;

      grid_layer->grid2_link[iadd2] = nlink;
      *grid_layer2 = grid_layer;
      grid_store->nlayers[iblk]++;
    }

  rv->num_links++;
  if ( rv->num_links >= rv->max_links )
    resize_remap_vars(rv, rv->resize_increment);

  rv->src_grid_add[nlink] = add1;
  rv->tgt_grid_add[nlink] = add2;

  for ( i = 0; i < num_wts; ++i ) rv->wts[num_wts*nlink+i] = weights[i];	      

}  /* store_link_cnsrv_fast */


/*
    This routine stores the address and weight for this link in the appropriate 
    address and weight arrays and resizes those arrays if necessary.
*/
static
void store_link_cnsrv(remapvars_t *rv, long add1, long add2, double *restrict weights,
		      int *link_add1[2], int *link_add2[2])
{
  /*
    Input variables:
    int  add1         ! address on source grid
    int  add2         ! address on target grid
    double weights[3] ! array of remapping weights for this link
  */
  /* Local variables */
  long nlink, min_link, max_link; /* link index */

  /*  If all weights are ZERO, do not bother storing the link */

  if ( IS_EQUAL(weights[0], 0) && IS_EQUAL(weights[1], 0) && IS_EQUAL(weights[2], 0) ) return;

  /*  Restrict the range of links to search for existing links */

  min_link = MIN(link_add1[0][add1], link_add2[0][add2]);
  max_link = MAX(link_add1[1][add1], link_add2[1][add2]);
  if ( min_link == -1 )
    {
      min_link = 0;
      max_link = -1;
    }

  /* If the link already exists, add the weight to the current weight arrays */

#if defined(SX)
#define STRIPED 1
#if STRIPED
#define STRIPLENGTH 4096
  {
    long ilink = max_link + 1;
    long strip, estrip;
    nlink = 0;
    for ( strip=min_link; strip <= max_link; strip+=STRIPLENGTH )
      {
	estrip = MIN(max_link-strip+1, STRIPLENGTH);
	for ( nlink = 0; nlink < estrip; ++nlink )
	  {
	    if ( add2 == rv->tgt_grid_add[strip+nlink] &&
		 add1 == rv->src_grid_add[strip+nlink] )
	      ilink = strip + nlink;
	  }
	if (ilink != (max_link + 1)) break;
      }
    nlink += strip;
    if (ilink != (max_link + 1)) nlink = ilink;
  }
#else
  {
    long ilink = max_link + 1;
    for ( nlink = min_link; nlink <= max_link; ++nlink )
      {
	if ( add2 == rv->tgt_grid_add[nlink] )
	  if ( add1 == rv->src_grid_add[nlink] ) ilink = nlink;
      }
    if ( ilink != (max_link + 1) ) nlink = ilink;
  }
#endif
#else
  for ( nlink = min_link; nlink <= max_link; ++nlink )
    {
      if ( add2 == rv->tgt_grid_add[nlink] )
	if ( add1 == rv->src_grid_add[nlink] ) break;
    }
#endif

  if ( nlink <= max_link )
    {
      rv->wts[3*nlink  ] += weights[0];
      rv->wts[3*nlink+1] += weights[1];
      rv->wts[3*nlink+2] += weights[2];

      return;
    }

  /*
     If the link does not yet exist, increment number of links and 
     check to see if remap arrays need to be increased to accomodate 
     the new link. Then store the link.
  */
  nlink = rv->num_links;

  rv->num_links++;
  if ( rv->num_links >= rv->max_links )
    resize_remap_vars(rv, rv->resize_increment);

  rv->src_grid_add[nlink] = add1;
  rv->tgt_grid_add[nlink] = add2;

  rv->wts[3*nlink  ] = weights[0];
  rv->wts[3*nlink+1] = weights[1];
  rv->wts[3*nlink+2] = weights[2];

  if ( link_add1[0][add1] == -1 ) link_add1[0][add1] = (int)nlink;
  if ( link_add2[0][add2] == -1 ) link_add2[0][add2] = (int)nlink;
  link_add1[1][add1] = (int)nlink;
  link_add2[1][add2] = (int)nlink;

}  /* store_link_cnsrv */

int rect_grid_search2(long *imin, long *imax, double xmin, double xmax, long nxm, const double *restrict xm);

static
long get_srch_cells_reg2d(const int *restrict src_grid_dims, 
			  const double *restrict src_corner_lat, const double *restrict src_corner_lon,
			  const double *restrict tgt_cell_bound_box, int *srch_add)
{
  long nx = src_grid_dims[0];
  long ny = src_grid_dims[1];
  long num_srch_cells;  /* num cells in restricted search arrays   */
  int lfound;
  long nxp1, nyp1;
  double src_lon_min, src_lon_max;
  int debug = 0;

  nxp1 = nx+1;
  nyp1 = ny+1;

  src_lon_min = src_corner_lon[0];
  src_lon_max = src_corner_lon[nx];

  double bound_lon1, bound_lon2;

  num_srch_cells = 0;

  long imin = nxp1, imax = -1, jmin = nyp1, jmax = -1;
  long im, jm;

  lfound = rect_grid_search2(&jmin, &jmax, tgt_cell_bound_box[0], tgt_cell_bound_box[1], nyp1, src_corner_lat);
  if ( jmin > 0 ) jmin--;
  if ( jmax < (ny-2) ) jmax++;
  bound_lon1 = tgt_cell_bound_box[2];
  bound_lon2 = tgt_cell_bound_box[3];
  if ( bound_lon1 <= src_lon_max && bound_lon2 >= src_lon_min )
    {
      if ( debug ) printf("  b1 %g %g\n", bound_lon1*RAD2DEG, bound_lon2*RAD2DEG);
      if ( bound_lon1 < src_lon_min && bound_lon2 > src_lon_min ) bound_lon1 = src_lon_min;
      if ( bound_lon2 > src_lon_max && bound_lon1 < src_lon_max ) bound_lon2 = src_lon_max;
      lfound = rect_grid_search2(&imin, &imax, bound_lon1, bound_lon2, nxp1, src_corner_lon);
      if ( imin != -1 )
	{
	  if ( debug )
	    printf("   %g %g imin %ld  imax %ld  jmin %ld jmax %ld\n", RAD2DEG*src_corner_lon[imin], RAD2DEG*src_corner_lon[imax+1], imin, imax, jmin, jmax);
	  for ( jm = jmin; jm <= jmax; ++jm )
	    for ( im = imin; im <= imax; ++im )
	      srch_add[num_srch_cells++] = jm*nx + im;
	}
    }

  bound_lon1 = tgt_cell_bound_box[2];
  bound_lon2 = tgt_cell_bound_box[3];
  if ( bound_lon1 <= src_lon_min && bound_lon2 >= src_lon_min )
    {
      long imin2 = nxp1, imax2 = -1;
      bound_lon1 += 2*M_PI;
      bound_lon2 += 2*M_PI;
      if ( debug ) printf("  b2 %g %g\n", bound_lon1*RAD2DEG, bound_lon2*RAD2DEG);
      if ( bound_lon1 < src_lon_min && bound_lon2 > src_lon_min ) bound_lon1 = src_lon_min;
      if ( bound_lon2 > src_lon_max && bound_lon1 < src_lon_max ) bound_lon2 = src_lon_max;
      lfound = rect_grid_search2(&imin2, &imax2, bound_lon1, bound_lon2, nxp1, src_corner_lon);
      if ( imin2 != -1 && imin2 == imax ) imin2 += 1;
      if ( imax2 != -1 && imax2 == imax ) imax2 += 1;
      if ( imin2 != -1 )
	{
	  if ( debug )
	    printf("   %g %g imin %ld  imax %ld  jmin %ld jmax %ld\n", RAD2DEG*src_corner_lon[imin2], RAD2DEG*src_corner_lon[imax2+1], imin2, imax2, jmin, jmax);
	  for ( jm = jmin; jm <= jmax; ++jm )
	    for ( im = imin2; im <= imax2; ++im )
	      srch_add[num_srch_cells++] = jm*nx + im;
	}
    }

  bound_lon1 = tgt_cell_bound_box[2];
  bound_lon2 = tgt_cell_bound_box[3];
  if ( bound_lon1 <= src_lon_max && bound_lon2 >= src_lon_max )
    {
      long imin3 = nxp1, imax3 = -1;
      bound_lon1 -= 2*M_PI;
      bound_lon2 -= 2*M_PI;
      if ( debug ) printf("  b3 %g %g\n", bound_lon1*RAD2DEG, bound_lon2*RAD2DEG);
      if ( bound_lon1 < src_lon_min && bound_lon2 > src_lon_min ) bound_lon1 = src_lon_min;
      if ( bound_lon2 > src_lon_max && bound_lon1 < src_lon_max ) bound_lon2 = src_lon_max;
      lfound = rect_grid_search2(&imin3, &imax3, bound_lon1, bound_lon2, nxp1, src_corner_lon);
      if ( imin3 != -1 && imin3 == imin ) imin3 -= 1;
      if ( imax3 != -1 && imax3 == imin ) imax3 -= 1;
      if ( imin3 != -1 )
	{
	  if ( debug )
	    printf("   %g %g imin %ld  imax %ld  jmin %ld jmax %ld\n", RAD2DEG*src_corner_lon[imin3], RAD2DEG*src_corner_lon[imax3+1], imin3, imax3, jmin, jmax);
	  for ( jm = jmin; jm <= jmax; ++jm )
	    for ( im = imin3; im <= imax3; ++im )
	      srch_add[num_srch_cells++] = jm*nx + im;
	}
    }

  return (num_srch_cells);
}

static
long get_srch_cells(long tgt_grid_add, long nbins, int *bin_addr1, int *bin_addr2,
		    restr_t *tgt_cell_bound_box, restr_t *src_cell_bound_box, long src_grid_size, int *srch_add)
{
  long num_srch_cells;  /* num cells in restricted search arrays   */
  long min_add;         /* addresses for restricting search of     */
  long max_add;         /* destination grid                        */
  long n, n2;           /* generic counters                        */
  long src_grid_add;    /* current linear address for src cell     */
  long src_grid_addm4;
  restr_t bound_box_lat1, bound_box_lat2, bound_box_lon1, bound_box_lon2;

  /* Restrict searches first using search bins */

  min_add = src_grid_size - 1;
  max_add = 0;

  for ( n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      if ( tgt_grid_add >= bin_addr1[n2] && tgt_grid_add <= bin_addr1[n2+1] )
	{
	  if ( bin_addr2[n2  ] < min_add ) min_add = bin_addr2[n2  ];
	  if ( bin_addr2[n2+1] > max_add ) max_add = bin_addr2[n2+1];
	}
    }

  /* Further restrict searches using bounding boxes */

  bound_box_lat1 = tgt_cell_bound_box[0];
  bound_box_lat2 = tgt_cell_bound_box[1];
  bound_box_lon1 = tgt_cell_bound_box[2];
  bound_box_lon2 = tgt_cell_bound_box[3];

  num_srch_cells = 0;
  for ( src_grid_add = min_add; src_grid_add <= max_add; ++src_grid_add )
    {
      src_grid_addm4 = src_grid_add<<2;
      if ( (src_cell_bound_box[src_grid_addm4+2] <= bound_box_lon2)  &&
	   (src_cell_bound_box[src_grid_addm4+3] >= bound_box_lon1) )
	{
	  if ( (src_cell_bound_box[src_grid_addm4  ] <= bound_box_lat2)  &&
	       (src_cell_bound_box[src_grid_addm4+1] >= bound_box_lat1) )
	    {
	      srch_add[num_srch_cells] = src_grid_add;
	      num_srch_cells++;
	    }
	}
    }

  if ( bound_box_lon1 < RESTR_SCALE(0.) || bound_box_lon2 > RESTR_SCALE(PI2) )
    {
      if ( bound_box_lon1 < RESTR_SCALE(0.) )
	{
	  bound_box_lon1 += RESTR_SCALE(PI2);
	  bound_box_lon2 += RESTR_SCALE(PI2);
	}
      else
	{
	  bound_box_lon1 -= RESTR_SCALE(PI2);
	  bound_box_lon2 -= RESTR_SCALE(PI2);
	}

      for ( src_grid_add = min_add; src_grid_add <= max_add; ++src_grid_add )
	{
	  src_grid_addm4 = src_grid_add<<2;
	  if ( (src_cell_bound_box[src_grid_addm4+2] <= bound_box_lon2)  &&
	       (src_cell_bound_box[src_grid_addm4+3] >= bound_box_lon1) )
	    {
	      if ( (src_cell_bound_box[src_grid_addm4  ] <= bound_box_lat2)  &&
		   (src_cell_bound_box[src_grid_addm4+1] >= bound_box_lat1) )
		{
		  long ii;
		  for ( ii = 0; ii < num_srch_cells; ++ii )
		    if ( srch_add[ii] == src_grid_add ) break;
		  
		  if ( ii == num_srch_cells )
		    {
		      srch_add[num_srch_cells] = src_grid_add;
		      num_srch_cells++;
		    }
		}
	    }
	}
    }

  return (num_srch_cells);
}


/*
  -----------------------------------------------------------------------

   This routine traces the perimeters of every grid cell on each
   grid checking for intersections with the other grid and computing
   line integrals for each subsegment.

  -----------------------------------------------------------------------
*/
void remap_conserv(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /* local variables */

  int lcheck = TRUE;

  long ioffset;
  long max_subseg = 100000; /* max number of subsegments per segment to prevent infinite loop */
                            /* 1000 is too small!!! */
  long src_grid_size;
  long tgt_grid_size;
  long src_num_cell_corners;
  long tgt_num_cell_corners;
  long src_grid_add;       /* current linear address for source grid cell   */
  long tgt_grid_add;       /* current linear address for target grid cell   */
  long n, n3, k;        /* generic counters                        */
  long corner;          /* corner of cell that segment starts from */
  long next_corn;       /* corner of cell that segment ends on     */
  long nbins, num_links;
  long num_subseg;      /* number of subsegments                   */

  int lcoinc;           /* flag for coincident segments            */
  int lrevers;          /* flag for reversing direction of segment */
  int lbegin;           /* flag for first integration of a segment */

  double intrsct_lat, intrsct_lon;         /* lat/lon of next intersect  */
  double beglat, endlat, beglon, endlon;   /* endpoints of current seg.  */
  double norm_factor = 0;                  /* factor for normalizing wts */

  double *tgt_centroid_lat, *tgt_centroid_lon;   /* centroid coords  */
  double *src_centroid_lat, *src_centroid_lon;   /* on each grid     */

  double begseg[2];         /* begin lat/lon for full segment */
  double weights[6];        /* local wgt array */
  long    num_wts;

  long    max_srch_cells;   /* num cells in restricted search arrays  */
  long    num_srch_cells;   /* num cells in restricted search arrays  */
  long    srch_corners;     /* num of corners of srch cells           */
  long    nsrch_corners;
  int*    srch_add;         /* global address of cells in srch arrays */
  int*    srch_add2[ompNumThreads];
  int     ompthID, i;
  double *srch_corner_lat;  /* lat of each corner of srch cells */
  double *srch_corner_lon;  /* lon of each corner of srch cells */

  int *link_add1[2];        /* min,max link add to restrict search */
  int *link_add2[2];        /* min,max link add to restrict search */

  /* Intersection */
  int last_loc = -1;        /* save location when crossing threshold  */
  int lthresh = FALSE;      /* flags segments crossing threshold bndy */
  double intrsct_lat_off = 0, intrsct_lon_off = 0; /* lat/lon coords offset for next search */

  /* Pole_intersection */
  /* Save last intersection to avoid roundoff during coord transformation */
  int luse_last = FALSE;
  double intrsct_x, intrsct_y;      /* x,y for intersection */
  /* Variables necessary if segment manages to hit pole */
  int avoid_pole_count = 0;         /* count attempts to avoid pole  */
  double avoid_pole_offset = TINY;  /* endpoint offset to avoid pole */
  grid_store_t *grid_store = NULL;
  double findex = 0;

  progressInit();

  nbins = src_grid->num_srch_bins;
  num_wts = rv->num_wts;

  if ( remap_store_link_fast )
    {
      grid_store = malloc(sizeof(grid_store_t));
      grid_store_init(grid_store, tgt_grid->size);
    }

  if ( cdoVerbose )
    {
      cdoPrint("north_thresh: %g", north_thresh);
      cdoPrint("south_thresh: %g", south_thresh);
    }

  if ( cdoTimer ) timer_start(timer_remap_con);

  src_grid_size = src_grid->size;
  tgt_grid_size = tgt_grid->size;

  src_num_cell_corners = src_grid->num_cell_corners;
  tgt_num_cell_corners = tgt_grid->num_cell_corners;

  if ( ! remap_store_link_fast )
    {
      link_add1[0] = malloc(src_grid_size*sizeof(int));
      link_add1[1] = malloc(src_grid_size*sizeof(int));
      link_add2[0] = malloc(tgt_grid_size*sizeof(int));
      link_add2[1] = malloc(tgt_grid_size*sizeof(int));

#if defined(SX)
#pragma vdir nodep
#endif
      for ( n = 0; n < src_grid_size; ++n )
	{
	  link_add1[0][n] = -1;
	  link_add1[1][n] = -1;
	}

#if defined(SX)
#pragma vdir nodep
#endif
      for ( n = 0; n < tgt_grid_size; ++n )
	{
	  link_add2[0][n] = -1;
	  link_add2[1][n] = -1;
	}
    }

  /* Initialize centroid arrays */

  src_centroid_lat = malloc(src_grid_size*sizeof(double));
  src_centroid_lon = malloc(src_grid_size*sizeof(double));
  tgt_centroid_lat = malloc(tgt_grid_size*sizeof(double));
  tgt_centroid_lon = malloc(tgt_grid_size*sizeof(double));

  for ( n = 0; n < src_grid_size; ++n )
    {
      src_centroid_lat[n] = 0;
      src_centroid_lon[n] = 0;
    }

  for ( n = 0; n < tgt_grid_size; ++n )
    {
      tgt_centroid_lat[n] = 0;
      tgt_centroid_lon[n] = 0;
    }

  double* srch_corner_lat2[ompNumThreads];
  double* srch_corner_lon2[ompNumThreads];
  long max_srch_cells2[ompNumThreads];

  /*  Integrate around each cell on source grid */

  for ( i = 0; i < ompNumThreads; ++i )
    {
      srch_corner_lat2[i] = NULL;
      srch_corner_lon2[i] = NULL;
    }

  for ( i = 0; i < ompNumThreads; ++i )
    max_srch_cells2[i] = 0;

  for ( i = 0; i < ompNumThreads; ++i )
    srch_add2[i] = malloc(tgt_grid_size*sizeof(int));

  srch_corners    = tgt_num_cell_corners;

  if ( cdoTimer ) timer_start(timer_remap_con_l1);

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoTimer, nbins, num_wts, src_centroid_lon, src_centroid_lat, \
         remap_store_link_fast, grid_store, link_add1, link_add2, rv, cdoVerbose, max_subseg, \
	 srch_corner_lat2, srch_corner_lon2, max_srch_cells2, 		\
	 src_num_cell_corners,	srch_corners, src_grid, tgt_grid, tgt_grid_size, src_grid_size, srch_add2, findex) \
  private(ompthID, srch_add, n, k, num_srch_cells, max_srch_cells, 	\
	  src_grid_add, tgt_grid_add, ioffset, nsrch_corners, corner, next_corn, beglat, beglon, \
	  endlat, endlon, lrevers, begseg, lbegin, num_subseg, srch_corner_lat, srch_corner_lon, \
	  weights, intrsct_lat, intrsct_lon, intrsct_lat_off, intrsct_lon_off, intrsct_x, intrsct_y, \
	  last_loc, lcoinc, lthresh, luse_last, avoid_pole_count, avoid_pole_offset)
#endif
  for ( src_grid_add = 0; src_grid_add < src_grid_size; ++src_grid_add )
    {
#if defined(_OPENMP)
      ompthID = omp_get_thread_num();
#else
      ompthID = 0;
#endif

      int lprogress = 1;
      if ( ompthID != 0 ) lprogress = 0;

#if defined(_OPENMP)
#pragma omp atomic
#endif
      findex++;
      if ( lprogress ) progressStatus(0, 0.5, findex/src_grid_size);

      srch_add = srch_add2[ompthID];

      lthresh   = FALSE;
      luse_last = FALSE;
      avoid_pole_count  = 0;
      avoid_pole_offset = TINY;

      /* Get search cells */
      num_srch_cells = get_srch_cells(src_grid_add, nbins, src_grid->bin_addr, tgt_grid->bin_addr,
				      src_grid->cell_bound_box+src_grid_add*4, tgt_grid->cell_bound_box, tgt_grid_size, srch_add);

      if ( num_srch_cells == 0 ) continue;

      /* Create search arrays */

      max_srch_cells  = max_srch_cells2[ompthID];
      srch_corner_lat = srch_corner_lat2[ompthID];
      srch_corner_lon = srch_corner_lon2[ompthID];

      if ( num_srch_cells > max_srch_cells )
	{
	  srch_corner_lat = realloc(srch_corner_lat, srch_corners*num_srch_cells*sizeof(double));
	  srch_corner_lon = realloc(srch_corner_lon, srch_corners*num_srch_cells*sizeof(double));

	  max_srch_cells  = num_srch_cells;

	  max_srch_cells2[ompthID]  = max_srch_cells;
	  srch_corner_lat2[ompthID] = srch_corner_lat;
	  srch_corner_lon2[ompthID] = srch_corner_lon;
	}

      /* gather1 */
      for ( n = 0; n < num_srch_cells; ++n )
	{
	  tgt_grid_add = srch_add[n];
	  ioffset = tgt_grid_add*srch_corners;

	  nsrch_corners = n*srch_corners;
	  for ( k = 0; k < srch_corners; k++ )
	    {
	      srch_corner_lat[nsrch_corners+k] = tgt_grid->cell_corner_lat[ioffset+k];
	      srch_corner_lon[nsrch_corners+k] = tgt_grid->cell_corner_lon[ioffset+k];
	    }
	}

      /* Integrate around this cell */

      ioffset = src_grid_add*src_num_cell_corners;

      for ( corner = 0; corner < src_num_cell_corners; ++corner )
	{
          next_corn = (corner+1)%src_num_cell_corners;

          /* Define endpoints of the current segment */

          beglat = src_grid->cell_corner_lat[ioffset+corner];
          beglon = src_grid->cell_corner_lon[ioffset+corner];
          endlat = src_grid->cell_corner_lat[ioffset+next_corn];
          endlon = src_grid->cell_corner_lon[ioffset+next_corn];
          lrevers = FALSE;

	  /*  To ensure exact path taken during both sweeps, always integrate segments in the same direction (SW to NE). */
          if ( (endlat < beglat) || (IS_EQUAL(endlat, beglat) && endlon < beglon) )
	    {
	      beglat = src_grid->cell_corner_lat[ioffset+next_corn];
	      beglon = src_grid->cell_corner_lon[ioffset+next_corn];
	      endlat = src_grid->cell_corner_lat[ioffset+corner];
	      endlon = src_grid->cell_corner_lon[ioffset+corner];
	      lrevers = TRUE;
	    }

          begseg[0] = beglat;
          begseg[1] = beglon;
          lbegin = TRUE;

          /*
	    If this is a constant-longitude segment, skip the rest 
	    since the line integral contribution will be ZERO.
          */
          if ( IS_NOT_EQUAL(endlon, beglon) )
	    {
	      num_subseg = 0;
	      /*
		Integrate along this segment, detecting intersections 
		and computing the line integral for each sub-segment
	      */
	      while ( IS_NOT_EQUAL(beglat, endlat) || IS_NOT_EQUAL(beglon, endlon) )
		{
		  /*  Prevent infinite loops if integration gets stuck near cell or threshold boundary */
		  num_subseg++;
		  if ( num_subseg >= max_subseg )
		    cdoAbort("Integration stalled: num_subseg exceeded limit (grid1[%d]: lon1=%g lon2=%g lat1=%g lat2=%g)!",
			     src_grid_add, beglon, endlon, beglat, endlat);

		  /* Uwe Schulzweida: skip very small regions */
		  if ( num_subseg%1000 == 0 )
		    {
		      if ( fabs(beglat-endlat) < 1.e-10 || fabs(beglon-endlon) < 1.e-10 )
			{
			  if ( cdoVerbose )
			    cdoPrint("Skip very small region (grid1[%d]): lon=%g dlon=%g lat=%g dlat=%g",
				     src_grid_add, beglon, endlon-beglon, beglat, endlat-beglat);
			  break;
			}
		    }

		  /* Find next intersection of this segment with a gridline on grid 2. */

		  intersection(&tgt_grid_add, &intrsct_lat, &intrsct_lon, &lcoinc,
			       beglat, beglon, endlat, endlon, begseg, 
			       lbegin, lrevers,
                               num_srch_cells, srch_corners, srch_add,
			       srch_corner_lat, srch_corner_lon,
			       &last_loc, &lthresh, &intrsct_lat_off, &intrsct_lon_off,
			       &luse_last, &intrsct_x, &intrsct_y,
			       &avoid_pole_count, &avoid_pole_offset);

		  lbegin = FALSE;

		  /* Compute line integral for this subsegment. */

		  if ( tgt_grid_add != -1 )
		    line_integral(weights, beglon, intrsct_lon, beglat, intrsct_lat,
				  src_grid->cell_center_lon[src_grid_add], tgt_grid->cell_center_lon[tgt_grid_add]);
		  else
		    line_integral(weights, beglon, intrsct_lon, beglat, intrsct_lat,
				  src_grid->cell_center_lon[src_grid_add], src_grid->cell_center_lon[src_grid_add]);

		  /* If integrating in reverse order, change sign of weights */

		  if ( lrevers ) for ( k = 0; k < 6; ++k ) weights[k] = -weights[k];

		  /*
		    Store the appropriate addresses and weights. 
		    Also add contributions to cell areas and centroids.
		  */
		  if ( tgt_grid_add != -1 )
		    if ( src_grid->mask[src_grid_add] )
		      {
#if defined(_OPENMP)
#pragma omp critical
#endif
			{
			  if ( remap_store_link_fast )
			    store_link_cnsrv_fast(rv, src_grid_add, tgt_grid_add, num_wts, weights, grid_store);
			  else
			    store_link_cnsrv(rv, src_grid_add, tgt_grid_add, weights, link_add1, link_add2);

			  tgt_grid->cell_frac[tgt_grid_add] += weights[3];
			}
			src_grid->cell_frac[src_grid_add] += weights[0];
		      }

		  src_grid->cell_area[src_grid_add] += weights[0];
		  src_centroid_lat[src_grid_add] += weights[1];
		  src_centroid_lon[src_grid_add] += weights[2];

		  /* Reset beglat and beglon for next subsegment. */
		  beglat = intrsct_lat;
		  beglon = intrsct_lon;
		}
	    }
          /* End of segment */
        }
    }

  if ( cdoTimer ) timer_stop(timer_remap_con_l1);

  /* Finished with all cells: deallocate search arrays */

  for ( i = 0; i < ompNumThreads; ++i )
    {
      free(srch_corner_lon2[i]);
      free(srch_corner_lat2[i]);
    }

  for ( i = 0; i < ompNumThreads; ++i )
    free(srch_add2[i]);

  /* Integrate around each cell on target grid */

  for ( i = 0; i < ompNumThreads; ++i )
    {
      srch_corner_lat2[i] = NULL;
      srch_corner_lon2[i] = NULL;
    }

  for ( i = 0; i < ompNumThreads; ++i )
    max_srch_cells2[i] = 0;

  for ( i = 0; i < ompNumThreads; ++i )
    srch_add2[i] = malloc(src_grid_size*sizeof(int));

  srch_corners    = src_num_cell_corners;
  max_srch_cells  = 0;
  srch_corner_lat = NULL;
  srch_corner_lon = NULL;

  if ( cdoTimer ) timer_start(timer_remap_con_l2);

  findex = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoTimer, nbins, num_wts, tgt_centroid_lon, tgt_centroid_lat, \
         remap_store_link_fast, grid_store, link_add1, link_add2, rv, cdoVerbose, max_subseg, \
	 srch_corner_lat2, srch_corner_lon2, max_srch_cells2, 		\
	 tgt_num_cell_corners, srch_corners, src_grid, tgt_grid, tgt_grid_size, src_grid_size, srch_add2, findex) \
  private(ompthID, srch_add, n, k, num_srch_cells, max_srch_cells,	\
	  src_grid_add, tgt_grid_add, ioffset, nsrch_corners, corner, next_corn, beglat, beglon, \
	  endlat, endlon, lrevers, begseg, lbegin, num_subseg, srch_corner_lat, srch_corner_lon, \
	  weights, intrsct_lat, intrsct_lon, intrsct_lat_off, intrsct_lon_off, intrsct_x, intrsct_y, \
	  last_loc, lcoinc, lthresh, luse_last, avoid_pole_count, avoid_pole_offset)
#endif
  for ( tgt_grid_add = 0; tgt_grid_add < tgt_grid_size; ++tgt_grid_add )
    {
#if defined(_OPENMP)
      ompthID = omp_get_thread_num();
#else
      ompthID = 0;
#endif

      int lprogress = 1;
      if ( ompthID != 0 ) lprogress = 0;

#if defined(_OPENMP)
#pragma omp atomic
#endif
      findex++;
      if ( lprogress ) progressStatus(0.5, 0.5, findex/tgt_grid_size);

      srch_add = srch_add2[ompthID];

      lthresh   = FALSE;
      luse_last = FALSE;
      avoid_pole_count  = 0;
      avoid_pole_offset = TINY;

      /* Get search cells */
      num_srch_cells = get_srch_cells(tgt_grid_add, nbins, tgt_grid->bin_addr, src_grid->bin_addr,
				      tgt_grid->cell_bound_box+tgt_grid_add*4, src_grid->cell_bound_box, src_grid_size, srch_add);

      if ( num_srch_cells == 0 ) continue;

      /* Create search arrays */
      
      max_srch_cells  = max_srch_cells2[ompthID];
      srch_corner_lat = srch_corner_lat2[ompthID];
      srch_corner_lon = srch_corner_lon2[ompthID];

      if ( num_srch_cells > max_srch_cells )
	{
	  srch_corner_lat = realloc(srch_corner_lat, srch_corners*num_srch_cells*sizeof(double));
	  srch_corner_lon = realloc(srch_corner_lon, srch_corners*num_srch_cells*sizeof(double));

	  max_srch_cells  = num_srch_cells;

	  max_srch_cells2[ompthID]  = max_srch_cells;
	  srch_corner_lat2[ompthID] = srch_corner_lat;
	  srch_corner_lon2[ompthID] = srch_corner_lon;
	}

      /* gather2 */
      for ( n = 0; n < num_srch_cells; ++n )
	{
	  src_grid_add = srch_add[n];
	  ioffset = src_grid_add*srch_corners;

	  nsrch_corners = n*srch_corners;
	  for ( k = 0; k < srch_corners; ++k )
	    {
	      srch_corner_lat[nsrch_corners+k] = src_grid->cell_corner_lat[ioffset+k];
	      srch_corner_lon[nsrch_corners+k] = src_grid->cell_corner_lon[ioffset+k];
	    }
	}

      /* Integrate around this cell */

      ioffset = tgt_grid_add*tgt_num_cell_corners;

      for ( corner = 0; corner < tgt_num_cell_corners; ++corner )
	{
          next_corn = (corner+1)%tgt_num_cell_corners;

          beglat = tgt_grid->cell_corner_lat[ioffset+corner];
          beglon = tgt_grid->cell_corner_lon[ioffset+corner];
          endlat = tgt_grid->cell_corner_lat[ioffset+next_corn];
          endlon = tgt_grid->cell_corner_lon[ioffset+next_corn];
          lrevers = FALSE;

	  /* To ensure exact path taken during both sweeps, always integrate in the same direction */
          if ( (endlat < beglat) || (IS_EQUAL(endlat, beglat) && endlon < beglon) )
	    {
	      beglat = tgt_grid->cell_corner_lat[ioffset+next_corn];
	      beglon = tgt_grid->cell_corner_lon[ioffset+next_corn];
	      endlat = tgt_grid->cell_corner_lat[ioffset+corner];
	      endlon = tgt_grid->cell_corner_lon[ioffset+corner];
	      lrevers = TRUE;
	    }

          begseg[0] = beglat;
          begseg[1] = beglon;
          lbegin = TRUE;

          /*
	    If this is a constant-longitude segment, skip the rest 
	    since the line integral contribution will be ZERO.
          */
          if ( IS_NOT_EQUAL(endlon, beglon) )
	    {
	      num_subseg = 0;
	      /*
		Integrate along this segment, detecting intersections 
		and computing the line integral for each sub-segment
	      */
	      while ( IS_NOT_EQUAL(beglat, endlat) || IS_NOT_EQUAL(beglon, endlon) )
		{
		  /*  Prevent infinite loops if integration gets stuck near cell or threshold boundary */
		  num_subseg++;
		  if ( num_subseg >= max_subseg )
		    cdoAbort("Integration stalled: num_subseg exceeded limit (grid2[%d]: lon1=%g lon2=%g lat1=%g lat2=%g)!",
			     tgt_grid_add, beglon, endlon, beglat, endlat);

		  /* Uwe Schulzweida: skip very small regions */
		  if ( num_subseg%1000 == 0 )
		    {
		      if ( fabs(beglat-endlat) < 1.e-10 || fabs(beglon-endlon) < 1.e-10 )
			{
			  if ( cdoVerbose )
			    cdoPrint("Skip very small region (grid2[%d]): lon=%g dlon=%g lat=%g dlat=%g",
				     tgt_grid_add, beglon, endlon-beglon, beglat, endlat-beglat);
			  break;
			}
		    }

		  /* Find next intersection of this segment with a gridline on grid 2. */

		  intersection(&src_grid_add, &intrsct_lat, &intrsct_lon, &lcoinc,
			       beglat, beglon, endlat, endlon, begseg,
			       lbegin, lrevers,
                               num_srch_cells, srch_corners, srch_add,
			       srch_corner_lat, srch_corner_lon,
			       &last_loc, &lthresh, &intrsct_lat_off, &intrsct_lon_off,
			       &luse_last, &intrsct_x, &intrsct_y,
			       &avoid_pole_count, &avoid_pole_offset);

		  lbegin = FALSE;

		  /* Compute line integral for this subsegment. */

		  if ( src_grid_add != -1 )
		    line_integral(weights, beglon, intrsct_lon, beglat, intrsct_lat,
				  src_grid->cell_center_lon[src_grid_add], tgt_grid->cell_center_lon[tgt_grid_add]);
		  else
		    line_integral(weights, beglon, intrsct_lon, beglat, intrsct_lat,
				  tgt_grid->cell_center_lon[tgt_grid_add], tgt_grid->cell_center_lon[tgt_grid_add]);

		  /* If integrating in reverse order, change sign of weights */

		  if ( lrevers ) for ( k = 0; k < 6; ++k ) weights[k] = -weights[k];

		  /*
		    Store the appropriate addresses and weights. 
		    Also add contributions to cell areas and centroids.
		    If there is a coincidence, do not store weights
		    because they have been captured in the previous loop.
		    The source grid mask is the master mask
		  */
		  if ( ! lcoinc && src_grid_add != -1 )
		    if ( src_grid->mask[src_grid_add] )
		      {
#if defined(_OPENMP)
#pragma omp critical
#endif
			{
			  if ( remap_store_link_fast )
			    store_link_cnsrv_fast(rv, src_grid_add, tgt_grid_add, num_wts, weights, grid_store);
			  else
			    store_link_cnsrv(rv, src_grid_add, tgt_grid_add, weights, link_add1, link_add2);

			  src_grid->cell_frac[src_grid_add] += weights[0];
			}
			tgt_grid->cell_frac[tgt_grid_add] += weights[3];
		      }

		  tgt_grid->cell_area[tgt_grid_add] += weights[3];
		  tgt_centroid_lat[tgt_grid_add] += weights[4];
		  tgt_centroid_lon[tgt_grid_add] += weights[5];

		  /* Reset beglat and beglon for next subsegment. */
		  beglat = intrsct_lat;
		  beglon = intrsct_lon;
		}
	    }
          /* End of segment */
	}
    }

  if ( cdoTimer ) timer_stop(timer_remap_con_l2);

  /* Finished with all cells: deallocate search arrays */

  for ( i = 0; i < ompNumThreads; ++i )
    {
      free(srch_corner_lon2[i]);
      free(srch_corner_lat2[i]);
    }

  for ( i = 0; i < ompNumThreads; ++i )
    free(srch_add2[i]);

  /*
     Correct for situations where N/S pole not explicitly included in
     grid (i.e. as a grid corner point). If pole is missing from only
     one grid, need to correct only the area and centroid of that 
     grid.  If missing from both, do complete weight calculation.
  */

  /* North Pole */
  weights[0] =  PI2;
  weights[1] =  PI*PI;
  weights[2] =  ZERO;
  weights[3] =  PI2;
  weights[4] =  PI*PI;
  weights[5] =  ZERO;

  src_grid_add = -1;
  /* pole_loop1 */
  for ( n = 0; n < src_grid_size; ++n )
    if ( src_grid->cell_area[n] < -THREE*PIH && src_grid->cell_center_lat[n] > ZERO )
      {
	src_grid_add = n;
#ifndef SX
	break;
#endif
      }

  tgt_grid_add = -1;
  /* pole_loop2 */
  for ( n = 0; n < tgt_grid_size; ++n )
    if ( tgt_grid->cell_area[n] < -THREE*PIH && tgt_grid->cell_center_lat[n] > ZERO )
      {
	tgt_grid_add = n;
#ifndef SX
	break;
#endif
      }

  if ( src_grid_add != -1 )
    {
      src_grid->cell_area[src_grid_add]     += weights[0];
      src_centroid_lat[src_grid_add] += weights[1];
      src_centroid_lon[src_grid_add] += weights[2];
    }

  if ( tgt_grid_add != -1 )
    {
      tgt_grid->cell_area[tgt_grid_add]     += weights[3];
      tgt_centroid_lat[tgt_grid_add] += weights[4];
      tgt_centroid_lon[tgt_grid_add] += weights[5];
    }

  if ( src_grid_add != -1 && tgt_grid_add != -1 )
    {
      if ( remap_store_link_fast )
	store_link_cnsrv_fast(rv, src_grid_add, tgt_grid_add, num_wts, weights, grid_store);
      else
	store_link_cnsrv(rv, src_grid_add, tgt_grid_add, weights, link_add1, link_add2);

      src_grid->cell_frac[src_grid_add] += weights[0];
      tgt_grid->cell_frac[tgt_grid_add] += weights[3];
    }

  /* South Pole */
  weights[0] =  PI2;
  weights[1] = -PI*PI;
  weights[2] =  ZERO;
  weights[3] =  PI2;
  weights[4] = -PI*PI;
  weights[5] =  ZERO;

  src_grid_add = -1;
  /* pole_loop3 */
  for ( n = 0; n < src_grid_size; ++n )
    if ( src_grid->cell_area[n] < -THREE*PIH && src_grid->cell_center_lat[n] < ZERO )
      {
	src_grid_add = n;
#ifndef SX
	break;
#endif
      }

  tgt_grid_add = -1;
  /* pole_loop4 */
  for ( n = 0; n < tgt_grid_size; ++n )
    if ( tgt_grid->cell_area[n] < -THREE*PIH && tgt_grid->cell_center_lat[n] < ZERO )
      {
	tgt_grid_add = n;
#ifndef SX
	break;
#endif
      }

  if ( src_grid_add != -1 )
    {
      src_grid->cell_area[src_grid_add]     += weights[0];
      src_centroid_lat[src_grid_add] += weights[1];
      src_centroid_lon[src_grid_add] += weights[2];
    }

  if ( tgt_grid_add != -1 )
    {
      tgt_grid->cell_area[tgt_grid_add]     += weights[3];
      tgt_centroid_lat[tgt_grid_add] += weights[4];
      tgt_centroid_lon[tgt_grid_add] += weights[5];
    }

  if ( src_grid_add != -1 && tgt_grid_add != -1 )
    {
      if ( remap_store_link_fast )
	store_link_cnsrv_fast(rv, src_grid_add, tgt_grid_add, num_wts, weights, grid_store);
      else
	store_link_cnsrv(rv, src_grid_add, tgt_grid_add, weights, link_add1, link_add2);

      src_grid->cell_frac[src_grid_add] += weights[0];
      tgt_grid->cell_frac[tgt_grid_add] += weights[3];
    }

  if ( remap_store_link_fast )
    {
      grid_store_delete(grid_store);
      free(grid_store);
    }


  /* Finish centroid computation */

  for ( n = 0; n < src_grid_size; ++n )
    if ( IS_NOT_EQUAL(src_grid->cell_area[n], 0) )
      {
        src_centroid_lat[n] /= src_grid->cell_area[n];
        src_centroid_lon[n] /= src_grid->cell_area[n];
      }

  for ( n = 0; n < tgt_grid_size; ++n )
    if ( IS_NOT_EQUAL(tgt_grid->cell_area[n], 0) )
      {
        tgt_centroid_lat[n] /= tgt_grid->cell_area[n];
        tgt_centroid_lon[n] /= tgt_grid->cell_area[n];
      }

  /* 2010-10-08 Uwe Schulzweida: remove all links with weights < 0 */

  /* 
  if ( 1 )
    {
      num_links = rv->num_links;

      if ( cdoVerbose )
	for ( n = 0; n < num_links; n++ )
	  printf("wts1: %d %g\n", n, rv->wts[3*n]);

      for ( n = 0; n < num_links; n++ )
	{
	  if ( rv->wts[3*n] < 0 )
	    {
	      int i, n2, nd;
     
	      for ( n2 = n+1; n2 < num_links; n2++ )
		if ( rv->wts[3*n2] >= 0 ) break;

	      nd = n2-n;
	      num_links -= nd;
	      for ( i = n; i < num_links; i++ )
		{
		  rv->wts[3*i]   = rv->wts[3*(i+nd)];
		  rv->wts[3*i+1] = rv->wts[3*(i+nd)+1];
		  rv->wts[3*i+2] = rv->wts[3*(i+nd)+2];
		  
		  rv->src_grid_add[i] = rv->src_grid_add[i+nd];
		  rv->tgt_grid_add[i] = rv->tgt_grid_add[i+nd];
		}
	    }
	}

     if ( cdoVerbose ) cdoPrint("Removed number of links = %ld", rv->num_links - num_links);

      rv->num_links = num_links;
    }
  */

  /* Include centroids in weights and normalize using destination area if requested */

  num_links = rv->num_links;

  if ( rv->norm_opt == NORM_OPT_DESTAREA )
    {
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(num_links, rv, tgt_grid, src_centroid_lat, src_centroid_lon)		\
  private(n, n3, src_grid_add, tgt_grid_add, weights, norm_factor)
#endif
      for ( n = 0; n < num_links; ++n )
	{
	  n3 = n*3;
	  src_grid_add = rv->src_grid_add[n]; tgt_grid_add = rv->tgt_grid_add[n];
	  weights[0] = rv->wts[n3]; weights[1] = rv->wts[n3+1]; weights[2] = rv->wts[n3+2];

          if ( IS_NOT_EQUAL(tgt_grid->cell_area[tgt_grid_add], 0) )
	    norm_factor = ONE/tgt_grid->cell_area[tgt_grid_add];
          else
            norm_factor = ZERO;

	  rv->wts[n3  ] =  weights[0]*norm_factor;
	  rv->wts[n3+1] = (weights[1] - weights[0]*src_centroid_lat[src_grid_add])*norm_factor;
	  rv->wts[n3+2] = (weights[2] - weights[0]*src_centroid_lon[src_grid_add])*norm_factor;
	}
    }
  else if ( rv->norm_opt == NORM_OPT_FRACAREA )
    {
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(num_links, rv, tgt_grid, src_centroid_lat, src_centroid_lon)		\
  private(n, n3, src_grid_add, tgt_grid_add, weights, norm_factor)
#endif
      for ( n = 0; n < num_links; ++n )
	{
	  n3 = n*3;
	  src_grid_add = rv->src_grid_add[n]; tgt_grid_add = rv->tgt_grid_add[n];
	  weights[0] = rv->wts[n3]; weights[1] = rv->wts[n3+1]; weights[2] = rv->wts[n3+2];

          if ( IS_NOT_EQUAL(tgt_grid->cell_frac[tgt_grid_add], 0) )
	    norm_factor = ONE/tgt_grid->cell_frac[tgt_grid_add];
          else
            norm_factor = ZERO;

	  rv->wts[n3  ] =  weights[0]*norm_factor;
	  rv->wts[n3+1] = (weights[1] - weights[0]*src_centroid_lat[src_grid_add])*norm_factor;
	  rv->wts[n3+2] = (weights[2] - weights[0]*src_centroid_lon[src_grid_add])*norm_factor;
	}
    }
  else if ( rv->norm_opt == NORM_OPT_NONE )
    {
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(num_links, rv, tgt_grid, src_centroid_lat, src_centroid_lon)	\
  private(n, n3, src_grid_add, tgt_grid_add, weights, norm_factor)
#endif
      for ( n = 0; n < num_links; ++n )
	{
	  n3 = n*3;
	  src_grid_add = rv->src_grid_add[n]; tgt_grid_add = rv->tgt_grid_add[n];
	  weights[0] = rv->wts[n3]; weights[1] = rv->wts[n3+1]; weights[2] = rv->wts[n3+2];

          norm_factor = ONE;

	  rv->wts[n3  ] =  weights[0]*norm_factor;
	  rv->wts[n3+1] = (weights[1] - weights[0]*src_centroid_lat[src_grid_add])*norm_factor;
	  rv->wts[n3+2] = (weights[2] - weights[0]*src_centroid_lon[src_grid_add])*norm_factor;
	}
    }

  if ( cdoVerbose )
    cdoPrint("Total number of links = %ld", rv->num_links);

  for ( n = 0; n < src_grid_size; ++n )
    if ( IS_NOT_EQUAL(src_grid->cell_area[n], 0) ) src_grid->cell_frac[n] /= src_grid->cell_area[n];

  for ( n = 0; n < tgt_grid_size; ++n )
    if ( IS_NOT_EQUAL(tgt_grid->cell_area[n], 0) ) tgt_grid->cell_frac[n] /= tgt_grid->cell_area[n];

  /* Perform some error checking on final weights  */

  if ( lcheck )
    {
      for ( n = 0; n < src_grid_size; ++n )
	{
	  if ( src_grid->cell_area[n] < -.01 )
	    cdoPrint("Source grid area error: %d %g", n, src_grid->cell_area[n]);

	  if ( src_centroid_lat[n] < -PIH-.01 || src_centroid_lat[n] > PIH+.01 )
	    cdoPrint("Source grid centroid lat error: %d %g", n, src_centroid_lat[n]);

	  src_centroid_lat[n] = 0;
	  src_centroid_lon[n] = 0;
	}

      for ( n = 0; n < tgt_grid_size; ++n )
	{
	  if ( tgt_grid->cell_area[n] < -.01 )
	    cdoPrint("Target grid area error: %d %g", n, tgt_grid->cell_area[n]);
	  if ( tgt_centroid_lat[n] < -PIH-.01 || tgt_centroid_lat[n] > PIH+.01 )
	    cdoPrint("Target grid centroid lat error: %d %g", n, tgt_centroid_lat[n]);

	  tgt_centroid_lat[n] = 0;
	  tgt_centroid_lon[n] = 0;
	}

      for ( n = 0; n < num_links; ++n )
	{
	  src_grid_add = rv->src_grid_add[n];
	  tgt_grid_add = rv->tgt_grid_add[n];

	  if ( rv->wts[3*n] < -0.01 )
	    cdoPrint("Map weight < 0! grid1idx=%d grid2idx=%d nlink=%d wts=%g",
		     src_grid_add, tgt_grid_add, n, rv->wts[3*n]);

	  if ( rv->norm_opt != NORM_OPT_NONE && rv->wts[3*n] > 1.01 )
	    cdoPrint("Map weight > 1! grid1idx=%d grid2idx=%d nlink=%d wts=%g",
		     src_grid_add, tgt_grid_add, n, rv->wts[3*n]);
	}

      for ( n = 0; n < num_links; ++n )
	{
	  tgt_grid_add = rv->tgt_grid_add[n];
	  tgt_centroid_lat[tgt_grid_add] += rv->wts[3*n];
	}

      /* 2012-01-24 Uwe Schulzweida: changed [tgt_grid_add] to [n] (bug fix) */
      for ( n = 0; n < tgt_grid_size; ++n )
	{
	  if ( rv->norm_opt == NORM_OPT_DESTAREA )
	    norm_factor = tgt_grid->cell_frac[n];
	  else if ( rv->norm_opt == NORM_OPT_FRACAREA )
	    norm_factor = ONE;
	  else if ( rv->norm_opt == NORM_OPT_NONE )
	    norm_factor = tgt_grid->cell_area[n];
	    
	  if ( tgt_centroid_lat[n] > 0 && fabs(tgt_centroid_lat[n] - norm_factor) > .01 )
	    cdoPrint("Error: sum of wts for map1 %d %g %g", n, tgt_centroid_lat[n], norm_factor);
	}
    } // lcheck

  free(src_centroid_lat);
  free(src_centroid_lon);
  free(tgt_centroid_lat);
  free(tgt_centroid_lon);

  if ( ! remap_store_link_fast )
    {
      free(link_add1[0]);
      free(link_add1[1]);
      free(link_add2[0]);
      free(link_add2[1]);
    }

  if ( cdoTimer ) timer_stop(timer_remap_con);

} /* remap_conserv */

/*
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
*/
static
void restrict_boundbox(const double *restrict grid_bound_box, double *restrict bound_box)
{
  if ( bound_box[0] < grid_bound_box[0] && bound_box[1] > grid_bound_box[0] ) bound_box[0] = grid_bound_box[0];
  if ( bound_box[1] > grid_bound_box[1] && bound_box[0] < grid_bound_box[1] ) bound_box[1] = grid_bound_box[1];

  if ( bound_box[2] >= grid_bound_box[3] && (bound_box[3]-2*M_PI) > grid_bound_box[2] ) { bound_box[2] -= 2*M_PI; bound_box[3] -= 2*M_PI; }
  if ( bound_box[3] <= grid_bound_box[2] && (bound_box[2]-2*M_PI) < grid_bound_box[3] ) { bound_box[2] += 2*M_PI; bound_box[3] += 2*M_PI; }
  //  if ( bound_box[2] < grid_bound_box[2] && bound_box[3] > grid_bound_box[2] ) bound_box[2] = grid_bound_box[2];
  //  if ( bound_box[3] > grid_bound_box[3] && bound_box[2] < grid_bound_box[3] ) bound_box[3] = grid_bound_box[3];
}

static
void boundbox_from_corners_reg2d(long grid_add, const int *restrict grid_dims, const double *restrict corner_lon,
				 const double *restrict corner_lat, double *restrict bound_box)
{
  long nx = grid_dims[0];
  long ix, iy;
  double clat1, clat2;

  iy = grid_add/nx;
  ix = grid_add - iy*nx;

  clat1 = corner_lat[iy  ];
  clat2 = corner_lat[iy+1];

  if ( clat2 > clat1 )
    {
      bound_box[0] = clat1;
      bound_box[1] = clat2;
    }
  else
    {
      bound_box[0] = clat2;
      bound_box[1] = clat1;
    }

  bound_box[2] = corner_lon[ix  ];
  bound_box[3] = corner_lon[ix+1];
}

static
void boundbox_from_corners1(long ic, long nc, const double *restrict corner_lon,
			    const double *restrict corner_lat, double *restrict bound_box)
{
  long inc, j;
  double clon, clat;

  inc = ic*nc;

  clat = corner_lat[inc];
  clon = corner_lon[inc];

  bound_box[0] = clat;
  bound_box[1] = clat;
  bound_box[2] = clon;
  bound_box[3] = clon;

  for ( j = 1; j < nc; ++j )
    {
      clat = corner_lat[inc+j];
      clon = corner_lon[inc+j];

      if ( clat < bound_box[0] ) bound_box[0] = clat;
      if ( clat > bound_box[1] ) bound_box[1] = clat;
      if ( clon < bound_box[2] ) bound_box[2] = clon;
      if ( clon > bound_box[3] ) bound_box[3] = clon;
    }

  if ( fabs(bound_box[3] - bound_box[2]) > PI )
    {
      bound_box[2] = 0;
      bound_box[3] = PI2;
    }

  /*
  double dlon = fabs(bound_box[3] - bound_box[2]);

  if ( dlon > PI )
    {
      if ( bound_box[3] > bound_box[2] && (bound_box[3]-PI2) < 0. )
	{
	  double tmp = bound_box[2];
	  bound_box[2] = bound_box[3] - PI2;
	  bound_box[3] = tmp;
	}
      else
	{
	  bound_box[2] = 0;
	  bound_box[3] = PI2;
	}
    }
  */
}

static
void boundbox_from_corners1r(long ic, long nc, const double *restrict corner_lon,
			     const double *restrict corner_lat, restr_t *restrict bound_box)
{
  long inc, j;
  restr_t clon, clat;

  inc = ic*nc;

  clat = RESTR_SCALE(corner_lat[inc]);
  clon = RESTR_SCALE(corner_lon[inc]);

  bound_box[0] = clat;
  bound_box[1] = clat;
  bound_box[2] = clon;
  bound_box[3] = clon;

  for ( j = 1; j < nc; ++j )
    {
      clat = RESTR_SCALE(corner_lat[inc+j]);
      clon = RESTR_SCALE(corner_lon[inc+j]);

      if ( clat < bound_box[0] ) bound_box[0] = clat;
      if ( clat > bound_box[1] ) bound_box[1] = clat;
      if ( clon < bound_box[2] ) bound_box[2] = clon;
      if ( clon > bound_box[3] ) bound_box[3] = clon;
    }

  if ( fabs(bound_box[3] - bound_box[2]) > RESTR_SCALE(PI) )
    {
      bound_box[2] = 0;
      bound_box[3] = RESTR_SCALE(PI2);
    }
  /*
  if ( RESTR_ABS(bound_box[3] - bound_box[2]) > RESTR_SCALE(PI) )
    {
      if ( bound_box[3] > bound_box[2] && (bound_box[3]-RESTR_SCALE(PI2)) < RESTR_SCALE(0.) )
	{
	  restr_t tmp = bound_box[2];
	  bound_box[2] = bound_box[3] - RESTR_SCALE(PI2);
	  bound_box[3] = tmp;
	}
    }
  */
}

//#if defined(HAVE_LIBYAC)
#include "clipping/clipping.h"
#include "clipping/area.h"

static
void cdo_compute_overlap_areas(unsigned N,
			       struct grid_cell *overlap_buffer,
			       struct grid_cell *source_cells,
			       struct grid_cell  target_cell,
			       double *partial_areas)
{
  for ( unsigned n = 0; n < N; n++ )
    {
      overlap_buffer[n].num_corners   = 0;
      overlap_buffer[n].edge_type     = NULL;
      overlap_buffer[n].coordinates_x = NULL;
      overlap_buffer[n].coordinates_y = NULL;
    }

  /* Do the clipping and get the cell for the overlapping area */

  cell_clipping(N, source_cells, target_cell, overlap_buffer);

  /* Get the partial areas for the overlapping regions */

  for ( unsigned n = 0; n < N; n++ )
    {
      partial_areas[n] = huiliers_area(overlap_buffer[n]);

      if ( overlap_buffer[n].coordinates_x != NULL ) free(overlap_buffer[n].coordinates_x);
      if ( overlap_buffer[n].coordinates_y != NULL ) free(overlap_buffer[n].coordinates_y);
      if ( overlap_buffer[n].edge_type     != NULL ) free(overlap_buffer[n].edge_type);
    }

#ifdef VERBOSE
  for ( unsigned n = 0; n < N; n++ )
    printf("overlap area : %lf\n", partial_areas[n]);
#endif
}
//#endif

static
int get_lonlat_circle_index(remapgrid_t *remap_grid)
{
  int lonlat_circle_index = -1;

  if ( remap_grid->num_cell_corners == 4 )
    {
      if ( remap_grid->remap_grid_type == REMAP_GRID_TYPE_REG2D )
	{
	  lonlat_circle_index = 1;
 	}
      else
	{
	  const double* cell_corner_lon = remap_grid->cell_corner_lon;
	  const double* cell_corner_lat = remap_grid->cell_corner_lat;
	  int gridsize = remap_grid->size;
	  int num_i = 0, num_eq0 = 0, num_eq1 = 0;
	  int iadd = gridsize/3-1;

	  if ( iadd == 0 ) iadd++;

	  for ( int i = 0; i < gridsize; i += iadd )
	    {
	      num_i++;

	      if ( IS_EQUAL(cell_corner_lon[i*4+1], cell_corner_lon[i*4+2]) &&
		   IS_EQUAL(cell_corner_lon[i*4+3], cell_corner_lon[i*4+0]) &&
		   IS_EQUAL(cell_corner_lat[i*4+0], cell_corner_lat[i*4+1]) &&
		   IS_EQUAL(cell_corner_lat[i*4+2], cell_corner_lat[i*4+3]) )
		{  
		  num_eq1++;
		}
	      else if ( IS_EQUAL(cell_corner_lon[i*4+0], cell_corner_lon[i*4+1]) &&
			IS_EQUAL(cell_corner_lon[i*4+2], cell_corner_lon[i*4+3]) &&
			IS_EQUAL(cell_corner_lat[i*4+1], cell_corner_lat[i*4+2]) &&
			IS_EQUAL(cell_corner_lat[i*4+3], cell_corner_lat[i*4+0]) )
		{
		  num_eq0++;
		}
	    }

	  if ( num_i == num_eq1 ) lonlat_circle_index = 1;
	  if ( num_i == num_eq0 ) lonlat_circle_index = 0;	      
	}
    }

  //printf("lonlat_circle_index %d\n", lonlat_circle_index);

  return(lonlat_circle_index);
}


void remap_consphere(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /* local variables */

  int    lcheck = TRUE;

  long   ioffset;
  long   src_grid_size;
  long   tgt_grid_size;
  long   src_num_cell_corners;
  long   tgt_num_cell_corners;
  long   src_grid_add;       /* current linear address for source grid cell   */
  long   tgt_grid_add;       /* current linear address for target grid cell   */
  long   n, k;               /* generic counters                        */
  long   nbins, num_links;
  double norm_factor = 0;    /* factor for normalizing wts */
  long   num_wts;
  long   max_srch_cells;     /* num cells in restricted search arrays  */
  long   num_srch_cells;     /* num cells in restricted search arrays  */
  long   srch_corners;       /* num of corners of srch cells           */
  int*   srch_add;           /* global address of cells in srch arrays */
  int    ompthID, i;

  /* Variables necessary if segment manages to hit pole */
  grid_store_t *grid_store = NULL;
  double findex = 0;
  long num_weights = 0;
  long nx = 0, ny = 0;
  int src_remap_grid_type = src_grid->remap_grid_type;
  int tgt_remap_grid_type = tgt_grid->remap_grid_type;
  double src_grid_bound_box[4];
  int lyac = FALSE;

  progressInit();

  nbins = src_grid->num_srch_bins;
  num_wts = rv->num_wts;

  if ( remap_store_link_fast )
    {
      grid_store = malloc(sizeof(grid_store_t));
      grid_store_init(grid_store, tgt_grid->size);
    }

  if ( cdoTimer ) timer_start(timer_remap_con);

  src_grid_size = src_grid->size;
  tgt_grid_size = tgt_grid->size;

  src_num_cell_corners = src_grid->num_cell_corners;
  tgt_num_cell_corners = tgt_grid->num_cell_corners;

  enum edge_type great_circle_type[] = {GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE};
  enum edge_type lonlat_circle_type[] = {LON_CIRCLE, LAT_CIRCLE, LON_CIRCLE, LAT_CIRCLE, LON_CIRCLE, LAT_CIRCLE, LON_CIRCLE, LAT_CIRCLE, LON_CIRCLE};

  enum edge_type *src_edge_type = great_circle_type;
  enum edge_type *tgt_edge_type = great_circle_type;

  if ( src_num_cell_corners == 4 )
    {
      int lonlat_circle_index = get_lonlat_circle_index(src_grid);
      if ( lonlat_circle_index >= 0 ) src_edge_type = &lonlat_circle_type[lonlat_circle_index];
    }

  if ( tgt_num_cell_corners == 4 )
    {
      int lonlat_circle_index = get_lonlat_circle_index(tgt_grid);
      if ( lonlat_circle_index >= 0 ) tgt_edge_type = &lonlat_circle_type[lonlat_circle_index];
    }

  double tgt_area;

  struct grid_cell* tgt_grid_cell;
  struct grid_cell* tgt_grid_cell2[ompNumThreads];  
  for ( i = 0; i < ompNumThreads; ++i )
    {
      tgt_grid_cell2[i] = malloc(sizeof(struct grid_cell));
      tgt_grid_cell2[i]->num_corners   = tgt_num_cell_corners;
      tgt_grid_cell2[i]->edge_type     = tgt_edge_type;
      tgt_grid_cell2[i]->coordinates_x = malloc(tgt_num_cell_corners*sizeof(double));
      tgt_grid_cell2[i]->coordinates_y = malloc(tgt_num_cell_corners*sizeof(double));
    }

  struct grid_cell* src_grid_cells;
  struct grid_cell* overlap_buffer;
  struct grid_cell* src_grid_cells2[ompNumThreads];
  struct grid_cell* overlap_buffer2[ompNumThreads];
  for ( i = 0; i < ompNumThreads; ++i )
    {
      src_grid_cells2[i] = NULL;
      overlap_buffer2[i] = NULL;
    }

  double* partial_areas;
  double* partial_weights;
  double* partial_areas2[ompNumThreads];
  double* partial_weights2[ompNumThreads];
  for ( i = 0; i < ompNumThreads; ++i )
    {
      partial_areas2[i]   = NULL;
      partial_weights2[i] = NULL;
    }

  long max_srch_cells2[ompNumThreads];
  for ( i = 0; i < ompNumThreads; ++i )
    max_srch_cells2[i] = 0;

  int* srch_add2[ompNumThreads];
  for ( i = 0; i < ompNumThreads; ++i )
    srch_add2[i] = malloc(src_grid_size*sizeof(int));

  srch_corners = src_num_cell_corners;

  if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
    {
      nx = src_grid->dims[0];
      ny = src_grid->dims[1];
     
      src_grid_bound_box[0] = src_grid->reg2d_corner_lat[0];
      src_grid_bound_box[1] = src_grid->reg2d_corner_lat[ny];
      if ( src_grid_bound_box[0] > src_grid_bound_box[1] )
	{
	  src_grid_bound_box[0] = src_grid->reg2d_corner_lat[ny];
	  src_grid_bound_box[1] = src_grid->reg2d_corner_lat[0];
	}
      src_grid_bound_box[2] = src_grid->reg2d_corner_lon[0];
      src_grid_bound_box[3] = src_grid->reg2d_corner_lon[nx];
      //printf("src_grid   lon: %g %g lat: %g %g\n", RAD2DEG*src_grid_bound_box[2],RAD2DEG*src_grid_bound_box[3],RAD2DEG*src_grid_bound_box[0],RAD2DEG*src_grid_bound_box[1] );
    }

  if ( cdoTimer ) timer_start(timer_remap_con_l2);

  findex = 0;

  int sum_srch_cells = 0;
  int sum_srch_cells2 = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoTimer, lyac, nbins, num_wts, nx, src_remap_grid_type, tgt_remap_grid_type, src_grid_bound_box,	\
	 src_edge_type, tgt_edge_type, partial_areas2, partial_weights2,  \
         remap_store_link_fast, grid_store, rv, cdoVerbose, max_srch_cells2, \
	 tgt_num_cell_corners, srch_corners, src_grid, tgt_grid, tgt_grid_size, src_grid_size, \
	 overlap_buffer2, src_grid_cells2, srch_add2, tgt_grid_cell2, findex, sum_srch_cells, sum_srch_cells2) \
  private(ompthID, srch_add, tgt_grid_cell, tgt_area, n, k, num_weights, num_srch_cells, max_srch_cells,  \
	  partial_areas, partial_weights, overlap_buffer, src_grid_cells, src_grid_add, tgt_grid_add, ioffset)
#endif
  for ( tgt_grid_add = 0; tgt_grid_add < tgt_grid_size; ++tgt_grid_add )
    {
#if defined(_OPENMP)
      ompthID = omp_get_thread_num();
#else
      ompthID = 0;
#endif

      int lprogress = 1;
      if ( ompthID != 0 ) lprogress = 0;

#if defined(_OPENMP)
#pragma omp atomic
#endif
      findex++;
      if ( lprogress ) progressStatus(0, 1, findex/tgt_grid_size);

      srch_add = srch_add2[ompthID];
      tgt_grid_cell = tgt_grid_cell2[ompthID];

      /* Get search cells */

      if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D && tgt_remap_grid_type == REMAP_GRID_TYPE_REG2D )
	{
	  double tgt_cell_bound_box[4];
	  boundbox_from_corners_reg2d(tgt_grid_add, tgt_grid->dims, tgt_grid->reg2d_corner_lon, tgt_grid->reg2d_corner_lat, tgt_cell_bound_box);
	  restrict_boundbox(src_grid_bound_box, tgt_cell_bound_box);
	  if ( cdoVerbose )
	    printf("bound_box %ld  lon: %g %g lat: %g %g\n",
		   tgt_grid_add, RAD2DEG*tgt_cell_bound_box[2],RAD2DEG*tgt_cell_bound_box[3],RAD2DEG*tgt_cell_bound_box[0],RAD2DEG*tgt_cell_bound_box[1] );
	  num_srch_cells = get_srch_cells_reg2d(src_grid->dims, src_grid->reg2d_corner_lat, src_grid->reg2d_corner_lon,
						tgt_cell_bound_box, srch_add);
	}
      else if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
	{
	  double tgt_cell_bound_box[4];
	  boundbox_from_corners1(tgt_grid_add, tgt_num_cell_corners, tgt_grid->cell_corner_lon, tgt_grid->cell_corner_lat, tgt_cell_bound_box);
	  restrict_boundbox(src_grid_bound_box, tgt_cell_bound_box);
	  if ( cdoVerbose )
	    printf("bound_box %ld  lon: %g %g lat: %g %g\n",
		   tgt_grid_add, RAD2DEG*tgt_cell_bound_box[2],RAD2DEG*tgt_cell_bound_box[3],RAD2DEG*tgt_cell_bound_box[0],RAD2DEG*tgt_cell_bound_box[1] );
	  num_srch_cells = get_srch_cells_reg2d(src_grid->dims, src_grid->reg2d_corner_lat, src_grid->reg2d_corner_lon,
						tgt_cell_bound_box, srch_add);
	}
      else
	{
	  restr_t tgt_cell_bound_box_r[4];
	  boundbox_from_corners1r(tgt_grid_add, tgt_num_cell_corners, tgt_grid->cell_corner_lon, tgt_grid->cell_corner_lat, tgt_cell_bound_box_r);

	  num_srch_cells = get_srch_cells(tgt_grid_add, nbins, tgt_grid->bin_addr, src_grid->bin_addr,
					  tgt_cell_bound_box_r, src_grid->cell_bound_box, src_grid_size, srch_add);
	}

      sum_srch_cells += num_srch_cells;

      if ( cdoVerbose )
	printf("tgt_grid_add %ld  num_srch_cells %ld\n", tgt_grid_add, num_srch_cells);

      if ( num_srch_cells == 0 ) continue;

      if ( tgt_remap_grid_type == REMAP_GRID_TYPE_REG2D )
	{
	  long nx = tgt_grid->dims[0];
	  long ix, iy;

	  iy = tgt_grid_add/nx;
	  ix = tgt_grid_add - iy*nx;

	  tgt_grid_cell->coordinates_x[0] = tgt_grid->reg2d_corner_lon[ix  ];
	  tgt_grid_cell->coordinates_y[0] = tgt_grid->reg2d_corner_lat[iy  ];
	  tgt_grid_cell->coordinates_x[1] = tgt_grid->reg2d_corner_lon[ix+1];
	  tgt_grid_cell->coordinates_y[1] = tgt_grid->reg2d_corner_lat[iy  ];
	  tgt_grid_cell->coordinates_x[2] = tgt_grid->reg2d_corner_lon[ix+1];
	  tgt_grid_cell->coordinates_y[2] = tgt_grid->reg2d_corner_lat[iy+1];
	  tgt_grid_cell->coordinates_x[3] = tgt_grid->reg2d_corner_lon[ix  ];
	  tgt_grid_cell->coordinates_y[3] = tgt_grid->reg2d_corner_lat[iy+1];
	}
      else
	{
	  for ( int ic = 0; ic < tgt_num_cell_corners; ++ic )
	    {
	      tgt_grid_cell->coordinates_x[ic] = tgt_grid->cell_corner_lon[tgt_grid_add*tgt_num_cell_corners+ic];
	      tgt_grid_cell->coordinates_y[ic] = tgt_grid->cell_corner_lat[tgt_grid_add*tgt_num_cell_corners+ic];
	    }
	}

      //printf("target: %ld\n", tgt_grid_add);
      if ( lyac )
	// if ( tgt_grid_add == 682 )
	for ( int n = 0; n < tgt_num_cell_corners; ++n )
	  {
	    printf("  TargetCell.coordinates_x[%d] = %g*rad;\n", n, tgt_grid_cell->coordinates_x[n]/DEG2RAD);
	    printf("  TargetCell.coordinates_y[%d] = %g*rad;\n", n, tgt_grid_cell->coordinates_y[n]/DEG2RAD);
	  }
      
      /* Create search arrays */

      max_srch_cells  = max_srch_cells2[ompthID];
      partial_areas   = partial_areas2[ompthID];
      partial_weights = partial_weights2[ompthID];
      overlap_buffer  = overlap_buffer2[ompthID];
      src_grid_cells  = src_grid_cells2[ompthID];

      if ( num_srch_cells > max_srch_cells )
	{
	  partial_areas   = realloc(partial_areas,   num_srch_cells*sizeof(double));
	  partial_weights = realloc(partial_weights, num_srch_cells*sizeof(double));

	  overlap_buffer = realloc(overlap_buffer, num_srch_cells*sizeof(struct grid_cell));
	  src_grid_cells = realloc(src_grid_cells, num_srch_cells*sizeof(struct grid_cell));

	  for ( n = max_srch_cells; n < num_srch_cells; ++n )
	    {
	      overlap_buffer[n].num_corners   = 0;
	      overlap_buffer[n].edge_type     = NULL;
	      overlap_buffer[n].coordinates_x = NULL;
	      overlap_buffer[n].coordinates_y = NULL;
	    }

	  for ( n = max_srch_cells; n < num_srch_cells; ++n )
	    {
	      src_grid_cells[n].num_corners   = srch_corners;
	      src_grid_cells[n].edge_type     = src_edge_type;
	      src_grid_cells[n].coordinates_x = malloc(srch_corners*sizeof(double));
	      src_grid_cells[n].coordinates_y = malloc(srch_corners*sizeof(double));
	    }

	  max_srch_cells = num_srch_cells;

	  max_srch_cells2[ompthID]  = max_srch_cells;
	  partial_areas2[ompthID]   = partial_areas;
	  partial_weights2[ompthID] = partial_weights;
	  overlap_buffer2[ompthID]  = overlap_buffer;
	  src_grid_cells2[ompthID]  = src_grid_cells;
	}

      // printf("  int ii = 0;\n");
      for ( n = 0; n < num_srch_cells; ++n )
	{
	  src_grid_add = srch_add[n];

	  if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
	    {
	      int ix, iy;

	      iy = src_grid_add/nx;
	      ix = src_grid_add - iy*nx;

	      src_grid_cells[n].coordinates_x[0] = src_grid->reg2d_corner_lon[ix  ];
	      src_grid_cells[n].coordinates_y[0] = src_grid->reg2d_corner_lat[iy  ];
	      src_grid_cells[n].coordinates_x[1] = src_grid->reg2d_corner_lon[ix+1];
	      src_grid_cells[n].coordinates_y[1] = src_grid->reg2d_corner_lat[iy  ];
	      src_grid_cells[n].coordinates_x[2] = src_grid->reg2d_corner_lon[ix+1];
	      src_grid_cells[n].coordinates_y[2] = src_grid->reg2d_corner_lat[iy+1];
	      src_grid_cells[n].coordinates_x[3] = src_grid->reg2d_corner_lon[ix  ];
	      src_grid_cells[n].coordinates_y[3] = src_grid->reg2d_corner_lat[iy+1];
	      /*
	      printf("source1: %ld %ld", num_srch_cells, n);
	      for ( k = 0; k < srch_corners; ++k )
		printf(" %g %g", src_grid_cells[n].coordinates_x[k]/DEG2RAD, src_grid_cells[n].coordinates_y[k]/DEG2RAD);
	      printf("\n");
	      */
	    }
	  else
	    {
	      ioffset = src_grid_add*srch_corners;

	      for ( k = 0; k < srch_corners; ++k )
		{
		  src_grid_cells[n].coordinates_x[k] = src_grid->cell_corner_lon[ioffset+k];
		  src_grid_cells[n].coordinates_y[k] = src_grid->cell_corner_lat[ioffset+k];
		}
	      /*
	      printf("source2: %ld %ld", num_srch_cells, n);
	      for ( k = 0; k < srch_corners; ++k )
		printf(" %g %g", src_grid_cells[n].coordinates_x[k]/DEG2RAD, src_grid_cells[n].coordinates_y[k]/DEG2RAD);
	      printf("\n");
	      */
	    }

	  if ( lyac )
	    //  if ( tgt_grid_add == 682 )
	    {
	      printf("n %d\n", (int)n);
	      for ( k = 0; k < srch_corners; ++k )
		{
		  printf("  SourceCell[ii].coordinates_x[%ld] = %g*rad;\n", k, src_grid_cells[n].coordinates_x[k]/DEG2RAD);
		  printf("  SourceCell[ii].coordinates_y[%ld] = %g*rad;\n", k, src_grid_cells[n].coordinates_y[k]/DEG2RAD);
		}
	      printf("  ii++;\n");
	    }
	}

      cdo_compute_overlap_areas(num_srch_cells, overlap_buffer, src_grid_cells, *tgt_grid_cell, partial_areas);

      tgt_area = huiliers_area(*tgt_grid_cell);
      // tgt_area = cell_area(tgt_grid_cell);

      for ( num_weights = 0, n = 0; n < num_srch_cells; ++n )
	{
	  if ( partial_areas[n] > 0 )
	    {
	      //printf(">>>>   %d %d %g %g\n", (int)tgt_grid_add, srch_add[n], tgt_area, partial_areas[n]);
	      partial_areas[num_weights] = partial_areas[n];
	      srch_add[num_weights] = srch_add[n];
	      num_weights++;
	    }
	}

      sum_srch_cells2 += num_weights;

      for ( n = 0; n < num_weights; ++n )
	partial_weights[n] = partial_areas[n] / tgt_area;

      correct_weights(num_weights, partial_weights);
      //#endif

      for ( n = 0; n < num_weights; ++n )
	{
	  src_grid_add = srch_add[n];

	  if ( cdoVerbose )
	    printf("tgt_grid_add %ld, src_grid_add %ld,  partial_weights[n] %g, tgt_area  %g\n", tgt_grid_add, src_grid_add, partial_weights[n], tgt_area);

	  // src_grid_add = n;
	  if ( partial_weights[n] <= 0. ) src_grid_add = -1;

	  /*
	    Store the appropriate addresses and weights. 
	    Also add contributions to cell areas.
	    The source grid mask is the master mask
	  */
	  if ( src_grid_add != -1 )
	    {
	      if ( src_grid->mask[src_grid_add] )
		{

#if defined(_OPENMP)
#pragma omp critical
#endif
		  {
		    store_link_cnsrv_fast(rv, src_grid_add, tgt_grid_add, num_wts, &partial_weights[n], grid_store);

		    src_grid->cell_frac[src_grid_add] += partial_weights[n];
		  }
		  tgt_grid->cell_frac[tgt_grid_add] += partial_weights[n];
		}
#if defined(_OPENMP)
#pragma omp critical
#endif
	      {
		src_grid->cell_area[src_grid_add] += partial_weights[n];
	      }
	    }

	  tgt_grid->cell_area[tgt_grid_add] += partial_weights[n];
	}
    }

  if ( cdoVerbose )
    {
      printf("sum_srch_cells : %d\n", sum_srch_cells);
      printf("sum_srch_cells2: %d\n", sum_srch_cells2);
    }

  if ( cdoTimer ) timer_stop(timer_remap_con_l2);

  /* Finished with all cells: deallocate search arrays */

  for ( i = 0; i < ompNumThreads; ++i )
    {
      for ( n = 0; n < max_srch_cells2[i]; n++ )
	{
	  free(src_grid_cells2[i][n].coordinates_x);
	  free(src_grid_cells2[i][n].coordinates_y);
	}
      free(src_grid_cells2[i]);
      free(overlap_buffer2[i]);

      free(partial_areas2[i]);
      free(partial_weights2[i]);

      free(tgt_grid_cell2[i]->coordinates_x);
      free(tgt_grid_cell2[i]->coordinates_y);
      free(tgt_grid_cell2[i]);

      free(srch_add2[i]);
    }

  if ( remap_store_link_fast )
    {
      grid_store_delete(grid_store);
      free(grid_store);
    }


  /* Normalize using destination area if requested */

  num_links = rv->num_links;

  if ( rv->norm_opt == NORM_OPT_DESTAREA )
    {
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(num_wts, num_links, rv, tgt_grid)	\
  private(n, tgt_grid_add, norm_factor)
#endif
      for ( n = 0; n < num_links; ++n )
	{
	  tgt_grid_add = rv->tgt_grid_add[n];

          if ( IS_NOT_EQUAL(tgt_grid->cell_area[tgt_grid_add], 0) )
	    norm_factor = ONE/tgt_grid->cell_area[tgt_grid_add];
          else
            norm_factor = ZERO;

	  rv->wts[n*num_wts] *= norm_factor;
	}
    }
  else if ( rv->norm_opt == NORM_OPT_FRACAREA )
    {
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(num_wts, num_links, rv, tgt_grid)	\
  private(n, tgt_grid_add, norm_factor)
#endif
      for ( n = 0; n < num_links; ++n )
	{
	  tgt_grid_add = rv->tgt_grid_add[n];

          if ( IS_NOT_EQUAL(tgt_grid->cell_frac[tgt_grid_add], 0) )
	    norm_factor = ONE/tgt_grid->cell_frac[tgt_grid_add];
          else
            norm_factor = ZERO;

	  rv->wts[n*num_wts] *= norm_factor;
	}
    }
  else if ( rv->norm_opt == NORM_OPT_NONE )
    {
    }

  if ( cdoVerbose )
    cdoPrint("Total number of links = %ld", rv->num_links);

  for ( n = 0; n < src_grid_size; ++n )
    if ( IS_NOT_EQUAL(src_grid->cell_area[n], 0) ) src_grid->cell_frac[n] /= src_grid->cell_area[n];

  for ( n = 0; n < tgt_grid_size; ++n )
    if ( IS_NOT_EQUAL(tgt_grid->cell_area[n], 0) ) tgt_grid->cell_frac[n] /= tgt_grid->cell_area[n];

  /* Perform some error checking on final weights  */

  if ( lcheck )
    {
      for ( n = 0; n < src_grid_size; ++n )
	{
	  if ( src_grid->cell_area[n] < -.01 )
	    cdoPrint("Source grid area error: %d %g", n, src_grid->cell_area[n]);
	}

      for ( n = 0; n < tgt_grid_size; ++n )
	{
	  if ( tgt_grid->cell_area[n] < -.01 )
	    cdoPrint("Target grid area error: %d %g", n, tgt_grid->cell_area[n]);
	}

      for ( n = 0; n < num_links; ++n )
	{
	  src_grid_add = rv->src_grid_add[n];
	  tgt_grid_add = rv->tgt_grid_add[n];

	  if ( rv->wts[n*num_wts] < -0.01 )
	    cdoPrint("Map weight < 0! grid1idx=%d grid2idx=%d nlink=%d wts=%g",
		     src_grid_add, tgt_grid_add, n, rv->wts[n*num_wts]);

	  if ( rv->norm_opt != NORM_OPT_NONE && rv->wts[n*num_wts] > 1.01 )
	    cdoPrint("Map weight > 1! grid1idx=%d grid2idx=%d nlink=%d wts=%g",
		     src_grid_add, tgt_grid_add, n, rv->wts[n*num_wts]);
	}
    } // lcheck

  if ( cdoTimer ) timer_stop(timer_remap_con);

} /* remap_consphere */

/*****************************************************************************/

void remap_stat(int remap_order, remapgrid_t src_grid, remapgrid_t tgt_grid, remapvars_t rv, const double *restrict array1, 
		const double *restrict array2, double missval)
{
  long n, ns, i;
  long idiff, imax, imin, icount;
  int *tgt_count;
  double minval, maxval, sum;
	  
  if ( remap_order == 2 )
    cdoPrint("Second order mapping from grid1 to grid2:");
  else
    cdoPrint("First order mapping from grid1 to grid2:");
  cdoPrint("----------------------------------------");

  ns = 0;
  sum = 0;
  minval =  DBL_MAX;
  maxval = -DBL_MAX;
  for ( n = 0; n < src_grid.size; ++n )
    {
      if ( !DBL_IS_EQUAL(array1[n], missval) )
	{
	  if ( array1[n] < minval ) minval = array1[n];
	  if ( array1[n] > maxval ) maxval = array1[n];
	  sum += array1[n];
	  ns++;
	}
    }
  if ( ns > 0 ) sum /= ns;
  cdoPrint("Grid1 min,mean,max: %g %g %g", minval, sum, maxval);

  ns = 0;
  sum = 0;
  minval =  DBL_MAX;
  maxval = -DBL_MAX;
  for ( n = 0; n < tgt_grid.size; ++n )
    {
      if ( !DBL_IS_EQUAL(array2[n], missval) )
	{
	  if ( array2[n] < minval ) minval = array2[n];
	  if ( array2[n] > maxval ) maxval = array2[n];
	  sum += array2[n];
	  ns++;
	}
    }
  if ( ns > 0 ) sum /= ns;
  cdoPrint("Grid2 min,mean,max: %g %g %g", minval, sum, maxval);

  /* Conservation Test */

  if ( src_grid.cell_area )
    {
      cdoPrint("Conservation:");
      sum = 0;
      for ( n = 0; n < src_grid.size; ++n )
	if ( !DBL_IS_EQUAL(array1[n], missval) )
	  sum += array1[n]*src_grid.cell_area[n]*src_grid.cell_frac[n];
      cdoPrint("Grid1 Integral = %g", sum);

      sum = 0;
      for ( n = 0; n < tgt_grid.size; ++n )
	if ( !DBL_IS_EQUAL(array2[n], missval) )
	  sum += array2[n]*tgt_grid.cell_area[n]*tgt_grid.cell_frac[n];
      cdoPrint("Grid2 Integral = %g", sum);
      /*
      for ( n = 0; n < src_grid.size; n++ )
       fprintf(stderr, "1 %d %g %g %g\n", n, array1[n], src_grid.cell_area[n], src_grid.cell_frac[n]);
      for ( n = 0; n < tgt_grid.size; n++ )
	fprintf(stderr, "2 %d %g %g %g\n", n, array2[n], tgt_grid.cell_area[n], tgt_grid.cell_frac[n]);
      */
    }

  cdoPrint("number of sparse matrix entries %d", rv.num_links);
  cdoPrint("total number of dest cells %d", tgt_grid.size);

  tgt_count = malloc(tgt_grid.size*sizeof(int));

  for ( n = 0; n < tgt_grid.size; ++n ) tgt_count[n] = 0;

#if defined(SX)
#pragma vdir nodep
#endif
  for ( n = 0; n < rv.num_links; ++n ) tgt_count[rv.tgt_grid_add[n]]++;

  imin = INT_MAX;
  imax = INT_MIN;
  for ( n = 0; n < tgt_grid.size; ++n )
    {
      if ( tgt_count[n] > 0 )
	if ( tgt_count[n] < imin ) imin = tgt_count[n];
      if ( tgt_count[n] > imax ) imax = tgt_count[n];
    }

  idiff =  (imax - imin)/10 + 1;
  icount = 0;
  for ( i = 0; i < tgt_grid.size; ++i )
    if ( tgt_count[i] > 0 ) icount++;

  cdoPrint("number of cells participating in remap %d", icount);

  if ( icount )
    {
      cdoPrint("min no of entries/row = %d", imin);
      cdoPrint("max no of entries/row = %d", imax);

      imax = imin + idiff;
      for ( n = 0; n < 10; ++n )
	{
	  icount = 0;
	  for ( i = 0; i < tgt_grid.size; ++i )
	    if ( tgt_count[i] >= imin && tgt_count[i] < imax ) icount++;

	  if ( icount )
	    cdoPrint("num of rows with entries between %d - %d  %d", imin, imax-1, icount);

	  imin = imin + idiff;
	  imax = imax + idiff;
	}
    }

  free(tgt_count);

  if ( rv.sort_add )
    cdoPrint("Sparse matrix entries are explicitly sorted.");

} /* remap_stat */

/*****************************************************************************/

void remap_gradients(remapgrid_t grid, const double *restrict array, double *restrict grad_lat,
		     double *restrict grad_lon, double *restrict grad_latlon)
{
  long n, nx, ny, grid_size;
  long i, j, ip1, im1, jp1, jm1, in, is, ie, iw, ine, inw, ise, isw;
  double delew, delns;
  double grad_lat_zero, grad_lon_zero;

  if ( grid.rank != 2 )
    cdoAbort("Internal problem (remap_gradients), grid rank = %d!", grid.rank);

  grid_size = grid.size;
  nx = grid.dims[0];
  ny = grid.dims[1];

#if defined(_OPENMP)
#pragma omp parallel for default(none)        \
  shared(grid_size, grad_lat, grad_lon, grad_latlon, grid, nx, ny, array) \
  private(n, i, j, ip1, im1, jp1, jm1, in, is, ie, iw, ine, inw, ise, isw, delew, delns, grad_lat_zero, grad_lon_zero)
#endif
  for ( n = 0; n < grid_size; ++n )
    {
      grad_lat[n] = ZERO;
      grad_lon[n] = ZERO;
      grad_latlon[n] = ZERO;

      if ( grid.mask[n] )
	{
	  delew = HALF;
	  delns = HALF;

	  j = n/nx + 1;
	  i = n - (j-1)*nx + 1;

	  ip1 = i + 1;
	  im1 = i - 1;
	  jp1 = j + 1;
	  jm1 = j - 1;

	  if ( ip1 > nx ) ip1 = ip1 - nx;
	  if ( im1 <  1 ) im1 = nx;
	  if ( jp1 > ny )
	    {
              jp1 = j;
              delns = ONE;
            }
	  if ( jm1 < 1 )
	    {
              jm1 = j;
              delns = ONE;
            }

	  in  = (jp1-1)*nx + i - 1;
	  is  = (jm1-1)*nx + i - 1;
	  ie  = (j  -1)*nx + ip1 - 1;
	  iw  = (j  -1)*nx + im1 - 1;

	  ine = (jp1-1)*nx + ip1 - 1;
	  inw = (jp1-1)*nx + im1 - 1;
	  ise = (jm1-1)*nx + ip1 - 1;
	  isw = (jm1-1)*nx + im1 - 1;

	  /* Compute i-gradient */

	  if ( ! grid.mask[ie] )
	    {
              ie = n;
              delew = ONE;
            }
	  if ( ! grid.mask[iw] )
	    {
              iw = n;
              delew = ONE;
            }
 
	  grad_lat[n] = delew*(array[ie] - array[iw]);

	  /* Compute j-gradient */

	  if ( ! grid.mask[in] )
	    {
              in = n;
              delns = ONE;
            }
	  if ( ! grid.mask[is] )
	    {
              is = n;
              delns = ONE;
            }
 
	  grad_lon[n] = delns*(array[in] - array[is]);

	  /* Compute ij-gradient */

	  delew = HALF;
	  if ( jp1 == j || jm1 == j )
	    delns = ONE;
	  else 
	    delns = HALF;

	  if ( ! grid.mask[ine] )
	    {
              if ( in != n )
		{
		  ine = in;
		  delew = ONE;
		}
              else if ( ie != n )
		{
		  ine = ie;
		  inw = iw;
		  if ( inw == n ) delew = ONE;
		  delns = ONE;
		}
              else
		{
		  ine = n;
		  inw = iw;
		  delew = ONE;
		  delns = ONE;
		}
	    }

	  if ( ! grid.mask[inw] )
	    {
              if ( in != n )
		{
		  inw = in;
		  delew = ONE;
		}
              else if ( iw != n )
		{
		  inw = iw;
		  ine = ie;
		  if ( ie == n ) delew = ONE;
		  delns = ONE;
		}
              else
		{
		  inw = n;
		  ine = ie;
		  delew = ONE;
		  delns = ONE;
		}
	    }

	  grad_lat_zero = delew*(array[ine] - array[inw]);

	  if ( ! grid.mask[ise] )
	    {
              if ( is != n )
		{
		  ise = is;
		  delew = ONE;
		}
              else if ( ie != n )
		{
		  ise = ie;
		  isw = iw;
		  if ( isw == n ) delew = ONE;
		  delns = ONE;
		}
              else
		{
		  ise = n;
		  isw = iw;
		  delew = ONE;
		  delns = ONE;
		}
	    }

	  if ( ! grid.mask[isw] )
	    {
              if ( is != n )
		{
		  isw = is;
		  delew = ONE;
		}
              else if ( iw != n )
		{
		  isw = iw;
		  ise = ie;
		  if ( ie == n ) delew = ONE;
		  delns = ONE;
		}
              else
		{
		  isw = n;
		  ise = ie;
		  delew = ONE;
		  delns = ONE;
		}
	    }

	  grad_lon_zero = delew*(array[ise] - array[isw]);

	  grad_latlon[n] = delns*(grad_lat_zero - grad_lon_zero);
	}
    }
} /* remap_gradients */

/*****************************************************************************/

void reorder_links(remapvars_t *rv)
{
  long j, nval = 0, num_blks = 0;
  long lastval;
  long nlinks;
  long max_links = 0;
  long n;
  long num_links;

  num_links = rv->num_links;

  printf("reorder_links\n");
  printf("  num_links %ld\n", num_links);
  rv->links.option = TRUE;

  lastval = -1;
  for ( n = 0; n < num_links; n++ )
    {
      if ( rv->tgt_grid_add[n] == lastval ) nval++;
      else
	{
	  if ( nval > num_blks ) num_blks = nval;
	  nval = 1;
	  max_links++;
	  lastval = rv->tgt_grid_add[n];
	}
    }

  if ( num_blks )
    {
      rv->links.max_links = max_links;
      rv->links.num_blks  = num_blks;

      printf("num_links %ld  max_links %ld  num_blks %ld\n", rv->num_links, max_links, num_blks);

      rv->links.num_links = malloc(num_blks*sizeof(int));
      rv->links.dst_add   = malloc(num_blks*sizeof(int *));
      rv->links.src_add   = malloc(num_blks*sizeof(int *));
      rv->links.w_index   = malloc(num_blks*sizeof(int *));
    }

  for ( j = 0; j < num_blks; j++ )
    {
      rv->links.dst_add[j] = malloc(max_links*sizeof(int));
      rv->links.src_add[j] = malloc(max_links*sizeof(int));
      rv->links.w_index[j] = malloc(max_links*sizeof(int));
    }

  for ( j = 0; j < num_blks; j++ )
    {
      nval = 0;
      lastval = -1;
      nlinks = 0;

      for ( n = 0; n < num_links; n++ )
	{
	  if ( rv->tgt_grid_add[n] == lastval ) nval++;
	  else
	    {
	      nval = 1;
	      lastval = rv->tgt_grid_add[n];
	    }
	  
	  if ( nval == j+1 )
	    {
	      rv->links.dst_add[j][nlinks] = rv->tgt_grid_add[n];
	      rv->links.src_add[j][nlinks] = rv->src_grid_add[n];
	      rv->links.w_index[j][nlinks] = n;
	      nlinks++;
	    }
	}

      rv->links.num_links[j] = nlinks;
      printf("loop %ld  nlinks %ld\n", j+1, nlinks);
    }
}
