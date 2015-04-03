#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      INTERPOLATION USING A DISTANCE-WEIGHTED AVERAGE                    */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

static
void get_restrict_add(remapgrid_t *src_grid, double plat, const int *restrict src_bin_add, int *minadd, int *maxadd)
{
  int n, n2;
  int min_add = 0, max_add = 0, nm1, np1;
  int nbins;
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

static
void nbr_store_distance(int nadd, double distance, int num_neighbors, int *restrict nbr_add, double *restrict nbr_dist)
{
  if ( num_neighbors == 1 )
    {
      if ( distance < nbr_dist[0] || (distance <= nbr_dist[0] && nadd < nbr_add[0]) )
	{
	  nbr_add[0]  = nadd;
	  nbr_dist[0] = distance;
	}
    }
  else
    {
      int n, nchk;
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
}

static
void nbr_check_distance(int num_neighbors, const int *restrict nbr_add, double *restrict nbr_dist)
{
  int nchk;
  double distance;

  /* Uwe Schulzweida: if distance is zero, set to small number */
  for ( nchk = 0; nchk < num_neighbors; ++nchk )
    {
      if ( nbr_add[nchk] >= 0 )
	{
	  distance = nbr_dist[nchk];
	  if ( IS_EQUAL(distance, 0.) ) distance = TINY;
	  nbr_dist[nchk] = distance;
	}
    }
}


static
double get_search_radius(void)
{
  double search_radius;
  extern double remap_search_radius;

  search_radius = remap_search_radius;

  if ( search_radius <    0. ) search_radius = 0.;
  if ( search_radius >  180. ) search_radius = 180.;

  search_radius = cos(search_radius*DEG2RAD);

  return (search_radius);
}

/*
   This routine finds the closest num_neighbor points to a search 
   point and computes a distance to each of the neighbors.
*/
static
void grid_search_nbr_reg2d(int num_neighbors, remapgrid_t *src_grid, int *restrict nbr_add, double *restrict nbr_dist, 
			   double plat, double plon, const int *restrict src_grid_dims,
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
  int n, nadd;
  int nx, nxm, ny;
  long ii, jj;
  int i, j, ix;
  int src_add[25];
  int num_add = 0;
  double distance;   //  Angular distance
  double search_radius = get_search_radius();
  double coslat_dst = cos(plat);  // cos(lat)  of the search point
  double coslon_dst = cos(plon);  // cos(lon)  of the search point
  double sinlat_dst = sin(plat);  // sin(lat)  of the search point
  double sinlon_dst = sin(plon);  // sin(lon)  of the search point
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
      int ix, iy;

      for ( int na = 0; na < num_add; ++na )
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
	  if ( distance >  1. ) distance =  1.;

	  if ( distance >= search_radius )
	    {
	      distance = acos(distance);

	      /* Store the address and distance if this is one of the smallest four so far */
	      nbr_store_distance(nadd, distance, num_neighbors, nbr_add, nbr_dist);
	    }
	}

      nbr_check_distance(num_neighbors, nbr_add, nbr_dist);
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
  int n, nadd;
  int min_add, max_add;
  double search_radius = get_search_radius();
  double coslat_dst = cos(plat);  // cos(lat)  of the search point
  double coslon_dst = cos(plon);  // cos(lon)  of the search point
  double sinlat_dst = sin(plat);  // sin(lat)  of the search point
  double sinlon_dst = sin(plon);  // sin(lon)  of the search point
  /* Loop over source grid and find nearest neighbors                         */
  /* restrict the search using search bins expand the bins to catch neighbors */

  get_restrict_add(src_grid, plat, src_bin_add, &min_add, &max_add);

  /* Initialize distance and address arrays */
  for ( n = 0; n < num_neighbors; ++n )
    {
      nbr_add[n]  = -1;
      nbr_dist[n] = BIGNUM;
    }

  int i, ndist = max_add - min_add + 1;

  if ( ndist <= 0 ) return;

  double distance;     /* Angular distance */
  double *dist = (double*) malloc(ndist*sizeof(double));
  int    *adds = (int*) malloc(ndist*sizeof(int));

  int j = 0;
  for ( i = 0; i < ndist; ++i )
    {
      nadd = min_add+i;
      /* Find distance to this point */
      distance =  sinlat_dst*sinlat[nadd] + coslat_dst*coslat[nadd]*
	         (coslon_dst*coslon[nadd] + sinlon_dst*sinlon[nadd]);
      /* 2008-07-30 Uwe Schulzweida: check that distance is inside the range of -1 to 1,
                                     otherwise the result of acos(distance) is NaN */
      if ( distance >  1. ) distance =  1.;

      if ( distance >= search_radius )
	{
	  dist[j] = distance;
	  adds[j] = nadd;
	  j++;
	}
    }
  ndist = j;

#if defined(_OPENMP) && _OPENMP >= OPENMP4
#pragma omp simd
#endif
  for ( j = 0; j < ndist; ++j )
    dist[j] = acos(dist[j]);

  for ( j = 0; j < ndist; ++j )
    nbr_store_distance(adds[j], dist[j], num_neighbors, nbr_add, nbr_dist);

  free(adds);
  free(dist);

  nbr_check_distance(num_neighbors, nbr_add, nbr_dist);

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

typedef struct
{
  int add;
  double wgts;
}
addwgts_t;

static
int cmpwgts(const void *s1, const void *s2)
{
  int cmp = 0;
  const addwgts_t *c1 = s1;
  const addwgts_t *c2 = s2;

  if      ( c1->add < c2->add ) cmp = -1;
  else if ( c1->add > c2->add ) cmp =  1;

  return (cmp);
}

static
void sort_dist_adds(int nadds, int src_add[], double wgts[])
{
  int n;
  addwgts_t addwgts[nadds];

  if ( nadds <= 1 ) return;

  for ( n = 1; n < nadds; ++n )
    if ( src_add[n] < src_add[n-1] ) break;
  if ( n == nadds ) return;

  for ( n = 0; n < nadds; ++n )
    {
      addwgts[n].add  = src_add[n];
      addwgts[n].wgts = wgts[n];
    }

  qsort(addwgts, nadds, sizeof(addwgts_t), cmpwgts);

  for ( n = 0; n < nadds; ++n )
    {
      src_add[n] = addwgts[n].add;
      wgts[n]    = addwgts[n].wgts;
    }  
}

/*
  -----------------------------------------------------------------------------------------

   This routine computes the inverse-distance weights for a nearest-neighbor interpolation.

  -----------------------------------------------------------------------------------------
*/
void scrip_remap_weights_distwgt(int num_neighbors, remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /*  Local variables */

  long src_grid_size;
  long tgt_grid_size;
  long n, nadds;
  long dst_add;                   /* destination address                         */
  double dist_tot;                /* sum of neighbor distances (for normalizing) */
  double *coslat, *sinlat;        /* cosine, sine of grid lats (for distance)    */
  double *coslon, *sinlon;        /* cosine, sine of grid lons (for distance)    */
  double plat, plon;              /* lat/lon coords of destination point         */
  double findex = 0;
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  /* Compute mappings from source to target grid */

  src_grid_size = src_grid->size;
  tgt_grid_size = tgt_grid->size;

  /* Compute cos, sin of lat/lon on source grid for distance calculations */

  if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
    {
      long nx = src_grid->dims[0];
      long ny = src_grid->dims[1];

      coslat = (double*) malloc(ny*sizeof(double));
      coslon = (double*) malloc(nx*sizeof(double));
      sinlat = (double*) malloc(ny*sizeof(double));
      sinlon = (double*) malloc(nx*sizeof(double));

      double *center_lon = src_grid->reg2d_center_lon;
      double *center_lat = src_grid->reg2d_center_lat;

      for ( n = 0; n < nx; ++n )
	{
	  double rlon = center_lon[n];
	  if ( rlon > PI2  ) rlon -= PI2;
	  if ( rlon < ZERO ) rlon += PI2;
	  coslon[n] = cos(rlon);
	  sinlon[n] = sin(rlon);
	}
      for ( n = 0; n < ny; ++n )
	{
	  coslat[n] = cos(center_lat[n]);
	  sinlat[n] = sin(center_lat[n]);
	}
    }
  else
    {
      coslat = (double*) malloc(src_grid_size*sizeof(double));
      coslon = (double*) malloc(src_grid_size*sizeof(double));
      sinlat = (double*) malloc(src_grid_size*sizeof(double));
      sinlon = (double*) malloc(src_grid_size*sizeof(double));

      double *center_lon = src_grid->cell_center_lon;
      double *center_lat = src_grid->cell_center_lat;

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(center_lon, center_lat, src_grid_size, coslat, coslon, sinlat, sinlon)
#endif
      for ( n = 0; n < src_grid_size; ++n )
	{
	  coslon[n] = cos(center_lon[n]);
	  sinlon[n] = sin(center_lon[n]);
	  coslat[n] = cos(center_lat[n]);
	  sinlat[n] = sin(center_lat[n]);
	}
    }

  /* Loop over destination grid  */
#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoTimer, num_neighbors, remap_grid_type, src_grid, tgt_grid, rv, tgt_grid_size, coslat, coslon, sinlat, sinlon, findex) \
  private(dst_add, n, nadds, dist_tot, plat, plon) \
  schedule(dynamic,1)
#endif
  for ( dst_add = 0; dst_add < tgt_grid_size; ++dst_add )
    {
      int nbr_mask[num_neighbors];    /* mask at nearest neighbors                   */
      int nbr_add[num_neighbors];     /* source address at nearest neighbors         */
      double nbr_dist[num_neighbors]; /* angular distance four nearest neighbors     */
      int lprogress = 1;
      if ( cdo_omp_get_thread_num() != 0 ) lprogress = 0;

#if defined(_OPENMP)
#pragma omp atomic
#endif
      findex++;
      if ( lprogress ) progressStatus(0, 1, findex/tgt_grid_size);

      if ( ! tgt_grid->mask[dst_add] ) continue;
	
      plat = tgt_grid->cell_center_lat[dst_add];
      plon = tgt_grid->cell_center_lon[dst_add];

      /* Find nearest grid points on source grid and distances to each point */
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	grid_search_nbr_reg2d(num_neighbors, src_grid, nbr_add, nbr_dist, 
			      plat, plon, src_grid->dims,
			      sinlat, coslat, sinlon, coslon,
			      src_grid->reg2d_center_lat, src_grid->reg2d_center_lon);
      else
	grid_search_nbr(num_neighbors, src_grid, nbr_add, nbr_dist, 
			plat, plon, src_grid->bin_addr,
			sinlat, coslat, sinlon, coslon);

      /* Compute weights based on inverse distance if mask is false, eliminate those points */

      dist_tot = 0.;
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

      nadds = 0;
      for ( n = 0; n < num_neighbors; ++n )
	{
          if ( nbr_mask[n] )
	    {
	      nbr_dist[nadds] = nbr_dist[n]/dist_tot;
	      nbr_add[nadds]  = nbr_add[n];
	      nadds++;

	      tgt_grid->cell_frac[dst_add] = ONE;
	    }
	}

      sort_dist_adds(nadds, nbr_add, nbr_dist);

      for ( n = 0; n < nadds; ++n )
	{
#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_nbr(rv, nbr_add[n], dst_add, nbr_dist[n]);
	}

    } /* for ( dst_add = 0; dst_add < tgt_grid_size; ++dst_add ) */

  free(coslat);
  free(coslon);
  free(sinlat);
  free(sinlon);

}  /* scrip_remap_weights_distwgt */
