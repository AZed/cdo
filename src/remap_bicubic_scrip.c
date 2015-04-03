#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"


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
void store_link_bicub(remapvars_t *rv, int dst_add, int src_add[4], double weights[4][4])
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

  for ( n = 0; n < 4; ++n )
    {
      rv->src_grid_add[nlink+n] = src_add[n];
      rv->tgt_grid_add[nlink+n] = dst_add;
      for ( k = 0; k < 4; ++k )
	rv->wts[4*(nlink+n)+k] = weights[k][n];
    }

} /* store_link_bicub */

static
void set_bicubic_weights(double iw, double jw, double wgts[4][4])
{
  wgts[0][0] = (1.-jw*jw*(3.-2.*jw))  * (1.-iw*iw*(3.-2.*iw));
  wgts[0][1] = (1.-jw*jw*(3.-2.*jw))  *     iw*iw*(3.-2.*iw);
  wgts[0][2] =     jw*jw*(3.-2.*jw)   *     iw*iw*(3.-2.*iw);
  wgts[0][3] =     jw*jw*(3.-2.*jw)   * (1.-iw*iw*(3.-2.*iw));
  wgts[1][0] = (1.-jw*jw*(3.-2.*jw))  *     iw*(iw-1.)*(iw-1.);
  wgts[1][1] = (1.-jw*jw*(3.-2.*jw))  *     iw*iw*(iw-1.);
  wgts[1][2] =     jw*jw*(3.-2.*jw)   *     iw*iw*(iw-1.);
  wgts[1][3] =     jw*jw*(3.-2.*jw)   *     iw*(iw-1.)*(iw-1.);
  wgts[2][0] =     jw*(jw-1.)*(jw-1.) * (1.-iw*iw*(3.-2.*iw));
  wgts[2][1] =     jw*(jw-1.)*(jw-1.) *     iw*iw*(3.-2.*iw);
  wgts[2][2] =     jw*jw*(jw-1.)      *     iw*iw*(3.-2.*iw);
  wgts[2][3] =     jw*jw*(jw-1.)      * (1.-iw*iw*(3.-2.*iw));
  wgts[3][0] =     iw*(iw-1.)*(iw-1.) *     jw*(jw-1.)*(jw-1.);
  wgts[3][1] =     iw*iw*(iw-1.)      *     jw*(jw-1.)*(jw-1.);
  wgts[3][2] =     iw*iw*(iw-1.)      *     jw*jw*(jw-1.);
  wgts[3][3] =     iw*(iw-1.)*(iw-1.) *     jw*jw*(jw-1.);
}

int num_src_points(const int* restrict mask, const int src_add[4], double src_lats[4]);

static
void renormalize_weights(const double src_lats[4], double wgts[4][4])
{
  int n;
  double sum_wgts = 0.0; /* sum of weights for normalization */
  /* 2012-05-08 Uwe Schulzweida: using absolute value of src_lats (bug fix) */
  for ( n = 0; n < 4; ++n ) sum_wgts  += fabs(src_lats[n]);
  for ( n = 0; n < 4; ++n ) wgts[0][n] = fabs(src_lats[n])/sum_wgts;
  for ( n = 0; n < 4; ++n ) wgts[1][n] = 0.;
  for ( n = 0; n < 4; ++n ) wgts[2][n] = 0.;
  for ( n = 0; n < 4; ++n ) wgts[3][n] = 0.;
}

static
void bicubic_warning(void)
{
  static int lwarn = TRUE;

  if ( cdoVerbose || lwarn )
    {
      lwarn = FALSE;
      // cdoWarning("Iteration for iw,jw exceed max iteration count of %d!", remap_max_iter);
      cdoWarning("Bicubic interpolation failed for some grid points - used a distance-weighted average instead!");
    }
}

typedef struct
{
  int add;
  double wgts[4];
}
addwgts_t;

static
int cmpwgts(const void *s1, const void *s2)
{
  int cmp = 0;
  const addwgts_t* c1 = (const addwgts_t*) s1;
  const addwgts_t* c2 = (const addwgts_t*) s2;

  if      ( c1->add < c2->add ) cmp = -1;
  else if ( c1->add > c2->add ) cmp =  1;

  return (cmp);
}

static
void sort_bicubic_adds(int src_add[4], double wgts[4][4])
{
  int n;
  addwgts_t addwgts[4];

  for ( n = 1; n < 4; ++n )
    if ( src_add[n] < src_add[n-1] ) break;
  if ( n == 4 ) return;

  for ( n = 0; n < 4; ++n )
    {
      addwgts[n].add     = src_add[n];
      addwgts[n].wgts[0] = wgts[0][n];
      addwgts[n].wgts[1] = wgts[1][n];
      addwgts[n].wgts[2] = wgts[2][n];
      addwgts[n].wgts[3] = wgts[3][n];
    }

  qsort(addwgts, 4, sizeof(addwgts_t), cmpwgts);

  for ( n = 0; n < 4; ++n )
    {
      src_add[n] = addwgts[n].add;
      wgts[0][n] = addwgts[n].wgts[0];
      wgts[1][n] = addwgts[n].wgts[1];
      wgts[2][n] = addwgts[n].wgts[2];
      wgts[3][n] = addwgts[n].wgts[3];
    }  
}

static
void bicubic_remap(double* restrict tgt_point, const double* restrict src_array, double wgts[4][4], const int src_add[4],
		   const double* restrict grad1, const double* restrict grad2, const double* restrict grad3)
{
  *tgt_point = 0.;
  for ( int n = 0; n < 4; ++n )
    *tgt_point += src_array[src_add[n]]*wgts[0][n] +
                      grad1[src_add[n]]*wgts[1][n] +
                      grad2[src_add[n]]*wgts[2][n] +
                      grad3[src_add[n]]*wgts[3][n];  
}

/*
  -----------------------------------------------------------------------

  This routine computes the weights for a bicubic interpolation.

  -----------------------------------------------------------------------
*/
void scrip_remap_weights_bicubic(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /*   Local variables */
  int  search_result;
  long tgt_grid_size;
  long dst_add;        /*  destination addresss                   */
  int src_add[4];      /*  address for the four source points     */
  double src_lats[4];  /*  latitudes  of four bilinear corners    */
  double src_lons[4];  /*  longitudes of four bilinear corners    */
  double wgts[4][4];   /*  bicubic weights for four corners       */
  double plat, plon;   /*  lat/lon coords of destination point    */
  double findex = 0;
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  tgt_grid_size = tgt_grid->size;

  /* Compute mappings from source to target grid */

  if ( src_grid->rank != 2 )
    cdoAbort("Can not do bicubic interpolation when source grid rank != 2"); 

  /* Loop over destination grid */

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoTimer, cdoVerbose, remap_grid_type, tgt_grid_size, src_grid, tgt_grid, rv, findex) \
  private(dst_add, src_add, src_lats, src_lons, wgts, plat, plon, search_result) \
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
	  for ( int n = 0; n < 4; ++n )
	    if ( ! src_grid->mask[src_add[n]] ) search_result = 0;
	}

      /* If point found, find local iw,jw coordinates for weights  */
      if ( search_result > 0 )
	{
	  double iw, jw;  /*  current guess for bilinear coordinate  */

          tgt_grid->cell_frac[dst_add] = 1.;

          if ( find_ij_weights(plon, plat, src_lats, src_lons, &iw, &jw) )
	    {
	      /* Successfully found iw,jw - compute weights */
	      set_bicubic_weights(iw, jw, wgts);

	      sort_bicubic_adds(src_add, wgts);

#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_bicub(rv, dst_add, src_add, wgts);
	    }
          else
	    {
	      bicubic_warning();

	      search_result = -1;
	    }
	}
	  
      /*
	Search for bicubic failed - use a distance-weighted average instead (this is typically near the pole)
	Distance was stored in src_lats!
      */
      if ( search_result < 0 )
	{
	  if ( num_src_points(src_grid->mask, src_add, src_lats) > 0 )
	    {
	      renormalize_weights(src_lats, wgts);

	      tgt_grid->cell_frac[dst_add] = 1.;

	      sort_bicubic_adds(src_add, wgts);

#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_bicub(rv, dst_add, src_add, wgts);
	    }
        }
    } /* grid_loop1 */

} /* scrip_remap_weights_bicubic */

/*
  -----------------------------------------------------------------------

  This routine computes ans apply the weights for a bicubic interpolation.

  -----------------------------------------------------------------------
*/
void scrip_remap_bicubic(remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double* restrict src_array, double* restrict tgt_array, double missval)
{
  /*   Local variables */
  int  search_result;
  long tgt_grid_size;
  long dst_add;        /*  destination addresss                 */
  int src_add[4];      /*  address for the four source points   */
  double src_lats[4];  /*  latitudes  of four bilinear corners  */
  double src_lons[4];  /*  longitudes of four bilinear corners  */
  double wgts[4][4];   /*  bicubic weights for four corners     */
  double plat, plon;   /*  lat/lon coords of destination point  */
  double findex = 0;
  double *grad1_lat, *grad1_lon, *grad1_latlon;
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  tgt_grid_size = tgt_grid->size;

  /* Compute mappings from source to target grid */

  if ( src_grid->rank != 2 )
    cdoAbort("Can not do bicubic interpolation when source grid rank != 2"); 

  grad1_lat    = (double*) malloc(src_grid->size*sizeof(double));
  grad1_lon    = (double*) malloc(src_grid->size*sizeof(double));
  grad1_latlon = (double*) malloc(src_grid->size*sizeof(double));

  remap_gradients(*src_grid, src_array, grad1_lat, grad1_lon, grad1_latlon);

  /* Loop over destination grid */

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoTimer, cdoVerbose, remap_grid_type, tgt_grid_size, src_grid, tgt_grid, src_array, tgt_array, missval, grad1_lat, grad1_lon, grad1_latlon, findex) \
  private(dst_add, src_add, src_lats, src_lons, wgts, plat, plon, search_result) \
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

      tgt_array[dst_add] = missval;

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
	  for ( int n = 0; n < 4; ++n )
	    if ( ! src_grid->mask[src_add[n]] ) search_result = 0;
	}

      /* If point found, find local iw,jw coordinates for weights  */
      if ( search_result > 0 )
	{
	  double iw, jw;  /*  current guess for bilinear coordinate  */

          tgt_grid->cell_frac[dst_add] = 1.;

          if ( find_ij_weights(plon, plat, src_lats, src_lons, &iw, &jw) )
	    {
	      /* Successfully found iw,jw - compute weights */
	      set_bicubic_weights(iw, jw, wgts);

	      sort_bicubic_adds(src_add, wgts);

	      bicubic_remap(&tgt_array[dst_add], src_array, wgts, src_add, grad1_lat, grad1_lon, grad1_latlon);
	    }
          else
	    {
	      bicubic_warning();

	      search_result = -1;
	    }
	}
	  
      /*
	Search for bicubic failed - use a distance-weighted average instead (this is typically near the pole)
	Distance was stored in src_lats!
      */
      if ( search_result < 0 )
	{
	  if ( num_src_points(src_grid->mask, src_add, src_lats) > 0 )
	    {
	      renormalize_weights(src_lats, wgts);

	      tgt_grid->cell_frac[dst_add] = 1.;

	      sort_bicubic_adds(src_add, wgts);

	      bicubic_remap(&tgt_array[dst_add], src_array, wgts, src_add, grad1_lat, grad1_lon, grad1_latlon);
	    }
        }
    } /* grid_loop1 */

  free(grad1_lat);
  free(grad1_lon);
  free(grad1_latlon);

} /* scrip_remap_bicubic */
