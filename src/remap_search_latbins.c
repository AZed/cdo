#include "cdo.h"
#include "remap.h"

void calc_bin_addr(long gridsize, long nbins, const restr_t* restrict bin_lats, const restr_t* restrict cell_bound_box, int* restrict bin_addr)
{
  long n, n2, nele, nele4;
  restr_t cell_bound_box_lat1, cell_bound_box_lat2;

  for ( n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      bin_addr[n2  ] = gridsize;
      bin_addr[n2+1] = 0;
    }

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  private(n, n2, nele4, cell_bound_box_lat1, cell_bound_box_lat2)  \
  shared(gridsize, nbins, bin_lats, cell_bound_box, bin_addr)
#endif
  for ( nele = 0; nele < gridsize; ++nele )
    {
      nele4 = nele<<2;
      cell_bound_box_lat1 = cell_bound_box[nele4  ];
      cell_bound_box_lat2 = cell_bound_box[nele4+1];
      for ( n = 0; n < nbins; ++n )
	{
	  n2 = n<<1;
	  if ( cell_bound_box_lat1 <= bin_lats[n2+1] &&
	       cell_bound_box_lat2 >= bin_lats[n2  ] )
	    {
	      /*
#if defined(_OPENMP)
	      if ( nele < bin_addr[n2  ] || nele > bin_addr[n2+1] )
#pragma omp critical
#endif
	      */
		{
		  bin_addr[n2  ] = MIN(nele, bin_addr[n2  ]);
		  bin_addr[n2+1] = MAX(nele, bin_addr[n2+1]);
		}
	    }
	}
    }
}
/*
static
void calc_bin_addr(long gridsize, long nbins, const restr_t* restrict bin_lats, const restr_t* restrict cell_bound_box, int* restrict bin_addr)
{
  long n, n2, nele, nele4;

  for ( n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      bin_addr[n2  ] = gridsize;
      bin_addr[n2+1] = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  private(nele4)	\
  shared(n2, gridsize, bin_lats, cell_bound_box, bin_addr)
#endif
      for ( nele = 0; nele < gridsize; ++nele )
	{
	  nele4 = nele<<2;

	  if ( cell_bound_box[nele4  ] <= bin_lats[n2+1] &&
	       cell_bound_box[nele4+1] >= bin_lats[n2  ] )
	    {
	      bin_addr[n2  ] = MIN(nele, bin_addr[n2  ]);
	      bin_addr[n2+1] = MAX(nele, bin_addr[n2+1]);
	    }
	}
    }
}
*/

void calc_lat_bins(remapgrid_t* src_grid, remapgrid_t* tgt_grid, int map_type)
{
  long nbins;
  long n;      /* Loop counter                  */
  long n2;
  double dlat;                /* lat/lon intervals for search bins  */
  restr_t *bin_lats = NULL;

  nbins = src_grid->num_srch_bins;
  dlat = PI/nbins;

  if ( cdoVerbose ) cdoPrint("Using %d latitude bins to restrict search.", nbins);

  if ( nbins > 0 )
    {
      bin_lats = src_grid->bin_lats = realloc(src_grid->bin_lats, 2*nbins*sizeof(restr_t));

      for ( n = 0; n < nbins; ++n )
	{
	  n2 = n<<1;
	  bin_lats[n2  ] = RESTR_SCALE((n  )*dlat - PIH);
	  bin_lats[n2+1] = RESTR_SCALE((n+1)*dlat - PIH);
	}

      src_grid->bin_addr = realloc(src_grid->bin_addr, 2*nbins*sizeof(int));

      calc_bin_addr(src_grid->size, nbins, bin_lats, src_grid->cell_bound_box, src_grid->bin_addr);

      if ( map_type == MAP_TYPE_CONSERV || map_type == MAP_TYPE_CONSPHERE )
	{
	  tgt_grid->bin_addr = realloc(tgt_grid->bin_addr, 2*nbins*sizeof(int));

	  calc_bin_addr(tgt_grid->size, nbins, bin_lats, tgt_grid->cell_bound_box, tgt_grid->bin_addr);

	  free(src_grid->bin_lats); src_grid->bin_lats = NULL;
	}
   }

  if ( map_type == MAP_TYPE_CONSPHERE )
    {
      free(tgt_grid->cell_bound_box); tgt_grid->cell_bound_box = NULL;
    }
 
  if ( map_type == MAP_TYPE_DISTWGT )
    {
      free(src_grid->cell_bound_box); src_grid->cell_bound_box = NULL;
    }
}
