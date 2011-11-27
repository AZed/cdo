/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2011 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

      Interpolate remapcon        First order conservative remapping
      Interpolate remapcon2       Second order conservative remapping
      Interpolate remapbil        Bilinear interpolation
      Interpolate remapbic        Bicubic interpolation
      Interpolate remapdis        Distance-weighted averaging
      Interpolate remapnn         Nearest neighbor remapping
      Interpolate remaplaf        Largest area fraction remapping
      Genweights  gencon          Generate first order conservative remap weights
      Genweights  gencon2         Generate second order conservative remap weights
      Genweights  genbil          Generate bilinear interpolation weights
      Genweights  genbic          Generate bicubic interpolation weights
      Genweights  gendis          Generate distance-weighted averaging weights
      Genweights  gennn           Generate nearest neighbor weights
      Genweights  genlaf          Generate largest area fraction weights
      Remap       remap           SCRIP grid remapping
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "remap.h"
#include "grid.h"


enum {REMAPCON, REMAPCON2, REMAPBIL, REMAPBIC, REMAPDIS, REMAPNN, REMAPLAF, REMAPSUM,
      GENCON, GENCON2, GENBIL, GENBIC, GENDIS, GENNN, GENLAF, REMAPXXX};

enum {HEAP_SORT, MERGE_SORT};

static
void get_map_type(int operfunc, int *map_type, int *submap_type, int *remap_order)
{
  switch ( operfunc )
    {
    case REMAPCON:
    case GENCON:
      *map_type = MAP_TYPE_CONSERV;
      *remap_order = 1;
      break;
    case REMAPCON2:
    case GENCON2:
      *map_type = MAP_TYPE_CONSERV;
      *remap_order = 2;
      break;
    case REMAPLAF:
    case GENLAF:
      *map_type = MAP_TYPE_CONSERV;
      *submap_type = SUBMAP_TYPE_LAF;
      break;
    case REMAPSUM:
      *map_type = MAP_TYPE_CONSERV;
      *submap_type = SUBMAP_TYPE_SUM;
      break;
    case REMAPBIL:
    case GENBIL:
      *map_type = MAP_TYPE_BILINEAR;
      break;
    case REMAPBIC:
    case GENBIC:
      *map_type = MAP_TYPE_BICUBIC;
      break;
    case REMAPDIS:
    case GENDIS:
      *map_type = MAP_TYPE_DISTWGT;
      break;
    case REMAPNN:
    case GENNN:
      *map_type = MAP_TYPE_DISTWGT1;
      break;
    default:
      cdoAbort("Unknown mapping method");
    }
}

static
int maptype2operfunc(int map_type, int submap_type, int remap_order)
{
  int operfunc = -1;

  if ( map_type == MAP_TYPE_CONSERV )
    {
      if ( submap_type == SUBMAP_TYPE_LAF )
	{
	  operfunc = REMAPLAF;
	  cdoPrint("Using remaplaf");
	}
      else
	{
	  if ( remap_order == 2 )
	    {
	      operfunc = REMAPCON2;
	      cdoPrint("Using remapcon2");
	    }
	  else
	    {
	      operfunc = REMAPCON;
	      cdoPrint("Using remapcon");
	    }
	}
    }
  else if ( map_type == MAP_TYPE_BILINEAR )
    {
      operfunc = REMAPBIL;
      cdoPrint("Using remapbil");
    }
  else if ( map_type == MAP_TYPE_BICUBIC )
    {
      operfunc = REMAPBIC;
      cdoPrint("Using remapbic");
    }
  else if ( map_type == MAP_TYPE_DISTWGT )
    {
      operfunc = REMAPDIS;
      cdoPrint("Using remapdis");
    }
  else if ( map_type == MAP_TYPE_DISTWGT1 )
    {
      operfunc = REMAPNN;
      cdoPrint("Using remapnn");
    }
  else
    cdoAbort("Unsupported mapping method (map_type = %d)", map_type);

  return (operfunc);
}

double remap_threshhold = 2;
int remap_test = 0;
int remap_order = 1;
int remap_restrict_type = RESTRICT_LATITUDE;
int remap_store_link_fast = TRUE;
int remap_non_global = FALSE;
int remap_num_srch_bins = 180;
int lremap_num_srch_bins = FALSE;
int remap_extrapolate = FALSE;
int lextrapolate = FALSE;
int max_remaps = 0;
int sort_mode = HEAP_SORT;

static
void get_remap_env(void)
{
  char *envstr;

  envstr = getenv("MAX_REMAPS");
  if ( envstr )
    {
      int ival;
      ival = atoi(envstr);
      if ( ival > 0 )
	{
	  max_remaps = ival;
	  if ( cdoVerbose )
	    cdoPrint("Set MAX_REMAPS to %d", max_remaps);
	}
    }

  envstr = getenv("REMAP_MAX_ITER");
  if ( envstr )
    {
      int ival;
      ival = atoi(envstr);
      if ( ival > 0 )
	{
	  remap_set_max_iter(ival);
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_MAX_ITER to %d", ival);
	}
    }
  /*
  envstr = getenv("REMAP_ORDER");
  if ( envstr )
    {
      int ival;
      ival = atoi(envstr);
      if ( ival > 0 )
	{
	  remap_order = ival;
	  if ( remap_order == 0 ) remap_order = 1;
	  if ( remap_order != 1 && remap_order != 2 )
	    cdoAbort("REMAP_ORDER must be 1 or 2");
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_ORDER to %d", remap_order);
	}
    }
  */
  envstr = getenv("REMAP_TEST");
  if ( envstr )
    {
      int ival;
      ival = atoi(envstr);
      if ( ival > 0 )
	{
	  remap_test = ival;
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_TEST to %d", remap_test);
	}
    }

#if defined (_OPENMP)
  if ( ompNumThreads == 1 )
    sort_mode = HEAP_SORT;
  else
    sort_mode = MERGE_SORT;
#endif

  envstr = getenv("REMAP_SORT_MODE");
  if ( envstr )
    {
      if      ( strcmp(envstr, "heap")  == 0 ) sort_mode = HEAP_SORT;
      else if ( strcmp(envstr, "merge") == 0 ) sort_mode = MERGE_SORT;

      if ( cdoVerbose )
	{
	  if      ( sort_mode == HEAP_SORT )
	    cdoPrint("Set sort_mode to HEAP_SORT");
	  else if ( sort_mode == MERGE_SORT )
	    cdoPrint("Set sort_mode to MERGE_SORT");
	}
    }

  envstr = getenv("REMAP_RESTRICT_TYPE");
  if ( envstr )
    {
      if      ( strcmp(envstr, "latitude") == 0 ) remap_restrict_type = RESTRICT_LATITUDE;
      else if ( strcmp(envstr, "latlon")   == 0 ) remap_restrict_type = RESTRICT_LATLON;

      if ( cdoVerbose )
	{
	  if      ( remap_restrict_type == RESTRICT_LATITUDE )
	    cdoPrint("Set REMAP_RESTRICT_TYPE to latitude");
	  else if ( remap_restrict_type == RESTRICT_LATLON )
	    cdoPrint("Set REMAP_RESTRICT_TYPE to latlon");
	}
    }

  envstr = getenv("REMAP_THRESHHOLD");
  if ( envstr )
    {
      double fval;
      fval = atof(envstr);
      if ( fval > 0 )
	{
	  remap_threshhold = fval;
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_THRESHHOLD to %g", remap_threshhold);
	}
    }

  envstr = getenv("REMAP_NUM_SRCH_BINS");
  if ( envstr )
    {
      int ival;
      ival = atoi(envstr);
      if ( ival > 0 )
	{
	  remap_num_srch_bins = ival;
	  lremap_num_srch_bins = TRUE;
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_NUM_SRCH_BINS to %d", remap_num_srch_bins);
	}
    }

  envstr = getenv("REMAP_NON_GLOBAL");
  if ( envstr )
    {
      int ival;
      ival = atoi(envstr);
      if ( ival >= 0 )
	{
	  remap_non_global = ival;
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_NON_GLOBAL to %d", remap_non_global);
	}
    }

  envstr = getenv("REMAP_STORE_LINK_FAST");
  if ( envstr )
    {
      int ival;
      ival = atoi(envstr);
      if ( ival >= 0 )
	{
	  remap_store_link_fast = ival;
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_STORE_LINK_FAST to %d", remap_store_link_fast);
	}
    }

  envstr = getenv("REMAP_EXTRAPOLATE");
  if ( envstr )
    {
      int ival;
      ival = atoi(envstr);
      if ( *envstr )
	{
	  if ( memcmp(envstr, "ON", 2) == 0 || memcmp(envstr, "on", 2) == 0 )
	    {
	      lextrapolate = TRUE;
	      remap_extrapolate = TRUE;
	    }
	  else if ( memcmp(envstr, "OFF", 3) == 0 || memcmp(envstr, "off", 3) == 0 )
	    {
	      lextrapolate = TRUE;
	      remap_extrapolate = FALSE;
	    }
	  else
	    cdoWarning("Environment variable REMAP_EXTRAPOLATE has wrong value!");

	  if ( cdoVerbose )
	    {
	      if ( remap_extrapolate == TRUE )
		cdoPrint("Extrapolation enabled!");
	      else if ( remap_extrapolate == FALSE )
		cdoPrint("Extrapolation disabled!");
	    }
	}
    }
}


void *Remap(void *argument)
{
  int operatorID;
  int operfunc;
  int streamID1, streamID2 = -1;
  int nrecs, ngrids;
  int nzaxis, zaxisID, zaxissize;
  int nvars;
  int index;
  int tsID, recID, varID, levelID;
  int gridsize, gridsize2;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int gridID1 = -1, gridID2;
  int gridtype;
  int nmiss1, nmiss2, i, j, r;
  int *imask = NULL;
  int nremaps = 0;
  int norm_opt = NORM_OPT_NONE;
  int map_type = -1;
  int submap_type = SUBMAP_TYPE_NONE;
  int need_gradiants = FALSE;
  int non_global;
  int lgridboxinfo = TRUE;
  int grid1sizemax;
  short *remapgrids = NULL;
  char varname[CDI_MAX_NAME];
  double missval;
  double *array1 = NULL, *array2 = NULL;
  double *grad1_lat = NULL, *grad1_lon = NULL, *grad1_latlon = NULL;
  remap_t *remaps;
  char *envstr;
  char *remap_file = NULL;
  int lwrite_remap;

  cdoInitialize(argument);

  cdoOperatorAdd("remapcon",    REMAPCON,    0, NULL);
  cdoOperatorAdd("remapcon2",   REMAPCON2,   0, NULL);
  cdoOperatorAdd("remapbil",    REMAPBIL,    0, NULL);
  cdoOperatorAdd("remapbic",    REMAPBIC,    0, NULL);
  cdoOperatorAdd("remapdis",    REMAPDIS,    0, NULL);
  cdoOperatorAdd("remapnn",     REMAPNN,     0, NULL);
  cdoOperatorAdd("remaplaf",    REMAPLAF,    0, NULL);
  cdoOperatorAdd("remapsum",    REMAPSUM,    0, NULL);
  cdoOperatorAdd("gencon",      GENCON,      1, NULL);
  cdoOperatorAdd("gencon2",     GENCON2,     1, NULL);
  cdoOperatorAdd("genbil",      GENBIL,      1, NULL);
  cdoOperatorAdd("genbic",      GENBIC,      1, NULL);
  cdoOperatorAdd("gendis",      GENDIS,      1, NULL);
  cdoOperatorAdd("gennn",       GENNN,       1, NULL);
  cdoOperatorAdd("genlaf",      GENLAF,      1, NULL);
  cdoOperatorAdd("remap",       REMAPXXX,    0, NULL);

  operatorID = cdoOperatorID();
  operfunc   = cdoOperatorF1(operatorID);
  lwrite_remap = cdoOperatorF2(operatorID);

  if ( operfunc == REMAPDIS || operfunc == GENDIS ||
       operfunc == REMAPNN  || operfunc == GENNN )
    remap_extrapolate = TRUE;

  get_remap_env();

  if ( cdoVerbose )
    {
      if ( remap_extrapolate == TRUE )
	cdoPrint("Extrapolation enabled!");
      else if ( remap_extrapolate == FALSE )
	cdoPrint("Extrapolation disabled!");
    }

  if ( operfunc == REMAPXXX )
    {
      operatorInputArg("grid description file or name, remap file (SCRIP netCDF)");
      operatorCheckArgc(2);
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
      remap_file = operatorArgv()[1];
    }
  else
    {
      operatorInputArg("grid description file or name");
      operatorCheckArgc(1);
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  remapgrids = (short *) malloc(ngrids*sizeof(short));
  for ( index = 0; index < ngrids; index++ )
    {
      remapgrids[index] = TRUE;

      gridID1 = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);

      if ( gridtype != GRID_LONLAT      &&
	   gridtype != GRID_GAUSSIAN    &&
	   gridtype != GRID_LCC         &&
	   gridtype != GRID_LAEA        &&
	   gridtype != GRID_SINUSOIDAL  &&
	   gridtype != GRID_GME         &&
	   gridtype != GRID_REFERENCE   &&
	   gridtype != GRID_CURVILINEAR &&
	   gridtype != GRID_UNSTRUCTURED )
	{
	  if ( gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED )
	    cdoAbort("Unsupported grid type: %s, use CDO option -R to convert reduced to regular grid!",
		     gridNamePtr(gridInqType(gridID1)));
	  else if ( gridInqType(gridID1) == GRID_GENERIC && gridInqSize(gridID1) == 1 )
	    remapgrids[index] = FALSE;
	  else
	    cdoAbort("Unsupported grid type: %s", gridNamePtr(gridInqType(gridID1)));
	}

      if ( remapgrids[index] )
	vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  for ( index = 0; index < ngrids; index++ )
    {
      if ( remapgrids[index] == TRUE ) break;
    }

  if ( index == ngrids )
    cdoAbort("No remappable grid found!");

  gridID1 = vlistGrid(vlistID1, index);

  if ( max_remaps == 0 )
    {
      nzaxis = vlistNzaxis(vlistID1);
      for ( index = 0; index < nzaxis; index++ )
        {
	  zaxisID = vlistZaxis(vlistID1, index);
	  zaxissize = zaxisInqSize(zaxisID);
          if ( zaxissize > max_remaps ) max_remaps = zaxissize;
	}

      nvars = vlistNvars(vlistID1);
      if ( nvars > max_remaps ) max_remaps = nvars;

      max_remaps++;

      if ( cdoVerbose )
        cdoPrint("Set max_remaps to %d", max_remaps);
    }

  remaps = (remap_t *) malloc(max_remaps*sizeof(remap_t));
  for ( r = 0; r < max_remaps; r++ )
    {
      remaps[r].gridID   = -1;
      remaps[r].gridsize = 0;
      remaps[r].nmiss    = 0;
    }

  if ( operfunc == REMAPXXX )
    {
      int gridsize2;

      read_remap_scrip(remap_file, gridID1, gridID2, &map_type, &submap_type, 
		       &remap_order, &remaps[0].grid, &remaps[0].vars);
      nremaps = 1;
      gridsize = remaps[0].grid.grid1_size;
      remaps[0].gridID = gridID1;
      remaps[0].gridsize = gridInqSize(gridID1);
      remaps[0].nmiss = 0;

      if ( map_type == MAP_TYPE_DISTWGT || map_type == MAP_TYPE_DISTWGT1 )
	{
	  if ( !lextrapolate ) remap_extrapolate = TRUE;
	}

      if ( gridIsCircular(gridID1) && !lextrapolate ) remap_extrapolate = TRUE;
      non_global = remap_non_global || !gridIsCircular(gridID1);
      if ( !remap_extrapolate && gridInqSize(gridID1) > 1 &&
	   (map_type == MAP_TYPE_DISTWGT || map_type == MAP_TYPE_DISTWGT1) &&
	   ((gridInqType(gridID1) == GRID_LONLAT && gridIsRotated(gridID1)) ||
	    (gridInqType(gridID1) == GRID_LONLAT && non_global) ||
	    (gridInqType(gridID1) == GRID_LCC) ||
	    (gridInqType(gridID1) == GRID_LAEA) ||
	    (gridInqType(gridID1) == GRID_SINUSOIDAL) ||
	    (gridInqType(gridID1) == GRID_CURVILINEAR && non_global)) )
	{
	  remaps[0].gridsize += 4*(gridInqXsize(gridID1)+2) + 4*(gridInqYsize(gridID1)+2);
	  remaps[0].grid.non_global = TRUE;
	}

      if ( gridInqType(gridID1) == GRID_GME ) gridsize = remaps[0].grid.grid1_nvgp;

      if ( gridsize != remaps[0].gridsize )
	cdoAbort("Size of source grid and weights from %s differ!", remap_file);

      if ( gridInqType(gridID1) == GRID_GME ) gridsize = remaps[0].grid.grid1_size;

      for ( i = 0; i < gridsize; i++ )
        if ( remaps[0].grid.grid1_mask[i] == FALSE )
          remaps[0].nmiss++;

      gridsize2 = gridInqSize(gridID2);
      if ( gridInqType(gridID2) == GRID_GME )
	{
	  int gridID2_gme;
	  int isize = 0;
	  remaps[0].grid.grid2_nvgp = gridInqSize(gridID2);
	  remaps[0].grid.grid2_vgpm = (int *) realloc(remaps[0].grid.grid2_vgpm,
						      gridInqSize(gridID2)*sizeof(int));
	  gridID2_gme = gridToUnstructured(gridID2);
	  gridInqMaskGME(gridID2_gme, remaps[0].grid.grid2_vgpm);
	  for ( i = 0; i < gridsize2; ++i )
	    if ( remaps[0].grid.grid2_vgpm[i] ) isize++;
	  gridsize2 = isize;
	}
      /*
      printf("grid2 %d %d %d\n", gridsize2, remaps[0].grid.grid2_nvgp, remaps[0].grid.grid2_size);
      */
      if ( remaps[0].grid.grid2_size != gridsize2 )
	cdoAbort("Size of target grid and weights from %s differ!", remap_file);

      operfunc = maptype2operfunc(map_type, submap_type, remap_order);

      if ( remap_test ) reorder_links(&remaps[0].vars);
    }

  get_map_type(operfunc, &map_type, &submap_type, &remap_order);

  if ( map_type == MAP_TYPE_CONSERV )
    {
      norm_opt = NORM_OPT_FRACAREA;

      envstr = getenv("NORMALIZE_OPT");

      if ( envstr )
        {
	  if      ( memcmp(envstr, "frac", 4) == 0 )
	    norm_opt = NORM_OPT_FRACAREA;
	  else if ( memcmp(envstr, "dest", 4) == 0 )
	    norm_opt = NORM_OPT_DESTAREA;
	  else if ( memcmp(envstr, "none", 4) == 0 )
	    norm_opt = NORM_OPT_NONE;
	  else
	    cdoWarning("NORMALIZE_OPT=%s unsupported!", envstr);
	}

      if ( cdoVerbose )
        {
	  if ( norm_opt == NORM_OPT_FRACAREA )
	    cdoPrint("Normalization option: frac");
	  else if ( norm_opt == NORM_OPT_DESTAREA )
	    cdoPrint("Normalization option: dest");
	  else
	    cdoPrint("Normalization option: none");
	}
    }

  grid1sizemax = vlistGridsizeMax(vlistID1);

  if ( map_type == MAP_TYPE_BICUBIC ) need_gradiants = TRUE;
  if ( map_type == MAP_TYPE_CONSERV && remap_order == 2 )
    {
      if ( cdoVerbose ) cdoPrint("Second order remapping");
      need_gradiants = TRUE;
    }
  else
    remap_order = 1;

  if ( need_gradiants )
    {
      grad1_lat    = (double *) malloc(grid1sizemax*sizeof(double));
      grad1_lon    = (double *) malloc(grid1sizemax*sizeof(double));
      grad1_latlon = (double *) malloc(grid1sizemax*sizeof(double));
    }

  array1 = (double *) malloc(grid1sizemax*sizeof(double));
  imask  = (int *) malloc(grid1sizemax*sizeof(int));

  gridsize = gridInqSize(gridID2);
  array2   = (double *) malloc(gridsize*sizeof(double));

  if ( ! lwrite_remap )
    {
      streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

      streamDefVlist(streamID2, vlistID2);
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      if ( ! lwrite_remap ) 
	streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss1);

	  gridID1 = vlistInqVarGrid(vlistID1, varID);

	  if ( remapgrids[vlistGridIndex(vlistID1, gridID1)] == FALSE )
	    {
	      if ( lwrite_remap ) continue;
	      else
		{
		  nmiss2 = nmiss1;
		  *array2 = *array1;
		  goto SKIPVAR;
		}
	    }

	  if ( map_type != MAP_TYPE_CONSERV && 
	       gridInqType(gridID1) == GRID_GME && gridInqType(gridID2) == GRID_GME )
	    cdoAbort("Only conservative remapping is available to remap between GME grids!");
	  /*
	  if ( gridIsRotated(gridID1) && map_type != MAP_TYPE_CONSERV )
	    cdoAbort("Only conservative remapping is available for rotated grids!");
	  */
	  missval = vlistInqVarMissval(vlistID1, varID);
	  gridsize = gridInqSize(gridID1);

	  if ( gridIsCircular(gridID1) && !lextrapolate ) remap_extrapolate = TRUE;
	  non_global = remap_non_global || !gridIsCircular(gridID1);
	  if ( !remap_extrapolate && gridInqSize(gridID1) > 1 &&
	       (map_type == MAP_TYPE_DISTWGT || map_type == MAP_TYPE_DISTWGT1) &&
	       ((gridInqType(gridID1) == GRID_LONLAT && gridIsRotated(gridID1)) ||
		(gridInqType(gridID1) == GRID_LONLAT && non_global) ||
		(gridInqType(gridID1) == GRID_LCC) ||
		(gridInqType(gridID1) == GRID_LAEA) ||
		(gridInqType(gridID1) == GRID_SINUSOIDAL) ||
		(gridInqType(gridID1) == GRID_CURVILINEAR && non_global)) )
	    {
	      int gridsize_new;
	      int nx, ny;
	      nx = gridInqXsize(gridID1);
	      ny = gridInqYsize(gridID1);
	      gridsize_new = gridsize + 4*(nx+2) + 4*(ny+2);
	      if ( gridsize_new > grid1sizemax )
		{
		  grid1sizemax = gridsize_new;
		  array1 = (double *) realloc(array1, grid1sizemax*sizeof(double));
		  imask  = (int *) realloc(imask, grid1sizemax*sizeof(int));

		  if ( need_gradiants )
		    {
		      grad1_lat    = (double *) realloc(grad1_lat, grid1sizemax*sizeof(double));
		      grad1_lon    = (double *) realloc(grad1_lon, grid1sizemax*sizeof(double));
		      grad1_latlon = (double *) realloc(grad1_latlon, grid1sizemax*sizeof(double));
		    }
		}
	      
	      for ( j = ny-1; j >= 0; j-- )
		for ( i = nx-1; i >= 0; i-- )
		  array1[(j+2)*(nx+4)+i+2] = array1[j*nx+i];

	      for ( j = 0; j < ny+4; j++ ) array1[j*(nx+4)+0]      = missval;
	      for ( j = 0; j < ny+4; j++ ) array1[j*(nx+4)+1]      = missval;
	      for ( j = 0; j < ny+4; j++ ) array1[j*(nx+4)+nx+2]   = missval;
	      for ( j = 0; j < ny+4; j++ ) array1[j*(nx+4)+nx+3]   = missval;
	      for ( i = 0; i < nx+4; i++ ) array1[     0*(nx+4)+i] = missval;
	      for ( i = 0; i < nx+4; i++ ) array1[     1*(nx+4)+i] = missval;
	      for ( i = 0; i < nx+4; i++ ) array1[(ny+2)*(nx+4)+i] = missval;
	      for ( i = 0; i < nx+4; i++ ) array1[(ny+3)*(nx+4)+i] = missval;

	      gridsize = gridsize_new;
	      nmiss1 += 4*(nx+2) + 4*(ny+2);
	    }

	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(array1[i], missval) )
	      imask[i] = FALSE;
	    else
	      imask[i] = TRUE;

	  for ( r = nremaps-1; r >= 0; r-- )
	    {
	      if ( gridID1 == remaps[r].gridID && nmiss1 == remaps[r].nmiss )
		{
		  if ( memcmp(imask, remaps[r].grid.grid1_mask, remaps[r].grid.grid1_size*sizeof(int)) == 0 )
		    break;
		}	      
	    }

	  if ( cdoVerbose && r >= 0 ) cdoPrint("Using remap %d", r);

	  if ( r < 0 )
	    {
	      if ( nremaps < max_remaps )
		{
		  r = nremaps;
		  nremaps++;
		}
	      else
		{
		  r = nremaps - 1;
		}

	      if ( remaps[r].gridID != gridID1 )
		{
		  if ( gridIsCircular(gridID1) && !lextrapolate ) remap_extrapolate = TRUE;
		  remaps[r].grid.non_global = FALSE;
		  non_global = remap_non_global || !gridIsCircular(gridID1);
		  if ( !remap_extrapolate && gridInqSize(gridID1) > 1 &&
		       (map_type == MAP_TYPE_DISTWGT || map_type == MAP_TYPE_DISTWGT1) &&
		       ((gridInqType(gridID1) == GRID_LONLAT && gridIsRotated(gridID1)) ||
			(gridInqType(gridID1) == GRID_LONLAT && non_global) ||
			(gridInqType(gridID1) == GRID_LCC) ||
			(gridInqType(gridID1) == GRID_LAEA) ||
			(gridInqType(gridID1) == GRID_SINUSOIDAL) ||
			(gridInqType(gridID1) == GRID_CURVILINEAR && non_global)) )
		    {
		      remaps[r].grid.non_global = TRUE;
		    }
		  /*
		    remaps[r].grid.luse_grid1_area = FALSE;
		    remaps[r].grid.luse_grid2_area = FALSE;
		  */
		  if ( gridInqType(gridID1) != GRID_UNSTRUCTURED && lremap_num_srch_bins == FALSE )
		    {
		      if ( !remap_extrapolate && (map_type == MAP_TYPE_DISTWGT || map_type == MAP_TYPE_DISTWGT1) )
			{
			  remap_num_srch_bins = 1;
			}
		      else
			{
			  int maxbins = 720;
			  int ysize1 = gridInqYsize(gridID1);
			  remap_num_srch_bins = ysize1/2 + ysize1%2;
			  if ( remap_num_srch_bins > maxbins ) remap_num_srch_bins = maxbins;
			  if ( remap_num_srch_bins < 1 )       remap_num_srch_bins = 1;
			}
		    }

		  remaps[r].grid.threshhold    = remap_threshhold;
		  remaps[r].grid.restrict_type = remap_restrict_type;
		  remaps[r].grid.num_srch_bins = remap_num_srch_bins;
		  remaps[r].grid.pinit = FALSE;

		  remaps[r].vars.norm_opt = norm_opt;
		  remaps[r].vars.pinit = FALSE;
		  
		  /* initialize grid information for both grids */
		  remapGridInit(map_type, remap_extrapolate, gridID1, gridID2, &remaps[r].grid);
		  remaps[r].grid.store_link_fast = remap_store_link_fast;
		}

	      remaps[r].gridID = gridID1;
	      remaps[r].nmiss  = nmiss1;

	      if ( gridInqType(gridID1) == GRID_GME )
		{
		  j = 0;
		  for ( i = 0; i < gridsize; i++ )
		    if ( remaps[r].grid.grid1_vgpm[i] ) imask[j++] = imask[i];
		}

	      memcpy(remaps[r].grid.grid1_mask, imask, remaps[r].grid.grid1_size*sizeof(int));

	      if ( map_type == MAP_TYPE_CONSERV )
		{
		  memset(remaps[r].grid.grid1_area, 0, remaps[r].grid.grid1_size*sizeof(double));
		  memset(remaps[r].grid.grid1_frac, 0, remaps[r].grid.grid1_size*sizeof(double));
		  memset(remaps[r].grid.grid2_area, 0, remaps[r].grid.grid2_size*sizeof(double));
		}
	      memset(remaps[r].grid.grid2_frac, 0, remaps[r].grid.grid2_size*sizeof(double));

	      /* initialize some remapping variables */
	      remapVarsInit(map_type, &remaps[r].grid, &remaps[r].vars);

	      if      ( map_type == MAP_TYPE_CONSERV  ) remap_conserv(&remaps[r].grid, &remaps[r].vars);
	      else if ( map_type == MAP_TYPE_BILINEAR ) remap_bilin(&remaps[r].grid, &remaps[r].vars);
	      else if ( map_type == MAP_TYPE_BICUBIC  ) remap_bicub(&remaps[r].grid, &remaps[r].vars);
	      else if ( map_type == MAP_TYPE_DISTWGT  ) remap_distwgt(&remaps[r].grid, &remaps[r].vars);
	      else if ( map_type == MAP_TYPE_DISTWGT1 ) remap_distwgt1(&remaps[r].grid, &remaps[r].vars);

	      if ( remaps[r].vars.num_links != remaps[r].vars.max_links )
		resize_remap_vars(&remaps[r].vars, remaps[r].vars.num_links-remaps[r].vars.max_links);

	      if ( sort_mode == MERGE_SORT )
		{ /* 
		  ** use a combination of the old sort_add and a split and merge approach.
                  ** The chunk size is determined by MERGE_SORT_LIMIT_SIZE in remaplib.c. 
		  ** OpenMP parallelism is supported
		  */   
		  sort_iter(remaps[r].vars.num_links, remaps[r].vars.num_wts,
			    remaps[r].vars.grid2_add, remaps[r].vars.grid1_add,
			    remaps[r].vars.wts, ompNumThreads);
		}
	      else
		{ /* use a pure heap sort without any support of parallelism */
		  sort_add(remaps[r].vars.num_links, remaps[r].vars.num_wts,
			   remaps[r].vars.grid2_add, remaps[r].vars.grid1_add,
			   remaps[r].vars.wts);
		}
	      	      
	      if ( lwrite_remap ) goto WRITE_REMAP;

	      if ( remap_test ) reorder_links(&remaps[r].vars);
	    }

	  if ( gridInqType(gridID1) == GRID_GME )
	    {
	      j = 0;
	      for ( i = 0; i < gridsize; i++ )
		if ( remaps[r].grid.grid1_vgpm[i] ) array1[j++] = array1[i];
	    }
	  
	  if ( need_gradiants )
	    {
	      if ( remaps[r].grid.grid1_rank != 2 && remap_order == 2 )
		cdoAbort("Second order remapping is only available for 2D grids!");

	      remap_gradients(remaps[r].grid, array1, grad1_lat, grad1_lon, grad1_latlon);
	    }

	  if ( operfunc == REMAPLAF )
	    remap_laf(array2, missval, gridInqSize(gridID2), remaps[r].vars.num_links, remaps[r].vars.wts,
		  remaps[r].vars.num_wts, remaps[r].vars.grid2_add, remaps[r].vars.grid1_add, array1);
	  else if ( operfunc == REMAPSUM )
	    remap_sum(array2, missval, gridInqSize(gridID2), remaps[r].vars.num_links, remaps[r].vars.wts,
		  remaps[r].vars.num_wts, remaps[r].vars.grid2_add, remaps[r].vars.grid1_add, array1);
	  else
	    remap(array2, missval, gridInqSize(gridID2), remaps[r].vars.num_links, remaps[r].vars.wts,
		  remaps[r].vars.num_wts, remaps[r].vars.grid2_add, remaps[r].vars.grid1_add,
		  array1, grad1_lat, grad1_lon, grad1_latlon, remaps[r].vars.links);

	  gridsize2 = gridInqSize(gridID2);

	  if ( operfunc == REMAPCON || operfunc == REMAPCON2 )
	    {
	      double grid2_err;

	      if ( remaps[r].vars.norm_opt == NORM_OPT_NONE )
		{
		  for ( i = 0; i < gridsize2; i++ )
		    {
		      grid2_err = remaps[r].grid.grid2_frac[i]*remaps[r].grid.grid2_area[i];
		      if ( fabs(grid2_err) > 0 )
			array2[i] = array2[i]/grid2_err;
		      else
			array2[i] = missval;
		    }
		}
	      else if ( remaps[r].vars.norm_opt == NORM_OPT_DESTAREA )
		{
		  for ( i = 0; i < gridsize2; i++ )
		    {
		      if ( fabs(remaps[r].grid.grid2_frac[i]) > 0 )
			array2[i] = array2[i]/remaps[r].grid.grid2_frac[i];
		      else
			array2[i] = missval;
		    }
		}
	    }

	  if ( operfunc == REMAPSUM )
	  {
	    double array1sum = 0;
	    double array2sum = 0;
   
	    for ( i = 0; i < gridsize; i++ )
	      printf("1 %d %g %g %g %g\n", i, array1[i], remaps[r].grid.grid1_frac[i],remaps[r].grid.grid1_area[i],remaps[r].grid.grid1_frac[i]);
	    for ( i = 0; i < gridsize; i++ )
	      array1sum += remaps[r].grid.grid1_area[i];

	    for ( i = 0; i < gridsize2; i++ )
	      printf("2 %d %g %g %g %g\n", i, array2[i], remaps[r].grid.grid2_frac[i],remaps[r].grid.grid2_area[i],remaps[r].grid.grid2_frac[i]);
	    for ( i = 0; i < gridsize2; i++ )
	      array2sum += remaps[r].grid.grid2_area[i];

	    printf("array1sum %g, array2sum %g\n", array1sum, array2sum);
	  }

	  vlistInqVarName(vlistID1, varID, varname);
	  if ( operfunc == REMAPCON || operfunc == REMAPCON2 )
	    if ( strcmp(varname, "gridbox_area") == 0 )
	      {
		double array1sum = 0;
		double array2sum = 0;

		for ( i = 0; i < gridsize; i++ )
		  array1sum += array1[i];

		for ( i = 0; i < gridsize2; i++ )
		  array2sum += remaps[r].grid.grid2_area[i];

		for ( i = 0; i < gridsize2; i++ )
		  array2[i] = remaps[r].grid.grid2_area[i]/array2sum*array1sum;

		if ( lgridboxinfo )
		  {
		    cdoPrint("%s replaced and scaled to %g", varname, array1sum);
		    lgridboxinfo = FALSE;
		  }
	      }

	  /* calculate some statistics */
	  if ( cdoVerbose )
	    remap_stat(remap_order, remaps[r].grid, remaps[r].vars, array1, array2, missval);

	  if ( gridInqType(gridID2) == GRID_GME )
	    {
	      int ni, nd;
	      ni = gridInqGMEni(gridID2);
	      nd = gridInqGMEnd(gridID2);
	      j = remaps[r].grid.grid2_size;

	      for ( i = gridsize2-1; i >=0 ; i-- )
		if ( remaps[r].grid.grid2_vgpm[i] ) array2[i] = array2[--j];

	      gme_grid_restore(array2, ni, nd);
	    }

	  nmiss2 = 0;
	  for ( i = 0; i < gridsize2; i++ )
	    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss2++;

	SKIPVAR:

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss2);
	}
      tsID++;
    }

  streamClose(streamID2);

  WRITE_REMAP:
 
  if ( lwrite_remap ) 
    write_remap_scrip(cdoStreamName(1), map_type, submap_type, remap_order, remaps[r].grid, remaps[r].vars);

  streamClose(streamID1);

  if ( remapgrids ) free(remapgrids);
  if ( imask )  free(imask);
  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  if ( grad1_latlon ) free(grad1_latlon);
  if ( grad1_lon ) free(grad1_lon);
  if ( grad1_lat ) free(grad1_lat);

  for ( r = 0; r < nremaps; r++ )
    {
      remapVarsFree(&remaps[r].vars);
      remapGridFree(&remaps[r].grid);
    }

  if ( remaps ) free(remaps);

  cdoFinish();

  return (0);
}
