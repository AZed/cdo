#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>

#include "dmemory.h"
#include "cdi.h"
#include "stream_int.h"
#include "file.h"
#include "varscan.h"
#include "datetime.h"
#include "vlist.h"
#include "stream_grb.h"


#if  defined  (HAVE_LIBGRIB_API)
#  include "cgribex.h"      /* gribGetSize, gribRead, gribGetZip */
#  include "gribapi.h"
#  include "grib_api.h"
#endif


extern int cdiInventoryMode;

typedef struct {
  int param;
  int level1;
  int level2;
  int ltype;
} compvar2_t; 


#if  defined  (HAVE_LIBGRIB_API)
static
int gribapiGetGridType(grib_handle *gh)
{
  int gridtype = GRID_GENERIC;
  int gribgridtype;
  long lpar;

    {
      GRIB_CHECK(grib_get_long(gh, "gridDefinitionTemplateNumber", &lpar), 0);
      gribgridtype = (int) lpar;

      switch (gribgridtype)
	{
	case  GRIB2_GTYPE_LATLON:     { GRIB_CHECK(grib_get_long(gh, "Ni", &lpar), 0);
	                                if ( lpar == (long) GRIB_MISSING_LONG ) break;
                                      }
	case  GRIB2_GTYPE_LATLON_ROT: { gridtype = GRID_LONLAT;    break; }
	case  GRIB2_GTYPE_LCC:        { gridtype = GRID_LCC;       break; }
	case  GRIB2_GTYPE_GAUSSIAN:   { GRIB_CHECK(grib_get_long(gh, "Ni", &lpar), 0);
	                                if ( lpar == (long) GRIB_MISSING_LONG )
					  gridtype = GRID_GAUSSIAN_REDUCED;
					else
					  gridtype = GRID_GAUSSIAN;
					break;
                                      }
	case  GRIB2_GTYPE_SPECTRAL:   { gridtype = GRID_SPECTRAL;  break; }
	case  GRIB2_GTYPE_GME:        { gridtype = GRID_GME;       break; }
	case  GRIB2_GTYPE_NUMBER:     { gridtype = GRID_REFERENCE; break; }
	}
    }

  return (gridtype);
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
int gribapiGetIsRotated(grib_handle *gh)
{
  int isRotated = 0;
  int gribgridtype;
  long lpar;

    {
      GRIB_CHECK(grib_get_long(gh, "gridDefinitionTemplateNumber", &lpar), 0);
      gribgridtype = (int) lpar;

      if ( gribgridtype == GRIB2_GTYPE_LATLON_ROT ) isRotated = 1;
    }

  return (isRotated);
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
int gribapiGetZaxisType(long editionNumber, int grib_ltype)
{
  int zaxistype = ZAXIS_GENERIC;

  if ( editionNumber <= 1 )
    {
      zaxistype = grib1ltypeToZaxisType(grib_ltype);
    }
  else
    {
      zaxistype = grib2ltypeToZaxisType(grib_ltype);
    }

  return (zaxistype);
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
int gribapiGetTimeUnits(grib_handle *gh)
{
  int timeunits = -1;
  long lpar;
  size_t len = 8;
  char stepunits[8];
  static int lprint = TRUE;

  GRIB_CHECK(grib_get_string(gh, "stepUnits", stepunits, &len), 0);

  len--;

  if      ( memcmp(stepunits, "s",   len) == 0 ) timeunits = TUNIT_SECOND;
  else if ( memcmp(stepunits, "m",   len) == 0 ) timeunits = TUNIT_MINUTE;
  else if ( memcmp(stepunits, "h",   len) == 0 ) timeunits = TUNIT_HOUR;
  else if ( memcmp(stepunits, "3h",  len) == 0 ) timeunits = TUNIT_3HOURS;
  else if ( memcmp(stepunits, "6h",  len) == 0 ) timeunits = TUNIT_6HOURS;
  else if ( memcmp(stepunits, "12h", len) == 0 ) timeunits = TUNIT_12HOURS;
  else if ( memcmp(stepunits, "D",   len) == 0 ) timeunits = TUNIT_DAY;
  else if ( memcmp(stepunits, "M",   len) == 0 ) timeunits = TUNIT_MONTH;
  else if ( memcmp(stepunits, "Y",   len) == 0 ) timeunits = TUNIT_YEAR;
  else if ( lprint )
    {
      Message("Step units >%s< unsupported!", stepunits);
      lprint = FALSE;
    }

  return (timeunits);
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
int gribapiTimeIsFC(grib_handle *gh)
{
  long editionNumber;
  int isFC = TRUE;

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  if ( editionNumber == 2 )
    {
      long sigofrtime;

      GRIB_CHECK(grib_get_long(gh, "significanceOfReferenceTime", &sigofrtime), 0);

      if ( sigofrtime == 3 ) isFC = FALSE;
    }

  return (isFC);
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
int gribapiGetTsteptype(grib_handle *gh)
{
  int tsteptype = 0;
  int timerange;
  long lpar;
  static int lprint = TRUE;

  if ( gribapiTimeIsFC(gh) )
    {
      int status;
      status = grib_get_long(gh, "stepType", &lpar);
      if ( status == 0 )
	{
	  timerange = (int) lpar;

	  // printf("timerange %d\n", timerange);

	  switch ( timerange )
	    {
	    case  0:  tsteptype = TSTEP_AVG;    break;
	    case  1:  tsteptype = TSTEP_ACCUM;  break;
	    case  2:  tsteptype = TSTEP_MIN;    break;
	    case  3:  tsteptype = TSTEP_MAX;    break;
	    case  4:  tsteptype = TSTEP_DIFF;   break;
	    default:
	      if ( lprint )
		{
		  Message("Time range %d unsupported", timerange);
		  lprint = FALSE;
		}
	    }
	}
    }

  return (tsteptype);
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
void gribapiGetValidityDateTime(grib_handle *gh, int *vdate, int *vtime)
{
  long lpar;
  long sigofrtime = 3;
  long editionNumber;

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  if ( editionNumber == 2 )
    {
      GRIB_CHECK(grib_get_long(gh, "significanceOfReferenceTime", &sigofrtime), 0);
    }

  if ( sigofrtime == 3 )
    {
      GRIB_CHECK(grib_get_long(gh, "dataDate", &lpar), 0);
      *vdate = (int) lpar;
      GRIB_CHECK(grib_get_long(gh, "dataTime", &lpar), 0);
      *vtime = (int) lpar*100;
    }
  else
    {
      GRIB_CHECK(grib_get_long(gh, "validityDate", &lpar), 0);
      *vdate = (int) lpar;
      GRIB_CHECK(grib_get_long(gh, "validityTime", &lpar), 0);
      *vtime = (int) lpar*100;
    }
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
void gribapiGetGrid(grib_handle *gh, grid_t *grid)
{
  long editionNumber;
  int gridtype;
  size_t datasize;
  long numberOfPoints;
  long lpar;

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  gridtype = gribapiGetGridType(gh);
  /*
  if ( streamptr->unreduced && gridtype == GRID_GAUSSIAN_REDUCED )
    {
      gridtype = GRID_GAUSSIAN;
      ISEC2_NumLon = 2*ISEC2_NumLat;
      ISEC4_NumValues = ISEC2_NumLon*ISEC2_NumLat;
    }
  */
  memset(grid, 0, sizeof(grid_t));

  GRIB_CHECK(grib_get_size(gh, "values", &datasize), 0);
  GRIB_CHECK(grib_get_long(gh, "numberOfPoints", &numberOfPoints), 0);

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
      {
	int nlon, nlat;

	GRIB_CHECK(grib_get_long(gh, "Ni", &lpar), 0);
	nlon = lpar;
	GRIB_CHECK(grib_get_long(gh, "Nj", &lpar), 0);
	nlat = lpar;

	if ( numberOfPoints != nlon*nlat )
	  Error("numberOfPoints (%d) and gridSize (%d) differ!",
		(int)numberOfPoints, nlon*nlat);
	grid->size  = numberOfPoints;
	grid->xsize = nlon;
	grid->ysize = nlat;
	grid->xinc  = 0;
	grid->yinc  = 0;
	grid->xdef  = 0;
	GRIB_CHECK(grib_get_double(gh, "longitudeOfFirstGridPointInDegrees", &grid->xfirst), 0);
	GRIB_CHECK(grib_get_double(gh, "longitudeOfLastGridPointInDegrees",  &grid->xlast), 0);
	GRIB_CHECK(grib_get_double(gh, "latitudeOfFirstGridPointInDegrees",  &grid->yfirst), 0);
	GRIB_CHECK(grib_get_double(gh, "latitudeOfLastGridPointInDegrees",   &grid->ylast), 0);
	GRIB_CHECK(grib_get_double(gh, "iDirectionIncrementInDegrees", &grid->xinc), 0);
	if ( gridtype == GRID_LONLAT )
	  GRIB_CHECK(grib_get_double(gh, "jDirectionIncrementInDegrees", &grid->yinc), 0);

	if ( IS_EQUAL(grid->xinc, GRIB_MISSING_DOUBLE) ) grid->xinc = 0;

	/* if ( IS_NOT_EQUAL(grid->xfirst, 0) || IS_NOT_EQUAL(grid->xlast, 0) ) */
	  {
	    if ( grid->xsize > 1 )
	      {
		if ( (grid->xfirst >= grid->xlast) && (grid->xfirst >= 180) ) grid->xfirst -= 360;

		if ( editionNumber <= 1 )
		  {
		    /* correct xinc if necessary */
		    if ( IS_EQUAL(grid->xfirst, 0) && grid->xlast > 354 )
		      {
			double xinc = 360. / grid->xsize;

			if ( fabs(grid->xinc-xinc) > 0.0 )
			  {
			    grid->xinc = xinc;
			    if ( CDI_Debug ) Message("set xinc to %g", grid->xinc);
			  }
		      }
		  }
	      }
	    grid->xdef = 2;
	  }
	grid->ydef = 0;
        /* if ( IS_NOT_EQUAL(grid->yfirst, 0) || IS_NOT_EQUAL(grid->ylast, 0) ) */
	  {
	    if ( grid->ysize > 1 )
	      {
		if ( editionNumber <= 1 )
		  {
		  }
	      }
	    grid->ydef = 2;
	  }
	break;
      }
    case GRID_GAUSSIAN_REDUCED:
      {
	int nlat, i;
	size_t dummy;
	long *pl;

	GRIB_CHECK(grib_get_long(gh, "Nj", &lpar), 0);
	nlat = lpar;

	grid->size   = numberOfPoints;

        grid->rowlon = (int *) malloc(nlat*sizeof(int));
        pl          = (long *) malloc(nlat*sizeof(long));
	dummy       = nlat;
	GRIB_CHECK(grib_get_long_array(gh, "pl", pl, &dummy), 0);
	for ( i = 0; i < nlat; ++i ) grid->rowlon[i] = pl[i];
	free(pl);

	grid->ysize  = nlat;
	grid->xinc   = 0;
	grid->yinc   = 0;
	grid->xdef   = 0;
	GRIB_CHECK(grib_get_double(gh, "longitudeOfFirstGridPointInDegrees", &grid->xfirst), 0);
	GRIB_CHECK(grib_get_double(gh, "longitudeOfLastGridPointInDegrees",  &grid->xlast), 0);
	GRIB_CHECK(grib_get_double(gh, "latitudeOfFirstGridPointInDegrees",  &grid->yfirst), 0);
	GRIB_CHECK(grib_get_double(gh, "latitudeOfLastGridPointInDegrees",   &grid->ylast), 0);
	GRIB_CHECK(grib_get_double(gh, "iDirectionIncrementInDegrees", &grid->xinc), 0);

	if ( IS_EQUAL(grid->xinc, GRIB_MISSING_DOUBLE) ) grid->xinc = 0;

	/* if ( IS_NOT_EQUAL(grid->xfirst, 0) || IS_NOT_EQUAL(grid->xlast, 0) ) */
	  {
	    if ( grid->xsize > 1 )
	      {
		if ( (grid->xfirst > grid->xlast) && (grid->xfirst >= 180) ) grid->xfirst -= 360;

		if ( editionNumber <= 1 )
		  {
		    /* correct xinc if necessary */
		    if ( IS_EQUAL(grid->xfirst, 0) && grid->xlast > 354 )
		      {
			double xinc = 360. / grid->xsize;
		      
			if ( fabs(grid->xinc-xinc) > 0.0 )
			  {
			    grid->xinc = xinc;
			    if ( CDI_Debug ) Message("set xinc to %g", grid->xinc);
			  }
		      }
		  }
	      }
	    grid->xdef   = 2;	    
	  }
	grid->ydef  = 0;
        /* if ( IS_NOT_EQUAL(grid->yfirst, 0) || IS_NOT_EQUAL(grid->ylast, 0) ) */
	  {
	    if ( grid->ysize > 1 )
	      {
		if ( editionNumber <= 1 )
		  {
		  }
	      }
	    grid->ydef   = 2;	    
	  }
	break;
      }
      /*
    case GRID_LCC:
      {
	if ( ISEC4_NumValues != ISEC2_NumLon*ISEC2_NumLat )
	  Error("numberOfPoints (%d) and gridSize (%d) differ!",
		ISEC4_NumValues, ISEC2_NumLon*ISEC2_NumLat);

	grid->size  = ISEC4_NumValues;
	grid->xsize = ISEC2_NumLon;
	grid->ysize = ISEC2_NumLat;

	grid->lcc_xinc      = ISEC2_Lambert_dx;
	grid->lcc_yinc      = ISEC2_Lambert_dy;
	grid->lcc_originLon = ISEC2_FirstLon * 0.001;
	grid->lcc_originLat = ISEC2_FirstLat * 0.001;
	grid->lcc_lonParY   = ISEC2_Lambert_Lov * 0.001;
	grid->lcc_lat1      = ISEC2_Lambert_LatS1 * 0.001;
	grid->lcc_lat2      = ISEC2_Lambert_LatS2 * 0.001;
	grid->lcc_projflag  = ISEC2_Lambert_ProjFlag;
	grid->lcc_scanflag  = ISEC2_ScanFlag;

	grid->xdef   = 0;	    
	grid->ydef   = 0;

	break;
      }
      */
    case GRID_SPECTRAL:
      {
	size_t len = 256;
	char typeOfPacking[256];
	GRIB_CHECK(grib_get_string(gh, "packingType", typeOfPacking, &len), 0);
	grid->lcomplex = 0;
	if ( strncmp(typeOfPacking, "spectral_complex", len) == 0 ) grid->lcomplex = 1;

	grid->size  = datasize;
	GRIB_CHECK(grib_get_long(gh, "J", &lpar), 0);
	grid->trunc = lpar;

	break;
      }
    case GRID_GME:
      {
	grid->size  = numberOfPoints;
	if ( grib_get_long(gh, "nd", &lpar) == 0 ) grid->nd  = lpar;
	if ( grib_get_long(gh, "Ni", &lpar) == 0 ) grid->ni  = lpar;
	if ( grib_get_long(gh, "n2", &lpar) == 0 ) grid->ni2 = lpar;
	if ( grib_get_long(gh, "n3", &lpar) == 0 ) grid->ni3 = lpar;

	break;
      }
    case GRID_REFERENCE:
      {
	char reference_link[8192];
	size_t len = sizeof(reference_link);
	reference_link[0] = 0;

	grid->size  = numberOfPoints;
	if ( grib_get_long(gh, "numberOfGridUsed", &lpar) == 0 )
	  {
	    grid->number   = lpar;
	    if ( grib_get_long(gh, "numberOfGridInReference", &lpar) == 0 ) grid->position = lpar;
	    if ( grib_get_string(gh, "gridDescriptionFile", reference_link, &len) == 0 )
	      {
		if ( strncmp(reference_link, "file://", 7) == 0 )
		  grid->reference = strdupx(reference_link);
	      }
	  }

	break;
      }
    case GRID_GENERIC:
      {
	int nlon = 0, nlat = 0;

	if ( grib_get_long(gh, "Ni", &lpar) == 0 ) nlon = lpar;
	if ( grib_get_long(gh, "Nj", &lpar) == 0 ) nlat = lpar;

	grid->size  = numberOfPoints;

	if ( nlon > 0 && nlat > 0 && nlon*nlat == grid->size )
	  {
	    grid->xsize = nlon;
	    grid->ysize = nlat;
	  }
	else
	  {
	    grid->xsize = 0;
	    grid->ysize = 0;
	  }

	break;
      }
    default:
      {
	Error("Unsupported grid type: %s", gridNamePtr(gridtype));
	break;
      }
    }

  grid->isRotated = FALSE;
  if ( gribapiGetIsRotated(gh) )
    {
      grid->isRotated = TRUE;
      GRIB_CHECK(grib_get_double(gh, "latitudeOfSouthernPoleInDegrees",  &grid->ypole), 0);
      GRIB_CHECK(grib_get_double(gh, "longitudeOfSouthernPoleInDegrees", &grid->xpole), 0);
      GRIB_CHECK(grib_get_double(gh, "angleOfRotation", &grid->angle), 0);
      /* change from south to north pole */
      grid->ypole = -grid->ypole;
      grid->xpole =  grid->xpole - 180;
    }

  grid->xvals = NULL;
  grid->yvals = NULL;
  grid->type  = gridtype;
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
double grib1GetLevel(grib_handle *gh, int leveltype)
{
  double dlevel;

  GRIB_CHECK(grib_get_double(gh, "level", &dlevel), 0);
  if ( leveltype == 100 ) dlevel *= 100;
  if ( dlevel < -2.e9 || dlevel > 2.e9 ) dlevel = 0;

  return (dlevel);
}

static
double grib2GetLevel(grib_handle *gh, int leveltype)
{
  double dlevel;

  GRIB_CHECK(grib_get_double(gh, "level", &dlevel), 0);
  if ( leveltype == 100 ) dlevel *= 100;
  if ( dlevel < -2.e9 || dlevel > 2.e9 ) dlevel = 0;

  return (dlevel);
}

static
void gribapiAddRecord(int streamID, int param, grib_handle *gh,
		      long recsize, off_t position, int prec, int comptype)
{
  long editionNumber;
  int zaxistype;
  int gridID = CDI_UNDEFID, varID;
  int levelID = 0;
  int tsID, recID;
  int level1, level2;
  int numavg;
  int tsteptype;
  int lbounds = 0;
  record_t *record;
  grid_t grid;
  int vlistID;
  stream_t *streamptr;
  int leveltype;
  double dlevel;
  long lpar;
  int status;
  char name[256], longname[256], units[256];
  size_t vlen; 

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  tsID    = streamptr->curTsID;
  recID   = recordNewEntry(streamID, tsID);
  record  = &streamptr->tsteps[tsID].records[recID];

  tsteptype = gribapiGetTsteptype(gh);
  // numavg  = ISEC1_AvgNum;
  numavg  = 0;

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  if ( editionNumber <= 1 )
    {
      status = grib_get_long(gh, "indicatorOfTypeOfLevel", &lpar);
      if ( status == 0 )
	{
	  leveltype = (int) lpar;
	  dlevel = grib1GetLevel(gh, leveltype);
	  if ( leveltype == 99 ) leveltype++;
	}
      else 
	{
	  leveltype = 0;
	  dlevel = 0;
	}
    }
  else
    {
      status = grib_get_long(gh, "typeOfFirstFixedSurface", &lpar);
      if ( status == 0 )
	{
	  leveltype = (int) lpar;
	  dlevel = grib2GetLevel(gh, leveltype);
	  if ( leveltype == 99 ) leveltype++;
	}
      else 
	{
	  leveltype = 0;
	  dlevel = 0;
	}
    }

  level1 = (int) dlevel;
  level2 = 0;

  // fprintf(stderr, "param %d %d %d %d\n", param, level1, level2, leveltype);

  (*record).size     = recsize;
  (*record).position = position;
  (*record).param    = param;
  (*record).ilevel   = level1;
  (*record).ilevel2  = level2;
  (*record).ltype    = leveltype;

  gribapiGetGrid(gh, &grid);

  gridID = varDefGrid(vlistID, grid, 0);

  zaxistype = gribapiGetZaxisType(editionNumber, leveltype);

  if ( zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF )
    {
      int vctsize;
      size_t dummy;
      double *vctptr;

      GRIB_CHECK(grib_get_long(gh, "NV", &lpar), 0);
      vctsize = lpar;
      if ( vctsize > 0 )
	{
	  vctptr = (double *) malloc(vctsize*sizeof(double));
	  dummy = vctsize;
	  GRIB_CHECK(grib_get_double_array(gh, "pv", vctptr, &dummy), 0);
	  varDefVCT(vctsize, vctptr);
	  free(vctptr);
	}
    }

  //lbounds = cgribexGetZaxisHasBounds(ISEC1_LevelType);

  if ( prec > 32 ) prec = DATATYPE_PACK32;
  if ( prec <  0 ) prec = DATATYPE_PACK;

  name[0] = 0;
  longname[0] = 0;
  units[0] = 0;

  vlen = 256;
  GRIB_CHECK(grib_get_string(gh, "shortName", name, &vlen), 0);
  if      ( vlen == 8 && memcmp(name, "unknown", vlen) == 0 ) name[0] = 0;
  else if ( vlen == 2 && memcmp(name, "~", vlen)       == 0 ) name[0] = 0;

  if ( name[0] != 0 )
    {
      vlen = 256;
      GRIB_CHECK(grib_get_string(gh, "name", longname, &vlen), 0);
      if ( vlen == 8 && memcmp(longname, "unknown", vlen) == 0 ) longname[0] = 0;
      vlen = 256;
      GRIB_CHECK(grib_get_string(gh, "units", units, &vlen), 0);
      if ( vlen == 8 && memcmp(units, "unknown", vlen) == 0 ) units[0] = 0;
    }
  // fprintf(stderr, "param %d name %s %s %s\n", param, name, longname, units); 

  varAddRecord(recID, param, gridID, zaxistype, lbounds, level1, level2,
	       prec, &varID, &levelID, 0, numavg, leveltype,
	       name, longname, units);

  (*record).varID   = varID;
  (*record).levelID = levelID;

  varDefCompType(varID, comptype);

  if ( varInqInst(varID) == CDI_UNDEFID )
    {
      long center, subcenter;
      int instID;
      GRIB_CHECK(grib_get_long(gh, "centre", &center), 0);
      GRIB_CHECK(grib_get_long(gh, "subCentre", &subcenter), 0);
      instID    = institutInq((int)center, (int)subcenter, NULL, NULL);
      if ( instID == CDI_UNDEFID )
	instID = institutDef((int)center, (int)subcenter, NULL, NULL);
      varDefInst(varID, instID);
    }

  if ( varInqModel(varID) == CDI_UNDEFID )
    {
      int modelID;
      long processID;
      status = grib_get_long(gh, "generatingProcessIdentifier", &processID);
      if ( status == 0 )
	{
	  modelID = modelInq(varInqInst(varID), processID, NULL);
	  if ( modelID == CDI_UNDEFID )
	    modelID = modelDef(varInqInst(varID), processID, NULL);
	  varDefModel(varID, modelID);
	}
    }

  if ( varInqTable(varID) == CDI_UNDEFID )
    {
      int pdis, pcat, pnum;

      cdiDecodeParam(param, &pnum, &pcat, &pdis);

      if ( pdis == 255 )
	{
	  int tableID;
	  int tabnum = pcat;

	  tableID = tableInq(varInqModel(varID), tabnum, NULL);

	  if ( tableID == CDI_UNDEFID )
	    tableID = tableDef(varInqModel(varID), tabnum, NULL);
	  varDefTable(varID, tableID);
	}
    }

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;

  if ( CDI_Debug )
    Message("varID = %d  param = %d  zaxistype = %d  gridID = %d  levelID = %d",
	    varID, param, zaxistype, gridID, levelID);
}
#endif

int gribapiScanTimestep1(int streamID)
{
#if  defined  (HAVE_LIBGRIB_API)
  off_t recpos = 0;
  unsigned char *gribbuffer = NULL;
  long buffersize = 0;
  int iret = 0, ipunp = 0, iword = 0;
  int rstatus;
  int status;
  int fileID;
  int rtabnum = 0;
  int rcode = 0, level1 = 0, level2 = 0;
  int vdate = 0, vtime = 0;
  int param = 0;
  DateTime datetime, datetime0;
  int tsID;
  int varID;
  size_t readsize;
  int nrecords, nrecs, recID;
  int prec;
  long recsize = 0;
  int warn_time = TRUE;
  int warn_numavg = TRUE;
  int taxisID = -1;
  int rdate = 0, rtime = 0, tunit = 0, fcast = 0;
  taxis_t *taxis;
  int vlistID;
  int comptype;
  long unzipsize;
  compvar2_t compVar, compVar0;
  stream_t *streamptr;
  grib_handle *gh = NULL;
  int leveltype;
  int pdis = 0, pcat = 0, pnum = 0;
  long editionNumber;
  long lpar;
  int bitsPerValue;
  double dlevel = 0;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  streamptr->curTsID = 0;

  tsID  = tstepsNewEntry(streamID);
  taxis = &streamptr->tsteps[tsID].taxis;

  if ( tsID != 0 )
    Error("Internal problem! tstepsNewEntry returns %d", tsID);

  fileID = streamInqFileID(streamID);

  nrecs = 0;
  while ( TRUE )
    {
      recsize = gribGetSize(fileID);
      recpos  = fileGetPos(fileID);

      if ( recsize == 0 )
	{
	  streamptr->ntsteps = 1;
	  break;
	}
      if ( recsize > buffersize )
	{
	  buffersize = recsize;
	  gribbuffer = (unsigned char *) realloc(gribbuffer, buffersize);
	}

      readsize = recsize;
      rstatus = gribRead(fileID, gribbuffer, &readsize);
      if ( rstatus ) break;

      comptype = COMPRESS_NONE;
      if ( gribGetZip(recsize, gribbuffer, &unzipsize) > 0 )
	{
	  comptype = COMPRESS_SZIP;
	  unzipsize += 100;
	  if ( (long) buffersize < unzipsize )
	    {
	      buffersize = unzipsize;
	      gribbuffer = (unsigned char *) realloc(gribbuffer, buffersize);
	    }
	}

      gh = grib_handle_new_from_message(NULL, (void *) gribbuffer, recsize);
      GRIB_CHECK(grib_set_double(gh, "missingValue", GRIBAPI_MISSVAL), 0);

      GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

      if ( editionNumber <= 1 )
	{
	  GRIB_CHECK(grib_get_long(gh, "table2Version", &lpar), 0);
	  rtabnum = (int) lpar;
	  GRIB_CHECK(grib_get_long(gh, "indicatorOfParameter", &lpar), 0);
	  rcode = (int) lpar;

	  param = cdiEncodeParam(rcode, rtabnum, 255);

	  status = grib_get_long(gh, "indicatorOfTypeOfLevel", &lpar);
	  if ( status == 0 )
	    {
	      leveltype = (int) lpar;
	      dlevel = grib1GetLevel(gh, leveltype);
	      if ( leveltype == 99 ) leveltype++;
	    }
	  else 
	    {
	      leveltype = 0;
	      dlevel = 0;
	    }
	}
      else
	{
	  size_t len = 256;
	  char typeOfPacking[256];

	  status = grib_get_string(gh, "packingType", typeOfPacking, &len);
	  if ( status == 0 )
	    {
	      // fprintf(stderr, "packingType %d %s\n", len, typeOfPacking);
	      if ( strncmp(typeOfPacking, "grid_jpeg", len) == 0 ) comptype = COMPRESS_JPEG;
	    }
	  
	  GRIB_CHECK(grib_get_long(gh, "discipline", &lpar), 0);
	  pdis = (int) lpar;

	  GRIB_CHECK(grib_get_long(gh, "parameterCategory", &lpar), 0);
	  pcat = (int) lpar;

	  GRIB_CHECK(grib_get_long(gh, "parameterNumber", &lpar), 0);
	  pnum = (int) lpar;

	  param = cdiEncodeParam(pnum, pcat, pdis);

	  status = grib_get_long(gh, "typeOfFirstFixedSurface", &lpar);
	  if ( status == 0 )
	    {
	      leveltype = (int) lpar;
	      dlevel = grib2GetLevel(gh, leveltype);
	      if ( leveltype == 99 ) leveltype++;
	    }
	  else 
	    {
	      leveltype = 0;
	      dlevel = 0;
	    }
	}

      level1 = (int) dlevel;
      level2 = 0;

      gribapiGetValidityDateTime(gh, &vdate, &vtime);
      /*
      printf("%d %d %d.%d.%d %d %g\n", vdate, vtime, pnum, pcat, pdis, leveltype, dlevel);
      */
      prec = DATATYPE_PACK;
      status = grib_get_long(gh,"bitsPerValue", &lpar);
      if ( status == 0 )
	{
	  bitsPerValue = (int) lpar;
	  if ( bitsPerValue > 0 && bitsPerValue <= 32 )
	    prec = bitsPerValue;
	}

      if ( nrecs == 0 )
	{
	  datetime0.date = vdate;
	  datetime0.time = vtime;
	  GRIB_CHECK(grib_get_long(gh, "dataDate", &lpar), 0);
	  rdate = (int) lpar;
	  GRIB_CHECK(grib_get_long(gh, "dataTime", &lpar), 0);
	  rtime = (int) lpar*100;
	  fcast = gribapiTimeIsFC(gh);
	  if ( fcast ) tunit = gribapiGetTimeUnits(gh);
	}
      else
	{
	  datetime.date  = vdate;
	  datetime.time  = vtime;

	  compVar.param  = param;
          compVar.level1 = level1;
          compVar.level2 = level2;
	  compVar.ltype  = leveltype;

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      compVar0.param  = streamptr->tsteps[0].records[recID].param;
	      compVar0.level1 = streamptr->tsteps[0].records[recID].ilevel;
	      compVar0.level2 = streamptr->tsteps[0].records[recID].ilevel2;
	      compVar0.ltype  = streamptr->tsteps[0].records[recID].ltype;
	      /*
	      printf("var0: %d %d %d %d %d\n", recID, compVar0.param, compVar0.level1, compVar0.level2, compVar0.ltype);
	      printf("var1: %d %d %d %d %d\n", recID, compVar.param, compVar.level1, compVar.level2, compVar.ltype);
	      */
	      if ( memcmp(&compVar0, &compVar, sizeof(compvar2_t)) == 0 ) break;
	    }

	  if ( cdiInventoryMode == 1 )
	    {
	      if ( recID < nrecs ) break;
	      if ( warn_time )
		if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 )
		  {
		    char paramstr[32];
		    cdiParamToString(param, paramstr, sizeof(paramstr));	    
		    Warning("Inconsistent verification time (param=%s level=%d)", paramstr, level1);
		    warn_time = FALSE;
		  }
	    }
	  else
	    {
	      if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 ) break;

	      if ( recID < nrecs )
		{
		  char paramstr[32];
		  cdiParamToString(param, paramstr, sizeof(paramstr));
		  Warning("Param=%s level=%d already exist, skipped!", paramstr, level1);
		  continue;
		}
	    }
	}
      /*
      if ( ISEC1_AvgNum )
	{
	  if (  taxis->numavg && warn_numavg && (taxis->numavg != ISEC1_AvgNum) )
	    {
	      Message("Change numavg from %d to %d not allowed!",
		      taxis->numavg, ISEC1_AvgNum);
	      warn_numavg = FALSE;
	    }
	  else
	    {
	      taxis->numavg = ISEC1_AvgNum;
	    }
	}
      */
      nrecs++;

      if ( CDI_Debug )
	Message("%4d %8d %4d  %8d %8d %6d", nrecs, (int)recpos, param, level1, vdate, vtime);

      gribapiAddRecord(streamID, param, gh, recsize, recpos, prec, comptype);

      grib_handle_delete(gh);
      gh = NULL;
    }

  if ( gh ) grib_handle_delete(gh);

  streamptr->rtsteps = 1;

  if ( nrecs == 0 ) return (CDI_EUFSTRUCT);

  cdiGenVars(streamID);

  if ( fcast )
    {
      taxisID = taxisCreate(TAXIS_RELATIVE);
      taxis->type  = TAXIS_RELATIVE;
      taxis->rdate = rdate;
      taxis->rtime = rtime;
      taxis->unit  = tunit;
    }
  else
    {
      taxisID = taxisCreate(TAXIS_ABSOLUTE);
      taxis->type  = TAXIS_ABSOLUTE;
    }

  taxis->vdate = datetime0.date;
  taxis->vtime = datetime0.time;

  vlistID = streamInqVlist(streamID);
  vlistDefTaxis(vlistID, taxisID);

  nrecords = streamptr->tsteps[0].nallrecs;
  if ( nrecords < streamptr->tsteps[0].recordSize )
    {
      streamptr->tsteps[0].recordSize = nrecords;
      streamptr->tsteps[0].records =
      (record_t *) realloc(streamptr->tsteps[0].records, nrecords*sizeof(record_t));
    }

  streamptr->tsteps[0].recIDs = (int *) malloc(nrecords*sizeof(int));
  streamptr->tsteps[0].nrecs = nrecords;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[0].recIDs[recID] = recID;

  streamptr->record->buffer     = gribbuffer;
  streamptr->record->buffersize = buffersize;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamID);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = TRUE;
      streamptr->tsteps[tsID].position = recpos;
    }

  if ( streamptr->ntsteps == 1 )
    {
      if ( taxis->vdate == 0 && taxis->vtime == 0 )
	{
	  streamptr->ntsteps = 0;
	  for ( varID = 0; varID < streamptr->nvars; varID++ )
	    {
	      vlistDefVarTime(vlistID, varID, TIME_CONSTANT);
	    }
	}
    }
#else
  Error("GRIB_API support not compiled in!");
#endif

  return (0);
}


int gribapiScanTimestep2(int streamID)
{
  int rstatus = 0;
#if  defined  (HAVE_LIBGRIB_API)
  int status;
  off_t recpos = 0;
  unsigned char *gribbuffer = NULL;
  long buffersize = 0;
  int iret = 0, ipunp = 0, iword = 0;
  int fileID;
  int rtabnum = 0;
  int rcode = 0, level1 = 0, level2 = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0;
  int tsID;
  int varID, gridID;
  size_t readsize;
  int nrecords, nrecs, recID, rindex;
  long recsize = 0;
  int warn_numavg = TRUE;
  int tsteptype;
  int taxisID = -1;
  taxis_t *taxis;
  int vlistID;
  long unzipsize;
  compvar2_t compVar, compVar0;
  stream_t *streamptr;
  grib_handle *gh = NULL;
  int leveltype;
  int pdis = 0, pcat = 0, pnum = 0;
  int param = 0;
  long editionNumber;
  long lpar;
  double dlevel = 0;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  streamptr->curTsID = 1;

  fileID  = streamInqFileID(streamID);
  vlistID = streamInqVlist(streamID);
  taxisID = vlistInqTaxis(vlistID);

  gribbuffer = (unsigned char *) streamptr->record->buffer;
  buffersize = streamptr->record->buffersize;
                                          
  tsID = streamptr->rtsteps;
  if ( tsID != 1 )
    Error("Internal problem! unexpeceted timestep %d", tsID+1);

  taxis = &streamptr->tsteps[tsID].taxis;

  fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

  cdiCreateRecords(streamID, tsID);

  nrecords = streamptr->tsteps[tsID].nallrecs;
  streamptr->tsteps[1].recIDs = (int *) malloc(nrecords*sizeof(int));
  streamptr->tsteps[1].nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[1].recIDs[recID] = -1;
      
  for ( recID = 0; recID < nrecords; recID++ )
    {
      varID = streamptr->tsteps[0].records[recID].varID;
      streamptr->tsteps[tsID].records[recID].position = 
	streamptr->tsteps[0].records[recID].position;
      streamptr->tsteps[tsID].records[recID].size     = 
	streamptr->tsteps[0].records[recID].size;
    }

  rindex = 0;
  while ( TRUE )
    {
      if ( rindex > nrecords ) break;

      recsize = gribGetSize(fileID);
      recpos  = fileGetPos(fileID);
      if ( recsize == 0 )
	{
	  streamptr->ntsteps = 2;
	  break;
	}
      if ( recsize > buffersize )
	{
	  buffersize = recsize;
	  gribbuffer = (unsigned char *) realloc(gribbuffer, (size_t)buffersize);
	}

      readsize = recsize;
      rstatus = gribRead(fileID, gribbuffer, &readsize);
      if ( rstatus ) break;

      if ( gribGetZip(recsize, gribbuffer, &unzipsize) > 0 )
	{
	  unzipsize += 100; /* need 0 to 1 bytes for rounding of bds */
	  if ( (long) buffersize < unzipsize )
	    {
	      buffersize = unzipsize;
	      gribbuffer = (unsigned char *) realloc(gribbuffer, buffersize);
	    }
	}

      gh = grib_handle_new_from_message(NULL, (void *) gribbuffer, recsize);
      GRIB_CHECK(grib_set_double(gh, "missingValue", GRIBAPI_MISSVAL), 0);

      GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

      if ( editionNumber <= 1 )
	{
	  GRIB_CHECK(grib_get_long(gh, "table2Version", &lpar), 0);
	  rtabnum = (int) lpar;
	  GRIB_CHECK(grib_get_long(gh, "indicatorOfParameter", &lpar), 0);
	  rcode = (int) lpar;

	  param = cdiEncodeParam(rcode, rtabnum, 255);

	  status = grib_get_long(gh, "indicatorOfTypeOfLevel", &lpar);
	  if ( status == 0 )
	    {
	      leveltype = (int) lpar;
	      dlevel = grib1GetLevel(gh, leveltype);
	      if ( leveltype == 99 ) leveltype++;
	    }
	  else 
	    {
	      leveltype = 0;
	      dlevel = 0;
	    }
	}
      else
	{
	  GRIB_CHECK(grib_get_long(gh, "discipline", &lpar), 0);
	  pdis = (int) lpar;

	  GRIB_CHECK(grib_get_long(gh, "parameterCategory", &lpar), 0);
	  pcat = (int) lpar;

	  GRIB_CHECK(grib_get_long(gh, "parameterNumber", &lpar), 0);
	  pnum = (int) lpar;

	  param = cdiEncodeParam(pnum, pcat, pdis);

	  status = grib_get_long(gh, "typeOfFirstFixedSurface", &lpar);
	  if ( status == 0 )
	    {
	      leveltype = (int) lpar;
	      dlevel = grib2GetLevel(gh, leveltype);
	      if ( leveltype == 99 ) leveltype++;
	    }
	  else 
	    {
	      leveltype = 0;
	      dlevel = 0;
	    }
	}

      level1 = (int) dlevel;
      level2 = 0;

      gribapiGetValidityDateTime(gh, &vdate, &vtime);

      if ( rindex == 0 )
	{
	  if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
	    {
	      taxis->type  = TAXIS_RELATIVE;
	      GRIB_CHECK(grib_get_long(gh, "dataDate", &lpar), 0);
	      taxis->rdate = (int) lpar;
	      GRIB_CHECK(grib_get_long(gh, "dataTime", &lpar), 0);
	      taxis->rtime = (int) lpar*100;
	      taxis->unit  = gribapiGetTimeUnits(gh);
	    }
	  else
	    {
	      taxis->type  = TAXIS_ABSOLUTE;
	    }
	  taxis->vdate = vdate;
	  taxis->vtime = vtime;

	  datetime0.date = vdate;
	  datetime0.time = vtime;
	}

      tsteptype = gribapiGetTsteptype(gh);
      /*
      if ( ISEC1_AvgNum )
	{
	  if (  taxis->numavg && warn_numavg &&
		(taxis->numavg != ISEC1_AvgNum) )
	    {
	      warn_numavg = FALSE;
	    }
	  else
	    {
	      taxis->numavg = ISEC1_AvgNum;
	    }
	}
      */
      datetime.date  = vdate;
      datetime.time  = vtime;
      compVar.param  = param;
      compVar.level1 = level1;
      compVar.level2 = level2;
      compVar.ltype  = leveltype;
      for ( recID = 0; recID < nrecords; recID++ )
	{
	  compVar0.param = streamptr->tsteps[tsID].records[recID].param;
	  compVar0.level1 = streamptr->tsteps[tsID].records[recID].ilevel;
	  compVar0.level2 = streamptr->tsteps[tsID].records[recID].ilevel2;
	  compVar0.ltype  = streamptr->tsteps[tsID].records[recID].ltype;

	  if ( memcmp(&compVar0, &compVar, sizeof(compvar2_t)) == 0 ) break;
	}

      if ( recID == nrecords )
	{
	  char paramstr[32];
	  cdiParamToString(param, paramstr, sizeof(paramstr));
	  Warning("Param=%s level=%d not defined at timestep 1!", paramstr, level1);
	  return (CDI_EUFSTRUCT);
	}

      if ( cdiInventoryMode == 1 )
	{
	  if ( streamptr->tsteps[tsID].records[recID].used )
	    {
	      break;
	    }
	  else
	    {
	      streamptr->tsteps[tsID].records[recID].used = TRUE;
	      streamptr->tsteps[tsID].recIDs[rindex] = recID;
	    }    
	}
      else
	{
	  if ( streamptr->tsteps[tsID].records[recID].used )
	    {
	      char paramstr[32];
	      cdiParamToString(param, paramstr, sizeof(paramstr));

	      if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 ) break;

	      Warning("Param=%s level=%d already exist, skipped!", paramstr, level1);
	      continue;
	    }
	  else
	    {
	      streamptr->tsteps[tsID].records[recID].used = TRUE;
	      streamptr->tsteps[tsID].recIDs[rindex] = recID;
	    }    
	}

      if ( CDI_Debug )
	Message("%4d %8d %4d %8d %8d %6d", rindex+1, (int)recpos, param, level1, vdate, vtime);

      streamptr->tsteps[tsID].records[recID].size = recsize;

      compVar0.param  = streamptr->tsteps[tsID].records[recID].param;
      compVar0.level1 = streamptr->tsteps[tsID].records[recID].ilevel;
      compVar0.level2 = streamptr->tsteps[tsID].records[recID].ilevel2;
      compVar0.ltype  = streamptr->tsteps[tsID].records[recID].ltype;

      if ( memcmp(&compVar0, &compVar, sizeof(compvar2_t)) != 0 )
	{
	  Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		  tsID, recID,
		  streamptr->tsteps[tsID].records[recID].param, param,
		  streamptr->tsteps[tsID].records[recID].ilevel, level1);
	  return (CDI_EUFSTRUCT);
	}

      streamptr->tsteps[1].records[recID].position = recpos;
      varID = streamptr->tsteps[tsID].records[recID].varID;
      gridID = vlistInqVarGrid(vlistID, varID);
      /*
      if ( gridInqSize(gridID) == 1 && gridInqType(gridID) == GRID_LONLAT )
	{
	  if ( IS_NOT_EQUAL(gridInqXval(gridID, 0),ISEC2_FirstLon*0.001) ||
	       IS_NOT_EQUAL(gridInqYval(gridID, 0),ISEC2_FirstLat*0.001) )
	    gridChangeType(gridID, GRID_TRAJECTORY);
	}
      */
      if ( tsteptype != vlistInqVarTsteptype(vlistID, varID) )
	vlistDefVarTsteptype(vlistID, varID, tsteptype);

      grib_handle_delete(gh);
      gh = NULL;
      
      rindex++;
    }
  
  if ( gh ) grib_handle_delete(gh);

  nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    {
      if ( ! streamptr->tsteps[tsID].records[recID].used )
	{
	  varID = streamptr->tsteps[tsID].records[recID].varID;
	  vlistDefVarTime(vlistID, varID, TIME_CONSTANT);
	}
      else
	{
	  nrecs++;
	}
    }
  streamptr->tsteps[tsID].nrecs = nrecs;

  streamptr->rtsteps = 2;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamID);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = TRUE;
      streamptr->tsteps[tsID].position = recpos;
    }

  streamptr->record->buffer     = gribbuffer;
  streamptr->record->buffersize = buffersize;
#endif

  return (rstatus);
}


int gribapiScanTimestep(int streamID)
{
  int rstatus = 0;
#if  defined  (HAVE_LIBGRIB_API)
  int status = 0;
  long recsize = 0;
  off_t recpos = 0;
  unsigned char *gribbuffer;
  long buffersize = 0;
  int iret = 0, ipunp = 0, iword = 0;
  int fileID;
  int rtabnum = 0;
  int rcode = 0, level1 = 0, level2 = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0;
  int tsID;
  int vrecID, recID;
  int warn_numavg = TRUE;
  size_t readsize;
  int taxisID = -1;
  taxis_t *taxis;
  int vlistID;
  int rindex, nrecs = 0;
  long unzipsize;
  compvar2_t compVar, compVar0;
  stream_t *streamptr;
  grib_handle *gh = NULL;
  int leveltype;
  int pdis = 0, pcat = 0, pnum = 0;
  int param = 0;
  long editionNumber;
  long lpar;
  double dlevel = 0;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  vlistID = streamInqVlist(streamID);

  if ( CDI_Debug )
    {
      Message("streamID = %d", streamID);
      Message("cts = %d", streamptr->curTsID);
      Message("rts = %d", streamptr->rtsteps);
      Message("nts = %d", streamptr->ntsteps);
    }

  tsID  = streamptr->rtsteps;
  taxis = &streamptr->tsteps[tsID].taxis;

  if ( streamptr->tsteps[tsID].recordSize == 0 )
    {
      gribbuffer = (unsigned char *) streamptr->record->buffer;
      buffersize = streamptr->record->buffersize;

      cdiCreateRecords(streamID, tsID);

      nrecs = streamptr->tsteps[1].nrecs;

      streamptr->tsteps[tsID].nrecs = nrecs;
      streamptr->tsteps[tsID].recIDs = (int *) malloc(nrecs*sizeof(int));
      for ( recID = 0; recID < nrecs; recID++ )
	streamptr->tsteps[tsID].recIDs[recID] = streamptr->tsteps[1].recIDs[recID];

      fileID = streamInqFileID(streamID);

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

      rindex = 0;
      while ( TRUE )
	{
	  if ( rindex > nrecs ) break;

	  recsize = gribGetSize(fileID);
	  recpos  = fileGetPos(fileID);
	  if ( recsize == 0 )
	    {
	      streamptr->ntsteps = streamptr->rtsteps + 1;
	      break;
	    }

	  if ( rindex >= nrecs ) break;

	  if ( recsize > buffersize )
	    {
	      buffersize = recsize;
	      gribbuffer = (unsigned char *) realloc(gribbuffer, buffersize);
	    }

	  readsize = recsize;
	  rstatus = gribRead(fileID, gribbuffer, &readsize);
	  if ( rstatus )
	    {
	      Error("Inconsistent timestep %d (GRIB record %d/%d)!", tsID+1, rindex+1,
		    streamptr->tsteps[tsID].recordSize);
	      break;
	    }

	  if ( gribGetZip(recsize, gribbuffer, &unzipsize) > 0 )
	    {
	      unzipsize += 100; /* need 0 to 1 bytes for rounding of bds */
	      if ( (long) buffersize < unzipsize )
		{
		  buffersize = unzipsize;
		  gribbuffer = (unsigned char *) realloc(gribbuffer, buffersize);
		}
	    }

	  gh = grib_handle_new_from_message(NULL, (void *) gribbuffer, recsize);
	  GRIB_CHECK(grib_set_double(gh, "missingValue", GRIBAPI_MISSVAL), 0);

	  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

	  if ( editionNumber <= 1 )
	    {
	      GRIB_CHECK(grib_get_long(gh, "table2Version", &lpar), 0);
	      rtabnum = (int) lpar;
	      GRIB_CHECK(grib_get_long(gh, "indicatorOfParameter", &lpar), 0);
	      rcode = (int) lpar;

	      param = cdiEncodeParam(rcode, rtabnum, 255);

	      status = grib_get_long(gh, "indicatorOfTypeOfLevel", &lpar);
	      if ( status == 0 )
		{
		  leveltype = (int) lpar;
		  dlevel = grib1GetLevel(gh, leveltype);
		  if ( leveltype == 99 ) leveltype++;
		}
	      else 
		{
		  leveltype = 0;
		  dlevel = 0;
		}
	    }
	  else
	    {
	      GRIB_CHECK(grib_get_long(gh, "discipline", &lpar), 0);
	      pdis = (int) lpar;

	      GRIB_CHECK(grib_get_long(gh, "parameterCategory", &lpar), 0);
	      pcat = (int) lpar;

	      GRIB_CHECK(grib_get_long(gh, "parameterNumber", &lpar), 0);
	      pnum = (int) lpar;

	      param = cdiEncodeParam(pnum, pcat, pdis);

	      status = grib_get_long(gh, "typeOfFirstFixedSurface", &lpar);
	      if ( status == 0 )
		{
		  leveltype = (int) lpar;
		  dlevel = grib2GetLevel(gh, leveltype);
		  if ( leveltype == 99 ) leveltype++;
		}
	      else 
		{
		  leveltype = 0;
		  dlevel = 0;
		}
	    }

	  level1 = (int) dlevel;
	  level2 = 0;

	  gribapiGetValidityDateTime(gh, &vdate, &vtime);

	  if ( rindex == nrecs ) break;

	  if ( rindex == 0 )
	    {
	      taxisID = vlistInqTaxis(vlistID);
	      if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
		{
		  taxis->type  = TAXIS_RELATIVE;
		  GRIB_CHECK(grib_get_long(gh, "dataDate", &lpar), 0);
		  taxis->rdate = (int) lpar;
		  GRIB_CHECK(grib_get_long(gh, "dataTime", &lpar), 0);
		  taxis->rtime = (int) lpar*100;
		  taxis->unit  = gribapiGetTimeUnits(gh);
		}
	      else
		{
		  taxis->type  = TAXIS_ABSOLUTE;
		}
	      taxis->vdate = vdate;
	      taxis->vtime = vtime;

	      datetime0.date = vdate;
	      datetime0.time = vtime;
	    }
	  /*
	  if ( ISEC1_AvgNum )
	    {
	      if (  taxis->numavg && warn_numavg &&
		   (taxis->numavg != ISEC1_AvgNum) )
		{
		  warn_numavg = FALSE;
		}
	      else
		{
		  taxis->numavg = ISEC1_AvgNum;
		}
	    }
	  */
	  datetime.date  = vdate;
	  datetime.time  = vtime;
	  compVar.param  = param;
          compVar.level1 = level1;
          compVar.level2 = level2;
          compVar.ltype  = leveltype;
	  for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	    {
	      recID   = streamptr->tsteps[1].recIDs[vrecID];
	      compVar0.param  = streamptr->tsteps[tsID].records[recID].param;
	      compVar0.level1 = streamptr->tsteps[tsID].records[recID].ilevel;
	      compVar0.level2 = streamptr->tsteps[tsID].records[recID].ilevel2;
	      compVar0.ltype  = streamptr->tsteps[tsID].records[recID].ltype;

	      if ( memcmp(&compVar0, &compVar, sizeof(compvar2_t)) == 0 ) break;
	    }

	  if ( vrecID == nrecs )
	    {
	      char paramstr[32];
	      cdiParamToString(param, paramstr, sizeof(paramstr));
	      Warning("Param=%s level=%d not available at timestep %d!", paramstr, level1, tsID+1);

	      if ( cdiInventoryMode == 1 )
		return (CDI_EUFSTRUCT);
	      else
		continue;
	    }

	  if ( cdiInventoryMode == 1 )
	    {
	      streamptr->tsteps[tsID].records[recID].used = TRUE;
	      streamptr->tsteps[tsID].recIDs[rindex] = recID;
	    }
	  else
	    {
	      if ( streamptr->tsteps[tsID].records[recID].used )
		{
		  char paramstr[32];
		  cdiParamToString(param, paramstr, sizeof(paramstr));

		  if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 ) break;
		  
		  if ( CDI_Debug )
		    Warning("Param=%s level=%d already exist, skipped!", paramstr, level1);

		  continue;
		}
	      else
		{
		  streamptr->tsteps[tsID].records[recID].used = TRUE;
		  streamptr->tsteps[tsID].recIDs[rindex] = recID;
		}    
	    }

	  if ( CDI_Debug )
	    Message("%4d %8d %4d %8d %8d %6d", rindex+1, (int)recpos, param, level1, vdate, vtime);

	  compVar0.param  = streamptr->tsteps[tsID].records[recID].param;
	  compVar0.level1 = streamptr->tsteps[tsID].records[recID].ilevel;
	  compVar0.level2 = streamptr->tsteps[tsID].records[recID].ilevel2;
	  compVar0.ltype  = streamptr->tsteps[tsID].records[recID].ltype;

	  if ( memcmp(&compVar0, &compVar, sizeof(compvar2_t)) != 0 )
	    {
	      Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		      tsID, recID,
		      streamptr->tsteps[tsID].records[recID].param, param,
		      streamptr->tsteps[tsID].records[recID].ilevel, level1);
	      Error("Invalid, unsupported or inconsistent record structure");
	    }
	  
	  streamptr->tsteps[tsID].records[recID].position = recpos;
	  streamptr->tsteps[tsID].records[recID].size = recsize;

	  if ( CDI_Debug )
	    Message("%4d %8d %4d %8d %8d %6d", rindex, (int)recpos, param, level1, vdate, vtime);

	  grib_handle_delete(gh);
	  gh = NULL;

	  rindex++;
	}

      if ( gh ) grib_handle_delete(gh);

      for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	{
	  recID   = streamptr->tsteps[tsID].recIDs[vrecID];
	  if ( ! streamptr->tsteps[tsID].records[recID].used ) break;
	}

      if ( vrecID < nrecs )
	{
	  char paramstr[32];
	  cdiParamToString(streamptr->tsteps[tsID].records[recID].param, paramstr, sizeof(paramstr));
	  Warning("Param %d level %d not found at timestep %d!",
		  paramstr, streamptr->tsteps[tsID].records[recID].ilevel, tsID+1);
	  return (CDI_EUFSTRUCT);
	}

      streamptr->rtsteps++;

      if ( streamptr->ntsteps != streamptr->rtsteps )
	{
	  tsID = tstepsNewEntry(streamID);
	  if ( tsID != streamptr->rtsteps )
	    Error("Internal error. tsID = %d", tsID);

	  streamptr->tsteps[tsID-1].next   = 1;
	  streamptr->tsteps[tsID].position = recpos;
	}

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);
      streamptr->tsteps[tsID].position = recpos;

      streamptr->record->buffer     = gribbuffer;
      streamptr->record->buffersize = buffersize;
    }

  if ( nrecs > 0 && nrecs < streamptr->tsteps[tsID].nrecs )
    {
      Warning("Incomplete timestep. Stop scanning at timestep %d.", tsID);
      streamptr->ntsteps = tsID;
    }

  rstatus = streamptr->ntsteps;
#else
  Error("GRIB_API support not compiled in!");
#endif

  return (rstatus);
}


int gribapiDecode(unsigned char *gribbuffer, int gribsize, double *data, int gridsize,
		  int unreduced, int *nmiss, int *zip, double missval)
{
  int status = 0;
#if  defined  (HAVE_LIBGRIB_API)
  long lpar;
  long editionNumber, numberOfPoints;
  size_t datasize, dummy, recsize;
  grib_handle *gh = NULL;

  if ( unreduced )
    {
      static int lwarn = 1;

      if ( lwarn )
	{
	  lwarn = 0;
	  Warning("Conversion of gaussian reduced grids unsupported!");
	}
    }

  recsize = gribsize;
  gh = grib_handle_new_from_message(NULL, (void *) gribbuffer, recsize);
  GRIB_CHECK(grib_set_double(gh, "missingValue", missval), 0);

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  /* get the size of the values array*/
  GRIB_CHECK(grib_get_size(gh, "values", &datasize), 0);
  GRIB_CHECK(grib_get_long(gh, "numberOfPoints", &numberOfPoints), 0);

  // printf("values_size = %d  numberOfPoints = %ld\n", datasize, numberOfPoints);

  if ( gridsize != (int) datasize )
    Error("Internal problem: gridsize(%d) != datasize(%d)!", gridsize, datasize);
  dummy = datasize;
  GRIB_CHECK(grib_get_double_array(gh, "values", data, &dummy), 0);

  int gridtype;
  GRIB_CHECK(grib_get_long(gh, "gridDefinitionTemplateNumber", &lpar), 0);
  gridtype = (int) lpar;

  *nmiss = 0;
  if ( gridtype < 50 || gridtype > 53 )
    {
      GRIB_CHECK(grib_get_long(gh, "numberOfMissing", &lpar), 0);
      *nmiss = (int) lpar;
      // printf("gridtype %d, nmiss %d\n", gridtype, nmiss);
    }

  grib_handle_delete(gh);

#else
  Error("GRIB_API support not compiled in!");
#endif

  return (status);
}

#if  defined  (HAVE_LIBGRIB_API)
static
void gribapiDefInstitut(grib_handle *gh, int vlistID, int varID)
{
  int instID;

  if ( vlistInqInstitut(vlistID) != CDI_UNDEFID )
    instID = vlistInqInstitut(vlistID);
  else
    instID = vlistInqVarInstitut(vlistID, varID);

  if ( instID != CDI_UNDEFID )
    {
      long center, subcenter;
      long center0, subcenter0;

      center    = institutInqCenter(instID);
      subcenter = institutInqSubcenter(instID);

      GRIB_CHECK(grib_get_long(gh, "centre", &center0), 0);
      GRIB_CHECK(grib_get_long(gh, "subCentre", &subcenter0), 0);
      
      if ( center != center0 )
	GRIB_CHECK(grib_set_long(gh, "centre", center), 0);
      if ( subcenter != subcenter0 )
	GRIB_CHECK(grib_set_long(gh, "subCentre", subcenter), 0);
    }
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
void gribapiDefModel(grib_handle *gh, int vlistID, int varID)
{
  int modelID;

  if ( vlistInqModel(vlistID) != CDI_UNDEFID )
    modelID = vlistInqModel(vlistID);
  else
    modelID = vlistInqVarModel(vlistID, varID);

  if ( modelID != CDI_UNDEFID )
    GRIB_CHECK(grib_set_long(gh, "generatingProcessIdentifier", modelInqGribID(modelID)), 0);
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
void gribapiDefParam(grib_handle *gh, int param, const char *name)
{
  int pdis, pcat, pnum;
  long editionNumber;

  cdiDecodeParam(param, &pnum, &pcat, &pdis);

  if ( pnum < 0 )
    {
      size_t len;
      int status;
      len = strlen(name);
      status = grib_set_string(gh, "shortName", name, &len);
      if ( status != 0 )
	Warning("grib_api: No match for shortName=%s", name);
    }
  else
    {
      if ( pnum < 0 ) pnum = -pnum;

      GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

      if ( editionNumber <= 1 )
	{
	  static int lwarn = TRUE;

	  if ( pdis != 255 && lwarn )
	    {
	      char paramstr[32];
	      cdiParamToString(param, paramstr, sizeof(paramstr));
	      Warning("Can not convert GRIB2 parameter (%s) to GRIB1!", paramstr);
	      lwarn = FALSE;
	    }

	  GRIB_CHECK(grib_set_long(gh, "table2Version",        pcat), 0);
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfParameter", pnum), 0);
	}
      else
	{
	  GRIB_CHECK(grib_set_long(gh, "discipline",        pdis), 0);
	  GRIB_CHECK(grib_set_long(gh, "parameterCategory", pcat), 0);
	  GRIB_CHECK(grib_set_long(gh, "parameterNumber",   pnum), 0);
	}
    }

  // printf("param: %d.%d.%d %s\n", pnum, pcat, pdis, name);
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
int gribapiDefTimerange(int tsteptype, int factor, int calendar,
			int rdate, int rtime, int vdate, int vtime, int *pip)
{
  int timerange = -1;
  int year, month, day, hour, minute, second;
  int julday1, secofday1, julday2, secofday2, days, secs;
  int ip = 0;

  cdiDecodeDate(rdate, &year, &month, &day);
  cdiDecodeTime(rtime, &hour, &minute, &second);
  encode_juldaysec(calendar, year, month, day, hour, minute, &julday1, &secofday1);

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);
  encode_juldaysec(calendar, year, month, day, hour, minute, &julday2, &secofday2);
  
  (void) julday_sub(julday1, secofday1, julday2, secofday2, &days, &secs);

  if ( !(int) fmod(days*86400.0 + secs, factor) )
    {
      ip = (int) ((days*86400.0 + secs)/factor);

      switch ( tsteptype )
	{
	case TSTEP_AVG:    timerange =  0;  break;
	case TSTEP_ACCUM:  timerange =  1;  break;
	case TSTEP_MIN:    timerange =  2;  break;
	case TSTEP_MAX:    timerange =  3;  break;
	case TSTEP_DIFF:   timerange =  4;  break;
	default:           timerange =  0;  break;
	}
    }

  *pip = ip;

  return (timerange);
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
int gribapiDefDateTime(grib_handle *gh, int timeunit, int date, int time)
{
  int factor = 1;
  char stepunits[8];
  size_t len;

  if ( date == 0 ) date = 10101;
  GRIB_CHECK(grib_set_long(gh, "dataDate", date), 0);
  GRIB_CHECK(grib_set_long(gh, "dataTime", time/100), 0);

  switch (timeunit)
    {
    case TUNIT_SECOND:  factor =     1; strcpy(stepunits, "s");   break;
    case TUNIT_MINUTE:  factor =    60; strcpy(stepunits, "m");   break;
    case TUNIT_HOUR:    factor =  3600; strcpy(stepunits, "h");   break;
    case TUNIT_3HOURS:  factor = 10800; strcpy(stepunits, "3h");  break;
    case TUNIT_6HOURS:  factor = 21600; strcpy(stepunits, "6h");  break;
    case TUNIT_12HOURS: factor = 43200; strcpy(stepunits, "12h"); break;
    case TUNIT_DAY:     factor = 86400; strcpy(stepunits, "D");   break;
    default:            factor =  3600; strcpy(stepunits, "h");   break;
    }

  len = strlen(stepunits) + 1;
  GRIB_CHECK(grib_set_string(gh, "stepUnits", stepunits, &len), 0);  

  return (factor);
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
void gribapiDefTime(grib_handle *gh , int vdate, int vtime, int tsteptype, int numavg, int taxisID)
{
  int timetype = -1;
  int timerange = 0;
  int timeunit;

  if ( taxisID != -1 ) timetype = taxisInqType(taxisID);

  timeunit = taxisInqTunit(taxisID);

  if ( timetype == TAXIS_RELATIVE )
    {
      int factor = 1;
      int rdate, rtime;
      int ip = 0;
      int calendar;
      
      calendar = taxisInqCalendar(taxisID);
      rdate    = taxisInqRdate(taxisID);
      rtime    = taxisInqRtime(taxisID);

      factor = gribapiDefDateTime(gh, timeunit, rdate, rtime);

      timerange = gribapiDefTimerange(tsteptype, factor, calendar,
				      rdate, rtime, vdate, vtime, &ip);
      // printf("timerange: %d %d\n", timerange, ip);

      if ( ip > 0 )
	{
	  GRIB_CHECK(grib_set_long(gh, "significanceOfReferenceTime", 1), 0);
	  //GRIB_CHECK(grib_set_long(gh, "stepType", timerange), 0);
	  GRIB_CHECK(grib_set_long(gh, "stepRange", ip), 0);
	}
    }

  if ( timetype == TAXIS_ABSOLUTE )
    {
      (void) gribapiDefDateTime(gh, timeunit, vdate, vtime);

      // GRIB_CHECK(grib_set_long(gh, "stepType", timerange), 0);
    }
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
void gribapiDefGrid(grib_handle *gh, int gridID, int ljpeg)
{
  int gridtype;
  int status;
  long editionNumber;
  static short lwarn = TRUE;
  size_t len;
  char *mesg;

  gridtype = gridInqType(gridID);

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);
  if ( editionNumber <= 1 )
    if ( gridtype == GRID_GME || gridtype == GRID_REFERENCE )
      gridtype = -1;

  if ( gridtype == GRID_GENERIC )
    {
      int xsize, ysize, gridsize;

      gridsize = gridInqSize(gridID);
      xsize = gridInqXsize(gridID);
      ysize = gridInqYsize(gridID);

      if ( (ysize ==  32 || ysize ==  48 || ysize ==  64 ||
	    ysize ==  96 || ysize == 160 || ysize == 192 ||
	    ysize == 240 || ysize == 320 || ysize == 384 ||
	    ysize == 480 || ysize == 768 ) &&
	   (xsize == 2*ysize || xsize == 1) )
	{
	  gridtype = GRID_GAUSSIAN;
	  gridChangeType(gridID, gridtype);
	}
      else if ( gridsize == 1 )
	{
	  gridtype = GRID_LONLAT;
	  gridChangeType(gridID, gridtype);
	}
      else if ( gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL) )
	{
	  gridtype = GRID_LONLAT;
	  gridChangeType(gridID, gridtype);
	}
    }
  else if ( gridtype == GRID_CURVILINEAR )
    {
      if ( lwarn && gridInqSize(gridID) > 1 )
	{
	  lwarn = FALSE;
	  Warning("Curvilinear grids are unsupported in GRIB format! Created wrong GDS!");
	}
      gridtype = GRID_LONLAT;
    }

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_TRAJECTORY:
      {
	int nlon = 0, nlat;
	double xfirst = 0, xlast = 0, xinc = 0;
	double yfirst = 0, ylast = 0, yinc = 0;
	double latIncr;

	if ( gridtype == GRID_GAUSSIAN )
	  {
	    mesg = "regular_gg"; len = strlen(mesg);
	    GRIB_CHECK(grib_set_string(gh, "gridType", mesg, &len), 0);
	  }
	else if ( gridtype == GRID_GAUSSIAN_REDUCED )
	  {
	    mesg = "reduced_gg"; len = strlen(mesg);
	    GRIB_CHECK(grib_set_string(gh, "gridType", mesg, &len), 0);
	  }
	else if ( gridtype == GRID_LONLAT && gridIsRotated(gridID) )
	  {
	    mesg = "rotated_ll"; len = strlen(mesg);
	    GRIB_CHECK(grib_set_string(gh, "gridType", mesg, &len), 0);
	  }
	else
	  {
	    mesg = "regular_ll"; len = strlen(mesg);
	    GRIB_CHECK(grib_set_string(gh, "gridType", mesg, &len), 0);
	  }

	nlon = gridInqXsize(gridID);
	nlat = gridInqYsize(gridID);

	if ( gridtype == GRID_GAUSSIAN_REDUCED )
	  {
	    int *rowlon, i;
	    long *pl = NULL;

	    nlon = 0;

	    rowlon = (int *) malloc(nlat*sizeof(int));
	    pl     = (long *) malloc(nlat*sizeof(long));
	    gridInqRowlon(gridID, rowlon);
	    for ( i = 0; i < nlat; ++i ) pl[i] = rowlon[i];

	    // GRIB_CHECK(grib_set_long_array(gh, "pl", pl, nlat), 0);

	    free(pl);
	    free(rowlon);
	  }
	else
	  {
	    if ( nlon == 0 )
	      {
		nlon = 1;
	      }
	    else
	      {
		xfirst = gridInqXval(gridID,      0);
		xlast  = gridInqXval(gridID, nlon-1);
		xinc   = gridInqXinc(gridID);
	      }
	  }

	if ( nlat == 0 )
	  {
	    nlat = 1;
	  }
	else
	  {
	    yfirst = gridInqYval(gridID,      0);
	    ylast  = gridInqYval(gridID, nlat-1);
	    yinc   = gridInqYinc(gridID);
	  }

	GRIB_CHECK(grib_set_long(gh, "Ni", nlon), 0);
	GRIB_CHECK(grib_set_long(gh, "Nj", nlat), 0);
	GRIB_CHECK(grib_set_double(gh, "longitudeOfFirstGridPointInDegrees", xfirst), 0);
	GRIB_CHECK(grib_set_double(gh, "longitudeOfLastGridPointInDegrees",  xlast), 0);
	GRIB_CHECK(grib_set_double(gh, "latitudeOfFirstGridPointInDegrees",  yfirst), 0);
	GRIB_CHECK(grib_set_double(gh, "latitudeOfLastGridPointInDegrees",   ylast), 0);
	GRIB_CHECK(grib_set_double(gh, "iDirectionIncrementInDegrees", xinc), 0);

        {
          long jscan = 0;
          if ( yfirst < ylast ) jscan = 1;
          GRIB_CHECK(grib_set_long(gh, "jScansPositively", jscan), 0);
        }
	/*
	if ( fabs(xinc*1000 - ISEC2_LonIncr) > FLT_EPSILON )
	  ISEC2_LonIncr = 0;
	*/
	if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
	  GRIB_CHECK(grib_set_long(gh, "numberOfParallelsBetweenAPoleAndTheEquator", nlat/2), 0);
	else
	  {
	    latIncr = yinc;
	    if ( latIncr < 0 ) latIncr = -latIncr;
	    GRIB_CHECK(grib_set_double(gh, "jDirectionIncrementInDegrees", latIncr), 0);
	    /*
	    if ( fabs(yinc*1000 - ISEC2_LatIncr) > FLT_EPSILON )
	      ISEC2_LatIncr = 0;
	    */
	  }
	/*
	if ( ISEC2_NumLon > 1 && ISEC2_NumLat == 1 ) 
	  if ( ISEC2_LonIncr != 0 && ISEC2_LatIncr == 0 ) ISEC2_LatIncr = ISEC2_LonIncr;

	if ( ISEC2_NumLon == 1 && ISEC2_NumLat > 1 ) 
	  if ( ISEC2_LonIncr == 0 && ISEC2_LatIncr != 0 ) ISEC2_LonIncr = ISEC2_LatIncr;

	if ( ISEC2_LatIncr == 0 || ISEC2_LonIncr == 0 )
	  ISEC2_ResFlag = 0;
	else
	  ISEC2_ResFlag = 128;
	*/
	if ( gridIsRotated(gridID) )
	  {
	    double xpole, ypole, angle;
	    xpole = gridInqXpole(gridID);
	    ypole = gridInqYpole(gridID);
	    angle = gridInqAngle(gridID);
	    /* change from noth to south pole */
	    ypole = -ypole;
	    xpole =  xpole + 180;
	    GRIB_CHECK(grib_set_double(gh, "latitudeOfSouthernPoleInDegrees",  ypole), 0);
	    GRIB_CHECK(grib_set_double(gh, "longitudeOfSouthernPoleInDegrees", xpole), 0);
	    GRIB_CHECK(grib_set_double(gh, "angleOfRotation", angle), 0);
	  }

	/* East -> West */
	//if ( ISEC2_LastLon < ISEC2_FirstLon ) ISEC2_ScanFlag += 128;

	/* South -> North */
	//if ( ISEC2_LastLat > ISEC2_FirstLat ) ISEC2_ScanFlag += 64;

	if ( ljpeg && editionNumber == 2 )
	  {
	    mesg = "grid_jpeg"; len = strlen(mesg);
	    GRIB_CHECK(grib_set_string(gh, "packingType", mesg, &len), 0);
	  }
	else
	  {
	    mesg = "grid_simple"; len = strlen(mesg);
	    GRIB_CHECK(grib_set_string(gh, "packingType", mesg, &len), 0);
	  }

	break;
      }
      /*
    case GRID_LCC:
      {
	double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
	int xsize, ysize;
	int projflag, scanflag;

	xsize = gridInqXsize(gridID);
	ysize = gridInqYsize(gridID);

	gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
		   &projflag, &scanflag);

	ISEC2_GridType = GRIB2_GTYPE_LCC;
	ISEC2_NumLon   = xsize;
	ISEC2_NumLat   = ysize;
	ISEC2_FirstLon = NINT(originLon * 1000);
	ISEC2_FirstLat = NINT(originLat * 1000);
	ISEC2_Lambert_Lov    = NINT(lonParY * 1000);
	ISEC2_Lambert_LatS1  = NINT(lat1 * 1000);
	ISEC2_Lambert_LatS2  = NINT(lat2 * 1000);
	ISEC2_Lambert_dx     = NINT(xincm);
	ISEC2_Lambert_dy     = NINT(yincm);
	ISEC2_Lambert_LatSP  = 0;
	ISEC2_Lambert_LatSP  = 0;
	ISEC2_Lambert_ProjFlag = projflag;
	ISEC2_ScanFlag = scanflag;

	break;
      }
      */
    case GRID_SPECTRAL:
      {
	int trunc = gridInqTrunc(gridID);

	mesg = "sh"; len = strlen(mesg);
	GRIB_CHECK(grib_set_string(gh, "gridType", mesg, &len), 0);

	GRIB_CHECK(grib_set_long(gh, "J", trunc), 0);
	GRIB_CHECK(grib_set_long(gh, "K", trunc), 0);
	GRIB_CHECK(grib_set_long(gh, "M", trunc), 0);

	// GRIB_CHECK(grib_set_long(gh, "numberOfDataPoints", gridInqSize(gridID)), 0);

	if ( gridInqComplexPacking(gridID) )
	  {
	    if ( editionNumber == 2 ) GRIB_CHECK(grib_set_long(gh, "numberOfValues", gridInqSize(gridID)), 0);
	    mesg = "spectral_complex"; len = strlen(mesg);
	    GRIB_CHECK(grib_set_string(gh, "packingType", mesg, &len), 0);
	    /*
	    GRIB_CHECK(grib_set_long(gh, "JS", 20), 0);
	    GRIB_CHECK(grib_set_long(gh, "KS", 20), 0);
	    GRIB_CHECK(grib_set_long(gh, "MS", 20), 0);
	    */
	  }
	else
	  {
	    mesg = "spectral_simple"; len = strlen(mesg);
	    GRIB_CHECK(grib_set_string(gh, "packingType", mesg, &len), 0);
	  }

	break;
      }
    case GRID_GME:
      {
	GRIB_CHECK(grib_set_long(gh, "gridDefinitionTemplateNumber", GRIB2_GTYPE_GME), 0);

	GRIB_CHECK(grib_set_long(gh, "nd", gridInqGMEnd(gridID)), 0);
	GRIB_CHECK(grib_set_long(gh, "Ni", gridInqGMEni(gridID)), 0);
	GRIB_CHECK(grib_set_long(gh, "n2", gridInqGMEni2(gridID)), 0);
	GRIB_CHECK(grib_set_long(gh, "n3", gridInqGMEni3(gridID)), 0);
	GRIB_CHECK(grib_set_long(gh, "latitudeOfThePolePoint", 90000), 0);
	GRIB_CHECK(grib_set_long(gh, "longitudeOfThePolePoint", 0), 0);

	GRIB_CHECK(grib_set_long(gh, "numberOfDataPoints", gridInqSize(gridID)), 0);
	GRIB_CHECK(grib_set_long(gh, "totalNumberOfGridPoints", gridInqSize(gridID)), 0);

	break;
      }
    case GRID_REFERENCE:
      {
	static int warning = 1;
	status = grib_set_long(gh, "gridDefinitionTemplateNumber", GRIB2_GTYPE_NUMBER);
	if ( status != 0 && warning )
	  {
	    warning = 0;
	    Warning("Can not write reference grid!");
	    Warning("gridDefinitionTemplateNumber %d not found (grib2/template.3.%d.def)!",
		    GRIB2_GTYPE_NUMBER, GRIB2_GTYPE_NUMBER);
	  }
	else
	  {
	    GRIB_CHECK(grib_set_long(gh, "numberOfGridUsed", gridInqNumber(gridID)), 0);
	    GRIB_CHECK(grib_set_long(gh, "numberOfGridInReference", gridInqPosition(gridID)), 0);	
	  }

	break;
      }
    default:
      {
	Error("Unsupported grid type: %s", gridNamePtr(gridtype));
      }
    }
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
static
void gribapiDefLevel(grib_handle *gh, int param, int zaxisID, int levelID)
{
  double level;
  int zaxistype, ltype;
  static int warning = 1;
  long editionNumber;

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  zaxistype = zaxisInqType(zaxisID);
  ltype = zaxisInqLtype(zaxisID);

  if ( zaxistype == ZAXIS_GENERIC && ltype == 0 )
    {
      Message("Changed zaxis type from %s to %s",
	      zaxisNamePtr(zaxistype),
	      zaxisNamePtr(ZAXIS_PRESSURE));
      zaxistype = ZAXIS_PRESSURE;
      zaxisChangeType(zaxisID, zaxistype);
      zaxisDefUnits(zaxisID, "Pa");
    }

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:
      {
	if ( editionNumber <= 1 )
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_SURFACE), 0);
	else
	  GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_SURFACE), 0);
	GRIB_CHECK(grib_set_long(gh, "level", (long) zaxisInqLevel(zaxisID, levelID)), 0);
	break;
      }
    case ZAXIS_TOA:
      {
	if ( editionNumber <= 1 )
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_TOA), 0);
	else
	  GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_TOA), 0);
	break;
      }
    case ZAXIS_SEA_BOTTOM:
      {
	if ( editionNumber <= 1 )
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_SEA_BOTTOM), 0);
	else
	  GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_SEA_BOTTOM), 0);
	break;
      }
    case ZAXIS_ATMOSPHERE:
      {
	if ( editionNumber <= 1 )
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_ATMOSPHERE), 0);
	else
	  GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_ATMOSPHERE), 0);
	break;
      }
    case ZAXIS_MEANSEA:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	if ( editionNumber <= 1 )
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_MEANSEA), 0);
	else
	  GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_MEANSEA), 0);
	GRIB_CHECK(grib_set_double(gh, "level", level), 0);

	break;
      }
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
      {
	int vctsize;

	if ( editionNumber <= 1 )
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_HYBRID), 0);
	else
	  GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_HYBRID), 0);

	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	  {
	    GRIB_CHECK(grib_set_long(gh, "level", (long) zaxisInqLevel(zaxisID, levelID)), 0);
	    Error("hybrid_half model level code missing!");
	    /*
	    ISEC1_Level1    = (int) zaxisInqLbound(zaxisID, levelID);
	    ISEC1_Level2    = (int) zaxisInqUbound(zaxisID, levelID);
	    */
	  }
	else
	  {
	    GRIB_CHECK(grib_set_long(gh, "level", (long) zaxisInqLevel(zaxisID, levelID)), 0);
	  }

	vctsize = zaxisInqVctSize(zaxisID);
	if ( vctsize == 0 && warning )
	  {
	    char paramstr[32];
	    cdiParamToString(param, paramstr, sizeof(paramstr));
	    Warning("VCT missing ( param = %s, zaxisID = %d )", paramstr, zaxisID);
	    warning = 0;
	  }
	GRIB_CHECK(grib_set_long(gh, "PVPresent", 1), 0);
	GRIB_CHECK(grib_set_double_array(gh, "pv", zaxisInqVctPtr(zaxisID), vctsize), 0);

	break;
      }
    case ZAXIS_PRESSURE:
      {
	double dum;
	char units[128];

	level = zaxisInqLevel(zaxisID, levelID);
	if ( level < 0 )
	  Warning("Pressure level of %f Pa is below zero!", level);

	zaxisInqUnits(zaxisID, units);
	if ( memcmp(units, "Pa", 2) != 0 ) level *= 100;

	if ( level < 32768 && (level < 100 || modf(level/100, &dum) > 0) )
	  {
	    if ( editionNumber <= 1 )
	      GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_99), 0);
	    else
	      GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB1_LTYPE_99), 0);
	    GRIB_CHECK(grib_set_double(gh, "level", level), 0);
	  }
	else
	  {
	    if ( editionNumber <= 1 )
	      GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_ISOBARIC), 0);
	    else
	      GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_ISOBARIC), 0);
	    //GRIB_CHECK(grib_set_double(gh, "scaledValueOfFirstFixedSurface", level), 0);
	    GRIB_CHECK(grib_set_double(gh, "level", level/100), 0);
	  }

	break;
      }
    case ZAXIS_HEIGHT:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	if ( editionNumber <= 1 )
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_HEIGHT), 0);
	else
	  GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_HEIGHT), 0);
	GRIB_CHECK(grib_set_double(gh, "level", level), 0);

	break;
      }
    case ZAXIS_ALTITUDE:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	if ( editionNumber <= 1 )
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_ALTITUDE), 0);
	else
	  GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_ALTITUDE), 0);
	GRIB_CHECK(grib_set_double(gh, "level", level), 0);

	break;
      }
    case ZAXIS_SIGMA:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	if ( editionNumber <= 1 )
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_SIGMA), 0);
	else
	  GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_SIGMA), 0);
	GRIB_CHECK(grib_set_double(gh, "level", level), 0);

	break;
      }
    case ZAXIS_DEPTH_BELOW_LAND:
      {
	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	  {
	    GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_LANDDEPTH), 0);
	    GRIB_CHECK(grib_set_double(gh, "level", zaxisInqLbound(zaxisID, levelID)), 0);
	    GRIB_CHECK(grib_set_long(gh, "typeOfSecondFixedSurface", GRIB2_LTYPE_LANDDEPTH), 0);
	    GRIB_CHECK(grib_set_double(gh, "level", zaxisInqUbound(zaxisID, levelID)), 0);
	    /*
	    ISEC1_LevelType = GRIB2_LTYPE_LANDDEPTH_LAYER;
	    ISEC1_Level1    = (int) zaxisInqLbound(zaxisID, levelID);
	    ISEC1_Level2    = (int) zaxisInqUbound(zaxisID, levelID);
	    */
	  }
	else
	  {
	    level = zaxisInqLevel(zaxisID, levelID);

	    if ( editionNumber <= 1 )
	      GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_LANDDEPTH), 0);
	    else
	      GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_LANDDEPTH), 0);
	    GRIB_CHECK(grib_set_double(gh, "level", level), 0);
	  }

	break;
      }
    case ZAXIS_DEPTH_BELOW_SEA:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	if ( editionNumber <= 1 )
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_SEADEPTH), 0);
	else
	  GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_SEADEPTH), 0);
	GRIB_CHECK(grib_set_double(gh, "level", level), 0);

	break;
      }
    case ZAXIS_ISENTROPIC:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	if ( editionNumber <= 1 )
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", GRIB1_LTYPE_ISENTROPIC), 0);
	else
	  GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", GRIB2_LTYPE_ISENTROPIC), 0);
	GRIB_CHECK(grib_set_double(gh, "level", level), 0);

	break;
      }
    case ZAXIS_GENERIC:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	if ( editionNumber <= 1 )
	  GRIB_CHECK(grib_set_long(gh, "indicatorOfTypeOfLevel", ltype), 0);
	else
	  GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", ltype), 0);
	GRIB_CHECK(grib_set_double(gh, "level", level), 0);

	break;
      }
    default:
      {
	Error("Unsupported zaxis type: %s", zaxisNamePtr(zaxistype));
	break;
      }
    }
}
#endif


void *gribHandleNew(int editionNumber)
{
  void *gh = NULL;

#if  defined  (HAVE_LIBGRIB_API)
  if ( editionNumber == 1 )
    gh = (void *) grib_handle_new_from_samples(NULL, "GRIB1");
  else
    gh = (void *) grib_handle_new_from_samples(NULL, "GRIB2");

  if ( gh == NULL ) Error("grib_handle_new_from_samples failed!");
#endif

  return (gh);
}


void gribHandleDelete(void *gh)
{
#if  defined  (HAVE_LIBGRIB_API)
  grib_handle_delete(gh);
#endif
}

/* #define GRIBAPIENCODETEST 1 */

size_t gribapiEncode(int varID, int levelID, int vlistID, int gridID, int zaxisID,
		     int vdate, int vtime, int tsteptype, int numavg, 
		     long datasize, const double *data, int nmiss, unsigned char *gribbuffer, size_t gribbuffersize,
		     int ljpeg, void *gribContainer)
{
  size_t nbytes = 0;
#if  defined  (HAVE_LIBGRIB_API)
  size_t recsize = 0;
  void *dummy = NULL;
  int datatype;
  int param;
  long bitsPerValue;
  char name[256];
  grib_handle *gh = NULL;
  gribContainer_t *gc = (gribContainer_t *) gribContainer;
  // extern unsigned char _grib_template_GRIB2[];

  param    = vlistInqVarParam(vlistID, varID);
  datatype = vlistInqVarDatatype(vlistID, varID);
  vlistInqVarName(vlistID, varID, name);

#if defined(GRIBAPIENCODETEST)
  gh = (grib_handle *) gribHandleNew();
#else
  gh = gc->gribHandle;
#endif
  
  if ( ! gc->init ) gribapiDefInstitut(gh, vlistID, varID);
  if ( ! gc->init ) gribapiDefModel(gh, vlistID, varID);

  if ( ! gc->init ) gribapiDefParam(gh, param, name);
  gribapiDefTime(gh, vdate, vtime, tsteptype, numavg, vlistInqTaxis(vlistID));

  /* bitsPerValue have to be defined befor call to DefGrid (complex packing) */
  bitsPerValue = grbBitsPerValue(datatype);
  GRIB_CHECK(grib_set_long(gh, "bitsPerValue", bitsPerValue), 0);

  if ( ! gc->init ) gribapiDefGrid(gh, gridID, ljpeg);
  gribapiDefLevel(gh, param, zaxisID, levelID);

  if ( nmiss > 0 )
    {
      GRIB_CHECK(grib_set_long(gh, "bitmapPresent", 1), 0);
      GRIB_CHECK(grib_set_double(gh, "missingValue", vlistInqVarMissval(vlistID, varID)), 0);
    }

  GRIB_CHECK(grib_set_double_array(gh, "values", data, datasize), 0);

  /* get the size of coded message  */
  GRIB_CHECK(grib_get_message(gh, (const void **)&dummy, &recsize), 0);
  //recsize += 512; /* add some space for possible filling */

  if ( recsize > gribbuffersize )
    Error("Internal problem: GRIB buffer too small (size = %ld, needed = %ld", gribbuffersize, recsize);

  /* get a copy of the coded message */
  GRIB_CHECK(grib_get_message_copy(gh, gribbuffer, &recsize), 0);

#if defined(GRIBAPIENCODETEST)
  gribHandleDelete(gh);
#endif

  gc->init = TRUE;

  nbytes = recsize;
#else
  Error("GRIB_API support not compiled in!");
#endif

  return (nbytes);
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
