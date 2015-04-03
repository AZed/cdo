#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "file.h"
#include "varscan.h"
#include "datetime.h"
#include "vlist.h"
#include "stream_grb.h"
#include "calendar.h"


#if  defined  (HAVE_LIBGRIB_API)
#  include "cgribex.h"      /* gribGetSize, gribRead, gribGetZip, GRIB1_LTYPE_99 */
#  include "gribapi.h"
#  include "grib_api.h"
#endif

#define  NINT(x)  ((x) < 0 ? (int)((x)-.5) : (int)((x)+.5))

extern int cdiInventoryMode;

typedef struct {
  int param;
  int level1;
  int level2;
  int ltype;
  char name[32];
} compvar2_t;


#if  defined  (HAVE_LIBGRIB_API)
static
int gribapiGetGridType(grib_handle *gh)
{
  int gridtype = GRID_GENERIC;
  int gribgridtype = -1;
  long lpar;

    {
      int status;
      status = grib_get_long(gh, "gridDefinitionTemplateNumber", &lpar);

      if ( status ==  0 ) gribgridtype = (int) lpar;

      switch (gribgridtype)
	{
	case  GRIB2_GTYPE_LATLON:        { GRIB_CHECK(grib_get_long(gh, "Ni", &lpar), 0);
	                                   if ( lpar == (long) GRIB_MISSING_LONG ) break;
                                         }
	case  GRIB2_GTYPE_LATLON_ROT:    { gridtype = GRID_LONLAT;    break; }
	case  GRIB2_GTYPE_LCC:           { gridtype = GRID_LCC;       break; }
	case  GRIB2_GTYPE_GAUSSIAN:      { GRIB_CHECK(grib_get_long(gh, "Ni", &lpar), 0);
	                                   if ( lpar == (long) GRIB_MISSING_LONG )
					     gridtype = GRID_GAUSSIAN_REDUCED;
					   else
					     gridtype = GRID_GAUSSIAN;
				  	   break;
                                         }
	case  GRIB2_GTYPE_SPECTRAL:      { gridtype = GRID_SPECTRAL;  break; }
	case  GRIB2_GTYPE_GME:           { gridtype = GRID_GME;       break; }
	case  GRIB2_GTYPE_UNSTRUCTURED:  { gridtype = GRID_UNSTRUCTURED; break; }
	}
    }

  return (gridtype);
}

static
int gribapiGetIsRotated(grib_handle *gh)
{
  int isRotated = 0;
  int gribgridtype = -1;
  long lpar;
  int status;

  status = grib_get_long(gh, "gridDefinitionTemplateNumber", &lpar);

  if ( status ==  0 ) gribgridtype = (int) lpar;

  if ( gribgridtype == GRIB2_GTYPE_LATLON_ROT ) isRotated = 1;

  return (isRotated);
}

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

static
int getTimeunits(long unitsOfTime)
{
  int timeunits = -1;

  switch (unitsOfTime)
    {
    case 13:  timeunits = TUNIT_SECOND;  break;
    case  0:  timeunits = TUNIT_MINUTE;  break;
    case  1:  timeunits = TUNIT_HOUR;    break;
    case 10:  timeunits = TUNIT_3HOURS;  break;
    case 11:  timeunits = TUNIT_6HOURS;  break;
    case 12:  timeunits = TUNIT_12HOURS; break;
    case  2:  timeunits = TUNIT_DAY;     break;
    default:  timeunits = TUNIT_HOUR;    break;
    }

  return (timeunits);
}

static
double timeunit_factor(int tu1, int tu2)
{
  double factor = 1;

  if ( tu2 == TUNIT_HOUR )
    {
      switch (tu1)
        {
        case TUNIT_SECOND:  factor = 3600;   break;
        case TUNIT_MINUTE:  factor = 60;     break;
        case TUNIT_HOUR:    factor = 1;      break;
        case TUNIT_3HOURS:  factor = 1./3;   break;
        case TUNIT_6HOURS:  factor = 1./6;   break;
        case TUNIT_12HOURS: factor = 1./12;  break;
        case TUNIT_DAY:     factor = 1./24;  break;
        }
    }

  return (factor);
}

static
int gribapiGetTimeUnits(grib_handle *gh)
{
  int timeunits = -1;
  long unitsOfTime = -1;
  int status;
  // size_t len = 8;
  //char stepunits[8];
  //static int lprint = TRUE;

  status = grib_get_long(gh, "indicatorOfUnitOfTimeRange", &unitsOfTime);

  timeunits = getTimeunits(unitsOfTime);

  /*
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
  */

  return (timeunits);
}

static
int gribapiGetEndStep(grib_handle *gh, int startStep, int timeunits)
{
  int endStep = startStep;
  int timeunits2 = timeunits;
  int status;
  long unitsOfTime;
  long lpar;

  status = grib_get_long(gh, "stepUnits", &unitsOfTime);
  if ( status == 0 ) timeunits2 = getTimeunits(unitsOfTime);
  //timeunits2 = gribapiGetTimeUnits(gh);

  status = grib_get_long(gh, "endStep", &lpar);

  if ( status == 0 )
    endStep = (int) ((lpar * timeunit_factor(timeunits, timeunits2)) + 0.5);
  // printf("%d %d %d %d %d %g\n", startStep, endStep, lpar, timeunits, timeunits2, timeunit_factor(timeunits, timeunits2));

  return (endStep);
}

static
int gribapiTimeIsFC(grib_handle *gh)
{
  long editionNumber;
  int isFC = TRUE;

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  if ( editionNumber > 1 )
    {
      long sigofrtime;

      GRIB_CHECK(grib_get_long(gh, "significanceOfReferenceTime", &sigofrtime), 0);

      if ( sigofrtime == 3 ) isFC = FALSE;
    }

  return (isFC);
}

static
int gribapiGetTsteptype(grib_handle *gh)
{
  int tsteptype = TSTEP_INSTANT;
  static int lprint = TRUE;

  if ( gribapiTimeIsFC(gh) )
    {
      int status;
      size_t len = 256;
      char stepType[256];

      status = grib_get_string(gh, "stepType", stepType, &len);
      if ( status == 0 && len > 1 && len < 256 )
	{
	  if      ( strncmp("instant", stepType, len) == 0 ) tsteptype = TSTEP_INSTANT;
	  else if ( strncmp("avg",     stepType, len) == 0 ) tsteptype = TSTEP_AVG;
	  else if ( strncmp("accum",   stepType, len) == 0 ) tsteptype = TSTEP_ACCUM;
	  else if ( strncmp("max",     stepType, len) == 0 ) tsteptype = TSTEP_MAX;
	  else if ( strncmp("min",     stepType, len) == 0 ) tsteptype = TSTEP_MIN;
	  else if ( strncmp("diff",    stepType, len) == 0 ) tsteptype = TSTEP_DIFF;
	  else if ( strncmp("rms",     stepType, len) == 0 ) tsteptype = TSTEP_RMS;
	  else if ( strncmp("sd",      stepType, len) == 0 ) tsteptype = TSTEP_SD;
	  else if ( strncmp("cov",     stepType, len) == 0 ) tsteptype = TSTEP_COV;
	  else if ( strncmp("ratio",   stepType, len) == 0 ) tsteptype = TSTEP_RATIO;
	  else if ( lprint )
	    {
	      Message("stepType %s unsupported, set to instant!", stepType);
	      lprint = FALSE;
	    }

	  // printf("stepType: %s %ld %d\n", stepType, len, tsteptype);
	}
    }

  return (tsteptype);
}

static
void gribapiGetDataDateTime(grib_handle *gh, int *datadate, int *datatime)
{
  long lpar;

  GRIB_CHECK(grib_get_long(gh, "dataDate", &lpar), 0);
  *datadate = (int) lpar;
  GRIB_CHECK(grib_get_long(gh, "dataTime", &lpar), 0);
  *datatime = (int) lpar*100;
}

static
void gribapiSetDataDateTime(grib_handle *gh, int datadate, int datatime)
{
  GRIB_CHECK(grib_set_long(gh, "dataDate", datadate), 0);
  GRIB_CHECK(grib_set_long(gh, "dataTime", datatime/100), 0);
}

static
int gribapiGetValidityDateTime(grib_handle *gh, int *vdate, int *vtime)
{
  int rdate, rtime;
  int timeUnits, startStep = 0, endStep;
  int tstepRange = 0;
  int range;
  int status;
  long lpar;
  long sigofrtime = 3;
  long editionNumber;

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  if ( editionNumber > 1 )
    {
      GRIB_CHECK(grib_get_long(gh, "significanceOfReferenceTime", &sigofrtime), 0);
    }
  else
    {
      GRIB_CHECK(grib_get_long(gh, "timeRangeIndicator", &sigofrtime), 0);
    }

  if ( sigofrtime == 3 )
    {
      gribapiGetDataDateTime(gh, vdate, vtime);
    }
  else
    {
      gribapiGetDataDateTime(gh, &rdate, &rtime);

      status = grib_get_long(gh, "forecastTime", &lpar);
      if ( status == 0 ) startStep = (int) lpar;
      timeUnits = gribapiGetTimeUnits(gh);
      endStep = gribapiGetEndStep(gh, startStep, timeUnits);

      range = endStep - startStep;

      if ( range > 0 )
	{
	  if ( startStep == 0 ) tstepRange = -1;
	  else                  tstepRange =  1;
	}

      {
	static int lprint = TRUE;
	extern int grib_calendar;
	int ryear, rmonth, rday, rhour, rminute, rsecond;
	int julday, secofday;
	int64_t time_period = endStep;
        int64_t addsec;

	cdiDecodeDate(rdate, &ryear, &rmonth, &rday);
	cdiDecodeTime(rtime, &rhour, &rminute, &rsecond);

	encode_caldaysec(grib_calendar, ryear, rmonth, rday, rhour, rminute, rsecond, &julday, &secofday);

	addsec = 0;
	switch ( timeUnits )
	  {
	  case TUNIT_SECOND:  addsec =         time_period; break;
	  case TUNIT_MINUTE:  addsec =    60 * time_period; break;
	  case TUNIT_HOUR:    addsec =  3600 * time_period; break;
	  case TUNIT_3HOURS:  addsec = 10800 * time_period; break;
	  case TUNIT_6HOURS:  addsec = 21600 * time_period; break;
	  case TUNIT_12HOURS: addsec = 43200 * time_period; break;
	  case TUNIT_DAY:     addsec = 86400 * time_period; break;
	  default:
	    if ( lprint )
	      {
	        Warning("Time unit %d unsupported", timeUnits);
		lprint = FALSE;
	      }
	    break;
	  }

	julday_add_seconds(addsec, &julday, &secofday);

	decode_caldaysec(grib_calendar, julday, secofday, &ryear, &rmonth, &rday, &rhour, &rminute, &rsecond);

	*vdate = cdiEncodeDate(ryear, rmonth, rday);
	*vtime = cdiEncodeTime(rhour, rminute, rsecond);
      }
    }

  return (tstepRange);
}

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

	if ( gridtype == GRID_GAUSSIAN )
          {
            GRIB_CHECK(grib_get_long(gh, "numberOfParallelsBetweenAPoleAndTheEquator", &lpar), 0);
            grid->np = lpar;
          }

	if ( numberOfPoints != nlon*nlat )
	  Error("numberOfPoints (%d) and gridSize (%d) differ!", (int)numberOfPoints, nlon*nlat);

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

        GRIB_CHECK(grib_get_long(gh, "numberOfParallelsBetweenAPoleAndTheEquator", &lpar), 0);
        grid->np = lpar;

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
	    grid->xdef = 2;
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
	    grid->ydef = 2;
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
    case GRID_UNSTRUCTURED:
      {
        char uuid[17];
    	char reference_link[8192];
        size_t len = sizeof(reference_link);
        reference_link[0] = 0;

    	grid->size  = numberOfPoints;
        if ( grib_get_long(gh, "numberOfGridUsed", &lpar) == 0 )
          {
            grid->number   = lpar;
            if ( grib_get_long(gh, "numberOfGridInReference", &lpar) == 0 ) grid->position = lpar;
            /*
            if ( grib_get_string(gh, "gridDescriptionFile", reference_link, &len) == 0 )
              {
                if ( strncmp(reference_link, "file://", 7) == 0 )
                  grid->reference = strdupx(reference_link);
              }
            */
            len = (size_t) 16;
            if ( grib_get_bytes(gh, "uuidOfHGrid", (unsigned char *) uuid, &len) == 0)
              {
                memcpy(grid->uuid, uuid, 16);
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

static
void grib1GetLevel(grib_handle *gh, int *leveltype, int *lbounds, int *level1, int *level2)
{
  int status;
  long lpar;
  double dlevel;

  *leveltype = 0;
  *lbounds = 0;
  *level1  = 0;
  *level2  = 0;

  status = grib_get_long(gh, "indicatorOfTypeOfLevel", &lpar);
  if ( status == 0 )
    {
      *leveltype = (int) lpar;

      switch (*leveltype)
	{
	case GRIB1_LTYPE_SIGMA_LAYER:
	case GRIB1_LTYPE_HYBRID_LAYER:
	case GRIB1_LTYPE_LANDDEPTH_LAYER:
	  { *lbounds = 1; break; }
	}

      if ( *lbounds == 0 )
	{
	  GRIB_CHECK(grib_get_double(gh, "level", &dlevel), 0);
	  if ( *leveltype == 100 ) dlevel *= 100;
	  if ( dlevel < -2.e9 || dlevel > 2.e9 ) dlevel = 0;
	  if ( *leveltype == GRIB1_LTYPE_99 ) *leveltype = 100;

	  *level1 = (int) dlevel;
	  *level2 = 0;
	}
      else
	{
	  GRIB_CHECK(grib_get_long(gh, "topLevel", &lpar), 0);
	  *level1 = lpar;
	  GRIB_CHECK(grib_get_long(gh, "bottomLevel", &lpar), 0);
	  *level2 = lpar;
	}
    }
}

static
double grib2ScaleFactor(long factor)
{
  double scaleFactor = 0;

  if      ( factor == 0 ) scaleFactor =    1;
  else if ( factor == 1 ) scaleFactor =    0.1;
  else if ( factor == 2 ) scaleFactor =    0.01;
  else if ( factor == 3 ) scaleFactor =    0.001;
  else if ( factor == 4 ) scaleFactor =    0.0001;

  return (scaleFactor);
}

static
void grib2GetLevel(grib_handle *gh, int *leveltype, int *lbounds, int *level1, int *level2, int *level_sf, int *level_unit)
{
  int status;
  int leveltype2 = -1;
  long lpar;
  long factor;

  *leveltype  = 0;
  *lbounds    = 0;
  *level1     = 0;
  *level2     = 0;
  *level_sf   = 0;
  *level_unit = 0;

  status = grib_get_long(gh, "typeOfFirstFixedSurface", &lpar);
  if ( status == 0 )
    {
      long llevel;
      double dlevel1 = 0, dlevel2 = 0;

      *leveltype = (int) lpar;

      status = grib_get_long(gh, "typeOfSecondFixedSurface", &lpar);
      if ( status == 0 ) leveltype2 = lpar;

      if ( *leveltype != 255 && leveltype2 != 255 && leveltype2 > 0 ) *lbounds = 1;
      if ( *leveltype == GRIB2_LTYPE_REFERENCE && leveltype2 == 1 ) *lbounds = 0;

      if ( *leveltype == GRIB2_LTYPE_LANDDEPTH )
        {
          *level_sf = 1000;
          *level_unit = CDI_UNIT_M;
        }
      else if ( *leveltype == GRIB2_LTYPE_ISOBARIC )
        {
          *level_sf = 1000;
          *level_unit = CDI_UNIT_PA;
        }

      GRIB_CHECK(grib_get_long(gh, "scaleFactorOfFirstFixedSurface", &factor), 0);
      GRIB_CHECK(grib_get_long(gh, "scaledValueOfFirstFixedSurface", &llevel), 0);
      if ( llevel != GRIB_MISSING_LONG )
        {
          if ( factor != GRIB_MISSING_LONG )
            dlevel1 = llevel*grib2ScaleFactor(factor);
          else
            dlevel1 = llevel;
        }

      if ( *level_sf != 0 ) dlevel1 *= (*level_sf);

      if ( *lbounds == 1 )
	{
          GRIB_CHECK(grib_get_long(gh, "scaleFactorOfSecondFixedSurface", &factor), 0);
          GRIB_CHECK(grib_get_long(gh, "scaledValueOfSecondFixedSurface", &llevel), 0);
          if ( llevel != GRIB_MISSING_LONG )
            {
              if ( factor != GRIB_MISSING_LONG )
                dlevel2 = llevel*grib2ScaleFactor(factor);
              else
                dlevel2 = llevel;
            }

          if ( *level_sf != 0 ) dlevel2 *= (*level_sf);
        }

      *level1 = (int) dlevel1;
      *level2 = (int) dlevel2;
    }
}

static
void gribapiGetString(grib_handle *gh, const char *key, char *string, size_t length)
{
  string[0] = 0;

  GRIB_CHECK(grib_get_string(gh, key, string, &length), 0);
  if      ( length == 8 && memcmp(string, "unknown", length) == 0 ) string[0] = 0;
  else if ( length == 2 && memcmp(string, "~", length)       == 0 ) string[0] = 0;
}

#if  defined  (HAVE_LIBGRIB_API)
static
void gribapiAddRecord(stream_t * streamptr, int param, grib_handle *gh,
		      long recsize, off_t position, int datatype, int comptype, size_t len, const char *varname,
                      int leveltype, int lbounds, int level1, int level2, int level_sf, int level_unit)
{
  long editionNumber;
  int zaxistype;
  int gridID = CDI_UNDEFID, varID;
  int levelID = 0;
  int tsID, recID;
  int numavg;
  int tsteptype;
  record_t *record;
  grid_t grid;
  int vlistID;
  long lpar;
  int status;
  char longname[256], units[256];
  size_t vlen;
  long ens_index = 0, ens_count = 0, ens_forecast_type = 0;

  vlistID = streamptr->vlistID;
  tsID    = streamptr->curTsID;
  recID   = recordNewEntry(streamptr, tsID);
  record  = &streamptr->tsteps[tsID].records[recID];

  tsteptype = gribapiGetTsteptype(gh);
  // numavg  = ISEC1_AvgNum;
  numavg  = 0;

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  // fprintf(stderr, "param %d %d %d %d\n", param, level1, level2, leveltype);

  (*record).size     = recsize;
  (*record).position = position;
  (*record).param    = param;
  (*record).ilevel   = level1;
  (*record).ilevel2  = level2;
  (*record).ltype    = leveltype;
  memcpy((*record).varname, varname, len);

  gribapiGetGrid(gh, &grid);

  gridID = varDefGrid(vlistID, grid, 0);

  zaxistype = gribapiGetZaxisType(editionNumber, leveltype);

  switch (zaxistype)
    {
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
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
        break;
      }
    case ZAXIS_REFERENCE:
      {
        size_t len;
        char uuid[17];
        long ltmp;
        long nhlev, nvgrid;

        GRIB_CHECK(grib_get_long(gh, "NV", &lpar), 0);
        if ( lpar != 6 )
          {
            fprintf(stderr, "Warning ...\n");
          }
        GRIB_CHECK(grib_get_long(gh, "nlev", &ltmp), 0);
        nhlev = ltmp;
        GRIB_CHECK(grib_get_long(gh, "numberOfVGridUsed", &ltmp), 0);
        nvgrid = ltmp;
        len = (size_t) 16;
        uuid[16] = 0;
        GRIB_CHECK(grib_get_bytes(gh, "uuidOfVGrid", (unsigned char *) uuid, &len), 0);
        varDefZAxisReference((int) nhlev, (int) nvgrid, uuid);
        break;
      }
    }

  // if ( datatype > 32 ) datatype = DATATYPE_PACK32;
  if ( datatype <  0 ) datatype = DATATYPE_PACK;

  longname[0] = 0;
  units[0] = 0;

  if ( varname[0] != 0 )
    {
      vlen = 256;
      gribapiGetString(gh, "name", longname, vlen);
      vlen = 256;
      gribapiGetString(gh, "units", units, vlen);
    }
  // fprintf(stderr, "param %d name %s %s %s\n", param, name, longname, units);

  varAddRecord(recID, param, gridID, zaxistype, lbounds, level1, level2, level_sf, level_unit,
	       datatype, &varID, &levelID, tsteptype, numavg, leveltype,
	       varname, longname, units);

  (*record).varID   = varID;
  (*record).levelID = levelID;

  varDefCompType(varID, comptype);

  /*
    Get the ensemble Info from the grib-2 Tables and update the intermediate datastructure.
    Further update to the "vlist" is handled in the same way as for GRIB-1 by "cdi_generate_vars"
  */
  status = grib_get_long(gh, "typeOfEnsembleForecast", &ens_forecast_type );
  if ( status == 0 )
    {
      GRIB_CHECK(grib_get_long(gh, "numberOfForecastsInEnsemble", &ens_count ), 0);
      GRIB_CHECK(grib_get_long(gh, "perturbationNumber", &ens_index ), 0);
    }

  if ( ens_index > 0 )
    varDefEnsembleInfo(varID, (int)ens_index, (int)ens_count, (int)ens_forecast_type);

  long typeOfGeneratingProcess = 0;
  status = grib_get_long(gh, "typeOfGeneratingProcess", &typeOfGeneratingProcess);
  if ( status == 0 )
    varDefTypeOfGeneratingProcess(varID, (int) typeOfGeneratingProcess);

  int    i;
  long   lval;
  double dval;

  /* we read the additional keys for the first variable record only. */
  int linitial_field = (varOptGribNentries(varID) == 0);

  for ( i = 0; i < cdiNAdditionalGRIBKeys; i++ )
    {
      if ( linitial_field )
	{
	  if ( grib_get_long(gh, cdiAdditionalGRIBKeys[i], &lval) == 0 )
            varDefOptGribInt(varID, lval, cdiAdditionalGRIBKeys[i]);
	}
      if ( linitial_field )
	{
	  if ( grib_get_double(gh, cdiAdditionalGRIBKeys[i], &dval) == 0 )
            varDefOptGribDbl(varID, dval, cdiAdditionalGRIBKeys[i]);
	}
      /* note: if the key is not defined, we do not throw an error! */
    }

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

static
int gribapiGetParam(grib_handle *gh)
{
  int pdis = 0, pcat = 0, pnum = 0;
  int param = 0;
  int status;
  long lpar;

  GRIB_CHECK(grib_get_long(gh, "discipline", &lpar), 0);
  pdis = (int) lpar;

  status = grib_get_long(gh, "parameterCategory", &lpar);
  if ( status == 0 ) pcat = (int) lpar;

  status = grib_get_long(gh, "parameterNumber", &lpar);
  if ( status == 0 ) pnum = (int) lpar;

  param = cdiEncodeParam(pnum, pcat, pdis);

  return (param);
}

static
compvar2_t gribapiVarSet(int param, int level1, int level2, int leveltype, char *name)
{
  compvar2_t compVar;
  size_t maxlen = sizeof(compVar.name);
  size_t len = strlen(name);
  if ( len > maxlen ) len = maxlen;

  compVar.param  = param;
  compVar.level1 = level1;
  compVar.level2 = level2;
  compVar.ltype  = leveltype;
  memset(compVar.name, 0, maxlen);
  memcpy(compVar.name, name, len);

  return (compVar);
}

static
int gribapiVarCompare(compvar2_t compVar, record_t record)
{
  int rstatus;
  compvar2_t compVar0;
  size_t maxlen = sizeof(compVar.name);

  compVar0.param  = record.param;
  compVar0.level1 = record.ilevel;
  compVar0.level2 = record.ilevel2;
  compVar0.ltype  = record.ltype;
  memcpy(compVar0.name, record.varname, maxlen);

  rstatus = memcmp(&compVar0, &compVar, sizeof(compvar2_t));

  return (rstatus);
}
#endif

int gribapiScanTimestep1(stream_t * streamptr)
{
#if  defined  (HAVE_LIBGRIB_API)
  off_t recpos = 0;
  unsigned char *gribbuffer = NULL;
  long buffersize = 0;
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
  int nrecs_scanned;
  int datatype;
  long recsize = 0;
  int warn_time = TRUE;
  // int warn_numavg = TRUE;
  int taxisID = -1;
  int rdate = 0, rtime = 0, tunit = 0, fcast = 0;
  taxis_t *taxis;
  int vlistID;
  int comptype;
  long unzipsize;
  compvar2_t compVar;
  grib_handle *gh = NULL;
  int leveltype;
  long editionNumber;
  long lpar;
  size_t len;
  int bitsPerValue;
  int lieee = FALSE;
  int lbounds;
  int level_sf, level_unit;
  char paramstr[32];
  char varname[256];

  streamptr->curTsID = 0;

  tsID  = tstepsNewEntry(streamptr);
  taxis = &streamptr->tsteps[tsID].taxis;

  if ( tsID != 0 )
    Error("Internal problem! tstepsNewEntry returns %d", tsID);

  fileID = streamptr->fileID;

  nrecs_scanned = 0;
  nrecs = 0;
  while ( TRUE )
    {
      level1 = 0;
      level2 = 0;
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

      lieee = FALSE;

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

      nrecs_scanned++;
      gh = grib_handle_new_from_message(NULL, (void *) gribbuffer, recsize);
      GRIB_CHECK(grib_set_double(gh, "missingValue", cdiDefaultMissval), 0);

      GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

      if ( editionNumber <= 1 )
	{
	  GRIB_CHECK(grib_get_long(gh, "table2Version", &lpar), 0);
	  rtabnum = (int) lpar;
	  GRIB_CHECK(grib_get_long(gh, "indicatorOfParameter", &lpar), 0);
	  rcode = (int) lpar;

	  param = cdiEncodeParam(rcode, rtabnum, 255);

	  grib1GetLevel(gh, &leveltype, &lbounds, &level1, &level2);
          level_sf = 0;
          level_unit = 0;
	}
      else
	{
	  size_t len = 256;
	  char typeOfPacking[256];

	  status = grib_get_string(gh, "packingType", typeOfPacking, &len);
	  if ( status == 0 )
	    {
	      // fprintf(stderr, "packingType %d %s\n", len, typeOfPacking);
	      if      ( strncmp(typeOfPacking, "grid_jpeg", len) == 0 ) comptype = COMPRESS_JPEG;
	      else if ( strncmp(typeOfPacking, "grid_ieee", len) == 0 ) lieee = TRUE;
	    }

	  param = gribapiGetParam(gh);

	  grib2GetLevel(gh, &leveltype, &lbounds, &level1, &level2, &level_sf, &level_unit);
	}

      cdiParamToString(param, paramstr, sizeof(paramstr));

      varname[0] = 0;
      gribapiGetString(gh, "shortName", varname, sizeof(varname));
      len = strlen(varname);
      if ( len > 32 ) len = 32;
      //printf("param = %s  name = %s   l1 = %d  l2 = %d\n", paramstr, varname, level1, level2);

      gribapiGetValidityDateTime(gh, &vdate, &vtime);
      /*
      printf("%d %d %d\n", vdate, vtime, leveltype);
      */
      if ( lieee )
        {
          datatype = DATATYPE_FLT64;
          status = grib_get_long(gh, "precision", &lpar);
          if ( status == 0 && lpar == 1 ) datatype = DATATYPE_FLT32;
        }
      else
        {
          datatype = DATATYPE_PACK;
          status = grib_get_long(gh, "bitsPerValue", &lpar);
          if ( status == 0 )
            {
              bitsPerValue = (int) lpar;
              if ( bitsPerValue > 0 && bitsPerValue <= 32 )
                datatype = bitsPerValue;
            }
        }

      if ( nrecs == 0 )
	{
	  datetime0.date = vdate;
	  datetime0.time = vtime;

          gribapiGetDataDateTime(gh, &rdate, &rtime);

	  fcast = gribapiTimeIsFC(gh);
	  if ( fcast ) tunit = gribapiGetTimeUnits(gh);
	}
      else
	{
	  datetime.date  = vdate;
	  datetime.time  = vtime;

	  compVar = gribapiVarSet(param, level1, level2, leveltype, varname);

	  for ( recID = 0; recID < nrecs; recID++ )
            if ( gribapiVarCompare(compVar, streamptr->tsteps[0].records[recID]) == 0 ) break;

	  if ( cdiInventoryMode == 1 )
	    {
	      if ( recID < nrecs ) break;
	      if ( warn_time )
		if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 )
		  {
                    if ( datetime0.date == 10101 && datetime0.time == 0 )
                      {
                        datetime0.date = datetime.date;
                        datetime0.time = datetime.time;

                        gribapiGetDataDateTime(gh, &rdate, &rtime);

                        fcast = gribapiTimeIsFC(gh);
                        if ( fcast ) tunit = gribapiGetTimeUnits(gh);
                      }
                    else
                      {
                        Warning("Inconsistent verification time (param=%s level=%d)", paramstr, level1);
                        warn_time = FALSE;
                      }
                  }
	    }
	  else
	    {
	      if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 ) break;

	      if ( recID < nrecs )
		{
		  Warning("Param=%s level=%d (record %d) already exist, skipped!", paramstr, level1, nrecs_scanned);
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

      gribapiAddRecord(streamptr, param, gh, recsize, recpos, datatype, comptype, len, varname,
                       leveltype, lbounds, level1, level2, level_sf, level_unit);

      grib_handle_delete(gh);
      gh = NULL;
    }

  if ( gh ) grib_handle_delete(gh);

  streamptr->rtsteps = 1;

  if ( nrecs == 0 ) return (CDI_EUFSTRUCT);

  cdi_generate_vars(streamptr);

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

  vlistID = streamptr->vlistID;
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
      tsID = tstepsNewEntry(streamptr);
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
	      vlistDefVarTsteptype(vlistID, varID, TSTEP_CONSTANT);
	    }
	}
    }
#else
  Error("GRIB_API support not compiled in!");
#endif

  return (0);
}


int gribapiScanTimestep2(stream_t * streamptr)
{
  int rstatus = 0;
#if  defined  (HAVE_LIBGRIB_API)
  off_t recpos = 0;
  unsigned char *gribbuffer = NULL;
  long buffersize = 0;
  int fileID;
  int rtabnum = 0;
  int rcode = 0, level1 = 0, level2 = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0;
  int tsID;
  int varID;
  // int gridID;
  size_t readsize;
  int nrecords, nrecs, recID, rindex;
  long recsize = 0;
  //  int warn_numavg = TRUE;
  int tsteptype;
  int taxisID = -1;
  taxis_t *taxis;
  int vlistID;
  long unzipsize;
  compvar2_t compVar;
  grib_handle *gh = NULL;
  int leveltype;
  int param = 0;
  long editionNumber;
  long lpar;
  int lbounds;
  int level_sf, level_unit;
  char paramstr[32];
  char varname[256];

  streamptr->curTsID = 1;

  fileID  = streamptr->fileID;
  vlistID = streamptr->vlistID;
  taxisID = vlistInqTaxis(vlistID);

  gribbuffer = (unsigned char *) streamptr->record->buffer;
  buffersize = streamptr->record->buffersize;

  tsID = streamptr->rtsteps;
  if ( tsID != 1 )
    Error("Internal problem! unexpeceted timestep %d", tsID+1);

  taxis = &streamptr->tsteps[tsID].taxis;

  fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

  cdi_create_records(streamptr, tsID);

  nrecords = streamptr->tsteps[tsID].nallrecs;
  streamptr->tsteps[1].recIDs = (int *) malloc(nrecords*sizeof(int));
  streamptr->tsteps[1].nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[1].recIDs[recID] = -1;

  for ( recID = 0; recID < nrecords; recID++ )
    {
      varID = streamptr->tsteps[0].records[recID].varID;
      streamptr->tsteps[tsID].records[recID].position = streamptr->tsteps[0].records[recID].position;
      streamptr->tsteps[tsID].records[recID].size     = streamptr->tsteps[0].records[recID].size;
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
      GRIB_CHECK(grib_set_double(gh, "missingValue", cdiDefaultMissval), 0);

      GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

      if ( editionNumber <= 1 )
	{
	  GRIB_CHECK(grib_get_long(gh, "table2Version", &lpar), 0);
	  rtabnum = (int) lpar;
	  GRIB_CHECK(grib_get_long(gh, "indicatorOfParameter", &lpar), 0);
	  rcode = (int) lpar;

	  param = cdiEncodeParam(rcode, rtabnum, 255);

	  grib1GetLevel(gh, &leveltype, &lbounds, &level1, &level2);
          level_sf = 0;
          level_unit = 0;
	}
      else
	{
	  param = gribapiGetParam(gh);

	  grib2GetLevel(gh, &leveltype, &lbounds, &level1, &level2, &level_sf, &level_unit);
	}

      cdiParamToString(param, paramstr, sizeof(paramstr));

      varname[0] = 0;
      gribapiGetString(gh, "shortName", varname, sizeof(varname));

      gribapiGetValidityDateTime(gh, &vdate, &vtime);

      if ( rindex == 0 )
	{
	  if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
	    {
	      taxis->type  = TAXIS_RELATIVE;

              gribapiGetDataDateTime(gh, &(taxis->rdate), &(taxis->rtime));

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

      compVar = gribapiVarSet(param, level1, level2, leveltype, varname);

      for ( recID = 0; recID < nrecords; recID++ )
        if ( gribapiVarCompare(compVar, streamptr->tsteps[tsID].records[recID]) == 0 ) break;

      if ( recID == nrecords )
	{
	  Warning("Param=%s (%s) l1=%d l2=%d not defined at timestep 1!", paramstr, varname, level1, level2);
	  return (CDI_EUFSTRUCT);
	}

      if ( streamptr->tsteps[tsID].records[recID].used )
        {
          if ( cdiInventoryMode == 1 ) break;
          else
	    {
	      if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 ) break;

	      Warning("Param=%s level=%d already exist, skipped!", paramstr, level1);
	      continue;
	    }
	}

      streamptr->tsteps[tsID].records[recID].used = TRUE;
      streamptr->tsteps[tsID].recIDs[rindex] = recID;

      if ( CDI_Debug )
	Message("%4d %8d %4d %8d %8d %6d", rindex+1, (int)recpos, param, level1, vdate, vtime);

      streamptr->tsteps[tsID].records[recID].size = recsize;

      if ( gribapiVarCompare(compVar, streamptr->tsteps[tsID].records[recID]) != 0 )
	{
	  Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		  tsID, recID,
		  streamptr->tsteps[tsID].records[recID].param, param,
		  streamptr->tsteps[tsID].records[recID].ilevel, level1);
	  return (CDI_EUFSTRUCT);
	}

      streamptr->tsteps[1].records[recID].position = recpos;
      varID = streamptr->tsteps[tsID].records[recID].varID;
      /*
      gridID = vlistInqVarGrid(vlistID, varID);
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
	  vlistDefVarTsteptype(vlistID, varID, TSTEP_CONSTANT);
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
      tsID = tstepsNewEntry(streamptr);
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


int gribapiScanTimestep(stream_t * streamptr)
{
  int rstatus = 0;
#if  defined  (HAVE_LIBGRIB_API)
  long recsize = 0;
  off_t recpos = 0;
  unsigned char *gribbuffer;
  long buffersize = 0;
  int fileID;
  int rtabnum = 0;
  int rcode = 0, level1 = 0, level2 = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0;
  int tsID;
  int vrecID, recID;
  //int warn_numavg = TRUE;
  size_t readsize;
  int taxisID = -1;
  taxis_t *taxis;
  int vlistID;
  int rindex, nrecs = 0;
  long unzipsize;
  compvar2_t compVar;
  grib_handle *gh = NULL;
  int leveltype;
  int param = 0;
  long editionNumber;
  long lpar;
  int lbounds;
  int level_sf, level_unit;
  char paramstr[32];
  char varname[256];

  vlistID = streamptr->vlistID;

  if ( CDI_Debug )
    {
      Message("streamID = %d", streamptr->self);
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

      cdi_create_records(streamptr, tsID);

      nrecs = streamptr->tsteps[1].nrecs;

      streamptr->tsteps[tsID].nrecs = nrecs;
      streamptr->tsteps[tsID].recIDs = (int *) malloc(nrecs*sizeof(int));
      for ( recID = 0; recID < nrecs; recID++ )
	streamptr->tsteps[tsID].recIDs[recID] = streamptr->tsteps[1].recIDs[recID];

      fileID = streamptr->fileID;

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
	      Warning("Inconsistent timestep %d (GRIB record %d/%d)!", tsID+1, rindex+1,
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
	  GRIB_CHECK(grib_set_double(gh, "missingValue", cdiDefaultMissval), 0);

	  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

	  if ( editionNumber <= 1 )
	    {
	      GRIB_CHECK(grib_get_long(gh, "table2Version", &lpar), 0);
	      rtabnum = (int) lpar;
	      GRIB_CHECK(grib_get_long(gh, "indicatorOfParameter", &lpar), 0);
	      rcode = (int) lpar;

	      param = cdiEncodeParam(rcode, rtabnum, 255);

	      grib1GetLevel(gh, &leveltype, &lbounds, &level1, &level2);
              level_sf = 0;
              level_unit = 0;
	    }
	  else
	    {
	      param = gribapiGetParam(gh);

	      grib2GetLevel(gh, &leveltype, &lbounds, &level1, &level2, &level_sf, &level_unit);
	    }

          cdiParamToString(param, paramstr, sizeof(paramstr));

          varname[0] = 0;
	  gribapiGetString(gh, "shortName", varname, sizeof(varname));

	  gribapiGetValidityDateTime(gh, &vdate, &vtime);

	  if ( rindex == nrecs ) break;

	  if ( rindex == 0 )
	    {
	      taxisID = vlistInqTaxis(vlistID);
	      if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
		{
		  taxis->type  = TAXIS_RELATIVE;

                  gribapiGetDataDateTime(gh, &(taxis->rdate), &(taxis->rtime));

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

	  compVar = gribapiVarSet(param, level1, level2, leveltype, varname);

	  for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	    {
	      recID   = streamptr->tsteps[1].recIDs[vrecID];
	      if ( gribapiVarCompare(compVar, streamptr->tsteps[tsID].records[recID]) == 0 ) break;
	    }

	  if ( vrecID == nrecs )
	    {
	      Warning("Param=%s level=%d not available at timestep %d!", paramstr, level1, tsID+1);

	      if ( cdiInventoryMode == 1 )
		return (CDI_EUFSTRUCT);
	      else
		continue;
	    }

	  if ( cdiInventoryMode != 1 )
	    {
	      if ( streamptr->tsteps[tsID].records[recID].used )
		{
		  if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 ) break;

		  if ( CDI_Debug )
		    Warning("Param=%s level=%d already exist, skipped!", paramstr, level1);

		  continue;
		}
	    }

          streamptr->tsteps[tsID].records[recID].used = TRUE;
          streamptr->tsteps[tsID].recIDs[rindex] = recID;

	  if ( CDI_Debug )
	    Message("%4d %8d %4d %8d %8d %6d", rindex+1, (int)recpos, param, level1, vdate, vtime);

	  if ( gribapiVarCompare(compVar, streamptr->tsteps[tsID].records[recID]) != 0 )
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
	  cdiParamToString(streamptr->tsteps[tsID].records[recID].param, paramstr, sizeof(paramstr));
	  Warning("Param %d level %d not found at timestep %d!",
		  paramstr, streamptr->tsteps[tsID].records[recID].ilevel, tsID+1);
	  return (CDI_EUFSTRUCT);
	}

      streamptr->rtsteps++;

      if ( streamptr->ntsteps != streamptr->rtsteps )
	{
	  tsID = tstepsNewEntry(streamptr);
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
		  int unreduced, int *nmiss, int *zip, double missval, int vlistID, int varID)
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

static
void gribapiDefParam(int editionNumber, grib_handle *gh, int param, const char *name)
{
  int pdis, pcat, pnum;

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

      if ( editionNumber <= 1 )
	{
	  if ( pdis != 255 )
	    {
	      char paramstr[32];
	      cdiParamToString(param, paramstr, sizeof(paramstr));
	      Warning("Can't convert GRIB2 parameter ID (%s) to GRIB1, set to %d.%d!", paramstr, pnum, pcat);
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

static
int gribapiDefStepUnits(grib_handle *gh, int timeunit, int gcinit)
{
  int factor = 1;
  long unitsOfTime;
  char stepunits[8];
  size_t len;

  switch (timeunit)
    {
    case TUNIT_SECOND:  factor =     1;  unitsOfTime = 13;  strcpy(stepunits, "s");   break;
    case TUNIT_MINUTE:  factor =    60;  unitsOfTime =  0;  strcpy(stepunits, "m");   break;
    case TUNIT_HOUR:    factor =  3600;  unitsOfTime =  1;  strcpy(stepunits, "h");   break;
    case TUNIT_3HOURS:  factor = 10800;  unitsOfTime = 10;  strcpy(stepunits, "3h");  break;
    case TUNIT_6HOURS:  factor = 21600;  unitsOfTime = 11;  strcpy(stepunits, "6h");  break;
    case TUNIT_12HOURS: factor = 43200;  unitsOfTime = 12;  strcpy(stepunits, "12h"); break;
    case TUNIT_DAY:     factor = 86400;  unitsOfTime =  2;  strcpy(stepunits, "D");   break;
    default:            factor =  3600;  unitsOfTime =  1;  strcpy(stepunits, "h");   break;
    }

  if ( !gcinit )
    {
      len = strlen(stepunits) + 1;
      GRIB_CHECK(grib_set_long(gh, "indicatorOfUnitOfTimeRange", unitsOfTime), 0);
      GRIB_CHECK(grib_set_string(gh, "stepUnits", stepunits, &len), 0);
    }

  return (factor);
}

static
int gribapiDefSteptype(int editionNumber, grib_handle *gh, int tsteptype, int gcinit)
{
  long proDefTempNum = 0;
  size_t len = 64;
  char stepType[64];

  switch ( tsteptype )
    {
    case TSTEP_AVG:      strcpy(stepType, "avg");     proDefTempNum = 8; break;
    case TSTEP_ACCUM:    strcpy(stepType, "accum");   proDefTempNum = 8; break;
    case TSTEP_MAX:      strcpy(stepType, "max");     proDefTempNum = 8; break;
    case TSTEP_MIN:      strcpy(stepType, "min");     proDefTempNum = 8; break;
    case TSTEP_DIFF:     strcpy(stepType, "diff");    proDefTempNum = 8; break;
    case TSTEP_RMS:      strcpy(stepType, "rms");     proDefTempNum = 8; break;
    case TSTEP_SD:       strcpy(stepType, "sd");      proDefTempNum = 8; break;
    case TSTEP_COV:      strcpy(stepType, "cov");     proDefTempNum = 8; break;
    case TSTEP_RATIO:    strcpy(stepType, "ratio");   proDefTempNum = 8; break;
    case TSTEP_INSTANT:  strcpy(stepType, "instant"); proDefTempNum = 0; break;
    default:             strcpy(stepType, "instant"); proDefTempNum = 0; break;
    }

  if ( !gcinit )
    {
      if ( editionNumber > 1 ) GRIB_CHECK(grib_set_long(gh, "productDefinitionTemplateNumber", proDefTempNum), 0);
      len = strlen(stepType);
      GRIB_CHECK(grib_set_string(gh, "stepType", stepType, &len), 0);
    }

  return ((int)proDefTempNum);
}

static
void gribapiDefDateTimeAbs(int editionNumber, grib_handle *gh, int date, int time, int tsteptype, int gcinit)
{
  if ( editionNumber > 1 ) GRIB_CHECK(grib_set_long(gh, "significanceOfReferenceTime", 0), 0);
  if ( editionNumber > 1 ) GRIB_CHECK(grib_set_long(gh, "stepRange", 0), 0);

  if ( date == 0 ) date = 10101;
  gribapiSetDataDateTime(gh, date, time);

  (void ) gribapiDefSteptype(editionNumber, gh, tsteptype, gcinit);
}

static
int gribapiDefDateTimeRel(int editionNumber, grib_handle *gh, int rdate, int rtime, int vdate, int vtime,
                          int tsteptype, int factor, int calendar, int gcinit)
{
  int status = -1;
  int year, month, day, hour, minute, second;
  int julday1, secofday1, julday2, secofday2, days, secs;
  long startStep = 0, endStep;
  long proDefTempNum = 0;

  cdiDecodeDate(rdate, &year, &month, &day);
  cdiDecodeTime(rtime, &hour, &minute, &second);
  encode_juldaysec(calendar, year, month, day, hour, minute, second, &julday1, &secofday1);

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);
  encode_juldaysec(calendar, year, month, day, hour, minute, second, &julday2, &secofday2);

  (void) julday_sub(julday1, secofday1, julday2, secofday2, &days, &secs);

  if ( !(int) fmod(days*86400.0 + secs, factor) )
    {
      endStep = (int) ((days*86400.0 + secs)/factor);

      if ( editionNumber > 1 ) GRIB_CHECK(grib_set_long(gh, "significanceOfReferenceTime", 1), 0);
      if ( editionNumber > 1 ) GRIB_CHECK(grib_set_long(gh, "stepRange", 0), 0);

      if ( rdate == 0 ) rdate = 10101;
      gribapiSetDataDateTime(gh, rdate, rtime);

      // printf(">>>>> tsteptype %d  startStep %ld  endStep %ld\n", tsteptype, startStep, endStep);

      proDefTempNum = gribapiDefSteptype(editionNumber, gh, tsteptype, gcinit);

      if ( proDefTempNum == 0 ) startStep = endStep;

      if ( editionNumber > 1 ) GRIB_CHECK(grib_set_long(gh, "forecastTime", startStep), 0);
      GRIB_CHECK(grib_set_long(gh, "endStep", endStep), 0);

      status = 0;
    }

  return (status);
}

static
void gribapiDefTime(int editionNumber, int typeOfGeneratingProcess, grib_handle *gh , int vdate, int vtime, int tsteptype, int numavg, int taxisID, int gcinit)
{
  int taxistype = -1;

  if ( taxisID != -1 ) taxistype = taxisInqType(taxisID);

  if ( typeOfGeneratingProcess == 196 )
    {
      vdate = 10101;
      vtime = 0;
      taxistype = TAXIS_ABSOLUTE;
    }
  /*
  else if ( typeOfGeneratingProcess == 9 )
    {
    }
  */

  if ( taxistype == TAXIS_RELATIVE )
    {
      int status;
      int calendar = taxisInqCalendar(taxisID);
      int rdate    = taxisInqRdate(taxisID);
      int rtime    = taxisInqRtime(taxisID);
      int timeunit = taxisInqTunit(taxisID);
      int factor   = gribapiDefStepUnits(gh, timeunit, gcinit);

      status = gribapiDefDateTimeRel(editionNumber, gh, rdate, rtime, vdate, vtime,
                                     tsteptype, factor, calendar, gcinit);

      if ( status != 0 ) taxistype = TAXIS_ABSOLUTE;
    }

  if ( taxistype == TAXIS_ABSOLUTE )
    {
      gribapiDefDateTimeAbs(editionNumber, gh, vdate, vtime, tsteptype, gcinit);
    }
}

static
void gribapiDefGrid(int editionNumber, grib_handle *gh, int gridID, int ljpeg, int lieee, int datatype, int nmiss, int gcinit)
{
  int gridtype;
  int status;
  static short lwarn = TRUE;
  size_t len;
  char *mesg;

  gridtype = gridInqType(gridID);

  if ( editionNumber <= 1 )
    if ( gridtype == GRID_GME || gridtype == GRID_UNSTRUCTURED )
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


  if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
    {
      if ( editionNumber != 2 || lieee ) { ljpeg = 0; }

      if ( ljpeg )
        {
          if ( nmiss > 0 ) ljpeg = 0;

          if ( ljpeg )
            {
              mesg = "grid_jpeg"; len = strlen(mesg);
              GRIB_CHECK(grib_set_string(gh, "packingType", mesg, &len), 0);
            }
          else
            {
              mesg = "grid_simple"; len = strlen(mesg);
              GRIB_CHECK(grib_set_string(gh, "packingType", mesg, &len), 0);
            }
        }
    }

  if ( gcinit ) return;

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
          {
            int np = gridInqNP(gridID);
            if ( np == 0 ) np = nlat/2;
            GRIB_CHECK(grib_set_long(gh, "numberOfParallelsBetweenAPoleAndTheEquator", np), 0);
          }
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

        if ( editionNumber != 2 ) { lieee = 0; ljpeg = 0; }

        if ( lieee )
          {
            mesg = "grid_ieee"; len = strlen(mesg);
            GRIB_CHECK(grib_set_string(gh, "packingType", mesg, &len), 0);

	    if ( datatype == DATATYPE_FLT64 )
	      GRIB_CHECK(grib_set_long(gh, "precision", 2), 0);
	    else
	      GRIB_CHECK(grib_set_long(gh, "precision", 1), 0);
          }
        else if ( ljpeg )
	  {
            if ( nmiss > 0 ) ljpeg = 0;

            if ( ljpeg )
              {
                mesg = "grid_jpeg"; len = strlen(mesg);
                GRIB_CHECK(grib_set_string(gh, "packingType", mesg, &len), 0);
              }
            else
              {
                mesg = "grid_simple"; len = strlen(mesg);
                GRIB_CHECK(grib_set_string(gh, "packingType", mesg, &len), 0);
              }
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
        /*
        if ( lieee )
          {
            printf("spectral_ieee\n");
            if ( editionNumber == 2 ) GRIB_CHECK(grib_set_long(gh, "numberOfValues", gridInqSize(gridID)), 0);
            mesg = "spectral_ieee"; len = strlen(mesg);
            GRIB_CHECK(grib_set_string(gh, "packingType", mesg, &len), 0);
          }
        else */ if ( gridInqComplexPacking(gridID) )
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
	GRIB_CHECK(grib_set_long(gh, "latitudeOfThePolePoint", 90000000), 0);
	GRIB_CHECK(grib_set_long(gh, "longitudeOfThePolePoint", 0), 0);

	GRIB_CHECK(grib_set_long(gh, "numberOfDataPoints", gridInqSize(gridID)), 0);
	GRIB_CHECK(grib_set_long(gh, "totalNumberOfGridPoints", gridInqSize(gridID)), 0);

	break;
      }
    case GRID_UNSTRUCTURED:
      {
	static int warning = 1;
	status = grib_set_long(gh, "gridDefinitionTemplateNumber", GRIB2_GTYPE_UNSTRUCTURED);
	if ( status != 0 && warning )
	  {
	    warning = 0;
	    Warning("Can't write reference grid!");
	    Warning("gridDefinitionTemplateNumber %d not found (grib2/template.3.%d.def)!",
		    GRIB2_GTYPE_UNSTRUCTURED, GRIB2_GTYPE_UNSTRUCTURED);
	  }
	else
	  {
            char uuid[17];
            int position = gridInqPosition(gridID);
            int number = gridInqNumber(gridID);
            if ( position < 0 ) position = 0;
            if ( number < 0 ) number = 0;
	    GRIB_CHECK(grib_set_long(gh, "numberOfGridUsed", number), 0);
	    GRIB_CHECK(grib_set_long(gh, "numberOfGridInReference", position), 0);
            len = 16;
            gridInqUUID(gridID, uuid);
	    if (grib_set_bytes(gh, "uuidOfHGrid", (unsigned char *) uuid, &len) != 0)
	      Warning("Can't write UUID!");
	  }

	break;
      }
    default:
      {
	Error("Unsupported grid type: %s", gridNamePtr(gridtype));
	break;
      }
    }
}

static
void getLevelFactor(double level, long *factor, long *out_scaled_value)
{
  double scaled_value  = level;
  long   iscaled_value = (long) round(scaled_value);
  long   i;

  const double eps = 1.e-8;
  for ( i=0; (fabs(scaled_value - (double) iscaled_value) >= eps) && i < 7; i++ )
    {
      scaled_value *= 10.;
      iscaled_value = round(scaled_value);
    }

  (*factor)           = i;
  (*out_scaled_value) = iscaled_value;
}

static
void gribapiDefLevelType(grib_handle *gh, int gcinit, const char *keyname, long leveltype)
{
  if ( !gcinit ) GRIB_CHECK(grib_set_long(gh, keyname, leveltype), 0);
}

static
void grib2DefLevel(grib_handle *gh, int gcinit, long leveltype, int lbounds, double level, double dlevel1, double dlevel2)
{
  long scaled_level;
  long factor;

  gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", leveltype);
  if ( lbounds ) gribapiDefLevelType(gh, gcinit, "typeOfSecondFixedSurface", leveltype);

  if ( !lbounds ) dlevel1 = level;

  getLevelFactor(dlevel1, &factor, &scaled_level);
  GRIB_CHECK(grib_set_long(gh, "scaleFactorOfFirstFixedSurface", factor), 0);
  GRIB_CHECK(grib_set_long(gh, "scaledValueOfFirstFixedSurface", scaled_level), 0);

  if ( lbounds )
    {
      getLevelFactor(dlevel2, &factor, &scaled_level);
      GRIB_CHECK(grib_set_long(gh, "scaleFactorOfSecondFixedSurface", factor), 0);
      GRIB_CHECK(grib_set_long(gh, "scaledValueOfSecondFixedSurface", scaled_level), 0);
    }
}

static
void gribapiDefLevel(int editionNumber, grib_handle *gh, int param, int zaxisID, int levelID, int gcinit)
{
  double level;
  int lbounds = 0;
  int zaxistype, ltype;
  static int warning = 1;
  char uuid[17];
  size_t len;
  double scalefactor;
  double dlevel1 = 0, dlevel2 = 0;


  zaxistype = zaxisInqType(zaxisID);
  ltype = zaxisInqLtype(zaxisID);
  level = zaxisInqLevel(zaxisID, levelID);

  if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
    {
      lbounds = 1;
      dlevel1 = zaxisInqLbound(zaxisID, levelID);
      dlevel2 = zaxisInqUbound(zaxisID, levelID);
    }
  else
    {
      dlevel1 = level;
      dlevel2 = 0;
    }

  if ( zaxistype == ZAXIS_GENERIC && ltype == 0 )
    {
      Message("Changed zaxis type from %s to %s", zaxisNamePtr(zaxistype), zaxisNamePtr(ZAXIS_PRESSURE));
      zaxistype = ZAXIS_PRESSURE;
      zaxisChangeType(zaxisID, zaxistype);
      zaxisDefUnits(zaxisID, "Pa");
    }

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:
    case ZAXIS_MEANSEA:
    case ZAXIS_HEIGHT:
    case ZAXIS_ALTITUDE:
    case ZAXIS_SIGMA:
    case ZAXIS_DEPTH_BELOW_SEA:
    case ZAXIS_ISENTROPIC:
      {
	if ( editionNumber <= 1 )
          gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", zaxisTypeToGrib1ltype(zaxistype));
        else
          gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", zaxisTypeToGrib2ltype(zaxistype));

        GRIB_CHECK(grib_set_long(gh, "level", level), 0);

	break;
      }
    case ZAXIS_CLOUD_BASE:
    case ZAXIS_CLOUD_TOP:
    case ZAXIS_ISOTHERM_ZERO:
    case ZAXIS_TOA:
    case ZAXIS_SEA_BOTTOM:
    case ZAXIS_LAKE_BOTTOM:
    case ZAXIS_SEDIMENT_BOTTOM:
    case ZAXIS_SEDIMENT_BOTTOM_TA:
    case ZAXIS_SEDIMENT_BOTTOM_TW:
    case ZAXIS_MIX_LAYER:
    case ZAXIS_ATMOSPHERE:
      {
        if ( lbounds )
          {
            if ( editionNumber <= 1 )
              gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", zaxisTypeToGrib1ltype(zaxistype));
            else
              {
                gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", zaxisTypeToGrib2ltype(zaxistype));
                gribapiDefLevelType(gh, gcinit, "typeOfSecondFixedSurface", zaxisTypeToGrib2ltype(zaxistype));
              }

            GRIB_CHECK(grib_set_long(gh, "topLevel", (long) dlevel1), 0);
            GRIB_CHECK(grib_set_long(gh, "bottomLevel", (long) dlevel2), 0);
          }
        else
          {
            if ( editionNumber <= 1 )
              gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", zaxisTypeToGrib1ltype(zaxistype));
            else
              gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", zaxisTypeToGrib2ltype(zaxistype));

            GRIB_CHECK(grib_set_long(gh, "level", (long) level), 0);
          }

        break;
      }
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
      {
	if ( lbounds )
	  {
	    if ( editionNumber <= 1 )
	      gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", GRIB1_LTYPE_HYBRID_LAYER);
            else
	      {
		gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", GRIB2_LTYPE_HYBRID);
		gribapiDefLevelType(gh, gcinit, "typeOfSecondFixedSurface", GRIB2_LTYPE_HYBRID);
	      }

	    GRIB_CHECK(grib_set_long(gh, "topLevel", (long) dlevel1), 0);
	    GRIB_CHECK(grib_set_long(gh, "bottomLevel", (long) dlevel2), 0);
	  }
	else
	  {
	    if ( editionNumber <= 1 )
              gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", GRIB1_LTYPE_HYBRID);
            else
              gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", GRIB2_LTYPE_HYBRID);

	    GRIB_CHECK(grib_set_long(gh, "level", (long) level), 0);
	  }

        if ( !gcinit )
          {
            int vctsize = zaxisInqVctSize(zaxisID);
            if ( vctsize == 0 && warning )
              {
                char paramstr[32];
                cdiParamToString(param, paramstr, sizeof(paramstr));
                Warning("VCT missing ( param = %s, zaxisID = %d )", paramstr, zaxisID);
                warning = 0;
              }
            GRIB_CHECK(grib_set_long(gh, "PVPresent", 1), 0);
            GRIB_CHECK(grib_set_double_array(gh, "pv", zaxisInqVctPtr(zaxisID), vctsize), 0);
          }

	break;
      }
    case ZAXIS_PRESSURE:
      {
	double dum;
	char units[128];

	if ( level < 0 ) Warning("Pressure level of %f Pa is below zero!", level);

	zaxisInqUnits(zaxisID, units);
	if ( memcmp(units, "Pa", 2) != 0 )
          {
            level   *= 100;
            dlevel1 *= 100;
            dlevel2 *= 100;
          }

        if ( editionNumber <= 1 )
          {
            long leveltype = GRIB1_LTYPE_ISOBARIC;

            if ( level < 32768 && (level < 100 || modf(level/100, &dum) > 0) )
              leveltype = GRIB1_LTYPE_99;
            else
              level /= 100;

            gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", leveltype);
            GRIB_CHECK(grib_set_double(gh, "level", level), 0);
	  }
	else
	  {
            grib2DefLevel(gh, gcinit, GRIB2_LTYPE_ISOBARIC, lbounds, level, dlevel1, dlevel2);
	  }

	break;
      }
    case ZAXIS_SNOW:
      {
        if ( editionNumber <= 1 )
          ; // not available
	else
          {
            grib2DefLevel(gh, gcinit, GRIB2_LTYPE_SNOW, lbounds, level, dlevel1, dlevel2);
          }

	break;
      }
    case ZAXIS_DEPTH_BELOW_LAND:
      {
	char units[128];

	zaxisInqUnits(zaxisID, units);

	if ( editionNumber <= 1 )
	  {
	    if      ( memcmp(units, "mm", 2) == 0 ) scalefactor =   0.1;
	    else if ( memcmp(units, "cm", 2) == 0 ) scalefactor =   1; // cm
	    else if ( memcmp(units, "dm", 2) == 0 ) scalefactor =  10;
	    else                                    scalefactor = 100;

	    gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", GRIB1_LTYPE_LANDDEPTH);
	    GRIB_CHECK(grib_set_double(gh, "level", level*scalefactor), 0);
	  }
	else
	  {
	    if      ( memcmp(units, "mm", 2) == 0 ) scalefactor = 0.001;
	    else if ( memcmp(units, "cm", 2) == 0 ) scalefactor = 0.01;
	    else if ( memcmp(units, "dm", 2) == 0 ) scalefactor = 0.1;
	    else                                    scalefactor = 1; // meter

            level   *= scalefactor;
            dlevel1 *= scalefactor;
            dlevel2 *= scalefactor;

            grib2DefLevel(gh, gcinit, GRIB2_LTYPE_LANDDEPTH, lbounds, level, dlevel1, dlevel2);
	  }

	break;
      }
    case ZAXIS_REFERENCE:
      {
        int number;

        if ( !gcinit )
          {
            GRIB_CHECK(grib_set_long(gh, "genVertHeightCoords", 1), 0);
          }

        if ( lbounds )
          {
            if ( editionNumber <= 1 )
              ; // not available
            else
              {
                number = zaxisInqNumber(zaxisID);
                gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", GRIB2_LTYPE_REFERENCE);
                gribapiDefLevelType(gh, gcinit, "typeOfSecondFixedSurface", GRIB2_LTYPE_REFERENCE);
                GRIB_CHECK(grib_set_long(gh, "NV", 6), 0);
                GRIB_CHECK(grib_set_long(gh, "nlev", zaxisInqNlevRef(zaxisID)), 0);
                GRIB_CHECK(grib_set_long(gh, "numberOfVGridUsed", number), 0);
                len = 16;
                zaxisInqUUID(zaxisID, uuid);
                if (grib_set_bytes(gh, "uuidOfVGrid", (unsigned char *) uuid, &len) != 0)
                  Warning("Can't write UUID!");
                GRIB_CHECK(grib_set_long(gh, "topLevel", (long) dlevel1), 0);
                GRIB_CHECK(grib_set_long(gh, "bottomLevel", (long) dlevel2), 0);
              }
          }
        else
          {
            if ( editionNumber <= 1 )
              ; // not available
            else
              {
                number = zaxisInqNumber(zaxisID);
                gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", GRIB2_LTYPE_REFERENCE);
                GRIB_CHECK(grib_set_long(gh, "NV", 6), 0);
                GRIB_CHECK(grib_set_long(gh, "nlev", zaxisInqNlevRef(zaxisID)), 0);
                GRIB_CHECK(grib_set_long(gh, "numberOfVGridUsed", number), 0);
                len = 16;
                zaxisInqUUID(zaxisID, uuid);
                if (grib_set_bytes(gh, "uuidOfVGrid", (unsigned char *) uuid, &len) != 0)
                  Warning("Can't write UUID!");
                GRIB_CHECK(grib_set_double(gh, "level", level), 0);
              }
          }

        break;
      }
    case ZAXIS_GENERIC:
      {
	if ( editionNumber <= 1 )
          gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", ltype);
        else
          gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", ltype);

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
		     long datasize, const double *data, int nmiss, unsigned char **gribbuffer, size_t *gribbuffersize,
		     int ljpeg, void *gribContainer)
{
  size_t nbytes = 0;
#if  defined  (HAVE_LIBGRIB_API)
  size_t recsize = 0;
  void *dummy = NULL;
  int datatype;
  int param;
  int lieee = FALSE;
  int ensID, ensCount, forecast_type; /* Ensemble Data */
  int typeOfGeneratingProcess;
  long bitsPerValue;
  long editionNumber = 2;
  char name[256];
  grib_handle *gh = NULL;
  gribContainer_t *gc = (gribContainer_t *) gribContainer;
  // extern unsigned char _grib_template_GRIB2[];

  param    = vlistInqVarParam(vlistID, varID);
  datatype = vlistInqVarDatatype(vlistID, varID);
  typeOfGeneratingProcess = vlistInqVarTypeOfGeneratingProcess(vlistID, varID);

  vlistInqVarName(vlistID, varID, name);

#if defined(GRIBAPIENCODETEST)
  gh = (grib_handle *) gribHandleNew(editionNumber);
#else
  gh = gc->gribHandle;
#endif

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  if ( typeOfGeneratingProcess == -1 ) typeOfGeneratingProcess = 0;
  if ( ! gc->init ) GRIB_CHECK(grib_set_long(gh, "typeOfGeneratingProcess", typeOfGeneratingProcess), 0);

  if ( ! gc->init ) gribapiDefInstitut(gh, vlistID, varID);
  if ( ! gc->init ) gribapiDefModel(gh, vlistID, varID);

  if ( ! gc->init ) gribapiDefParam(editionNumber, gh, param, name);
  /*
  if( vlistInqVarEnsemble( vlistID,  varID, &ensID, &ensCount, &forecast_type ) )
    {
      GRIB_CHECK(grib_set_long(gh, "typeOfEnsembleForecast", forecast_type ), 0);
      GRIB_CHECK(grib_set_long(gh, "numberOfForecastsInEnsemble", ensCount ), 0);
      GRIB_CHECK(grib_set_long(gh, "perturbationNumber", ensID ), 0);
    }
  */

  gribapiDefTime(editionNumber, typeOfGeneratingProcess, gh, vdate, vtime, tsteptype, numavg, vlistInqTaxis(vlistID), gc->init);

  if ( editionNumber == 2 && (datatype == DATATYPE_FLT32 || datatype == DATATYPE_FLT64) ) lieee = TRUE;

  /* bitsPerValue have to be defined before call to DefGrid (complex packing) */
  //  if ( lieee == FALSE )
    {
      bitsPerValue = grbBitsPerValue(datatype);
      GRIB_CHECK(grib_set_long(gh, "bitsPerValue", bitsPerValue), 0);
    }

  gribapiDefGrid(editionNumber, gh, gridID, ljpeg, lieee, datatype, nmiss, gc->init);

  gribapiDefLevel(editionNumber, gh, param, zaxisID, levelID, gc->init);

  /* ---------------------------------- */
  /* Local change: 2013-01-28, FP (DWD) */
  /* ---------------------------------- */

  vlist_t *vlistptr;
  vlistptr = vlist_to_pointer(vlistID);
  if (!gc->init) {
    int i;
    for (i=0; i<vlistptr->vars[varID].opt_grib_dbl_nentries; i++)
      {
	int ret = grib_set_double(gh, vlistptr->vars[varID].opt_grib_dbl_keyword[i],
                                  vlistptr->vars[varID].opt_grib_dbl_val[i]);
	if (ret != 0) {
	    fprintf(stderr, "key \"%s\"  :   value = %g\n",
                    vlistptr->vars[varID].opt_grib_dbl_keyword[i],
		    vlistptr->vars[varID].opt_grib_dbl_val[i]);
	}
	GRIB_CHECK(ret, 0);
      }
    for (i=0; i<vlistptr->vars[varID].opt_grib_int_nentries; i++)
      {
	int ret = grib_set_long(gh, vlistptr->vars[varID].opt_grib_int_keyword[i],
	                        vlistptr->vars[varID].opt_grib_int_val[i]);
	if (ret != 0) {
	    fprintf(stderr, "key \"%s\"  :   value = %d\n",
		    vlistptr->vars[varID].opt_grib_int_keyword[i],
		    vlistptr->vars[varID].opt_grib_int_val[i]);
	}
	GRIB_CHECK(ret, 0);
      }
  }


  if ( nmiss > 0 )
    {
      GRIB_CHECK(grib_set_long(gh, "bitmapPresent", 1), 0);
      GRIB_CHECK(grib_set_double(gh, "missingValue", vlistInqVarMissval(vlistID, varID)), 0);
    }

  GRIB_CHECK(grib_set_double_array(gh, "values", data, datasize), 0);

  /* get the size of coded message  */
  GRIB_CHECK(grib_get_message(gh, (const void **)&dummy, &recsize), 0);
  recsize += 512; /* add some space for possible filling */
  *gribbuffersize = recsize;
  *gribbuffer = (unsigned char *) malloc(*gribbuffersize);

  /* get a copy of the coded message */
  GRIB_CHECK(grib_get_message_copy(gh, *gribbuffer, &recsize), 0);

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
