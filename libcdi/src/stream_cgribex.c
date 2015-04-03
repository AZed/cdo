#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <float.h>  /* FLT_EPSILON */

#include "dmemory.h"
#include "cdi.h"
#include "stream_int.h"
#include "file.h"
#include "varscan.h"
#include "datetime.h"
#include "vlist.h"
#include "stream_grb.h"

#if  defined  (HAVE_LIBCGRIBEX)
#  include "cgribex.h"
#endif

extern int cdiInventoryMode;

typedef struct {
  int param;
  int level1;
  int level2;
  int ltype;
} compvar_t; 


#if  defined  (HAVE_LIBCGRIBEX)
static
int cgribexGetGridType(int *isec2)
{
  int gridtype = GRID_GENERIC;

  switch (ISEC2_GridType)
    {
    case  GRIB1_GTYPE_LATLON:     { if ( ISEC2_Reduced )      break; }
    case  GRIB1_GTYPE_LATLON_ROT: { gridtype = GRID_LONLAT;   break; }
    case  GRIB1_GTYPE_LCC:        { gridtype = GRID_LCC;      break; }
    case  GRIB1_GTYPE_GAUSSIAN:   { if ( ISEC2_Reduced )
	                              gridtype = GRID_GAUSSIAN_REDUCED;
                         	    else
				      gridtype = GRID_GAUSSIAN;
          	                    break;
                                  }
    case  GRIB1_GTYPE_SPECTRAL:   { gridtype = GRID_SPECTRAL; break; }
    case  GRIB1_GTYPE_GME:        { gridtype = GRID_GME;      break; }
    }

  return (gridtype);
}
#endif

#if  defined  (HAVE_LIBCGRIBEX)
static
int cgribexGetIsRotated(int *isec2)
{
  int isRotated = 0;

  if ( ISEC2_GridType == GRIB1_GTYPE_LATLON_ROT )
    {
      isRotated = 1;
    }

  return (isRotated);
}
#endif

#if  defined  (HAVE_LIBCGRIBEX)
static
int cgribexGetZaxisHasBounds(int grb_ltype)
{
  int lbounds = 0;

  switch (grb_ltype)
    {
    case GRIB1_LTYPE_SIGMA_LAYER:
    case GRIB1_LTYPE_HYBRID_LAYER:
    case GRIB1_LTYPE_LANDDEPTH_LAYER:
      {
	lbounds = 1;
	break;
      }
    }

  return (lbounds);
}
#endif

#if  defined  (HAVE_LIBCGRIBEX)
static
int cgribexGetTimeUnit(int *isec1)
{
  int timeunit = TUNIT_HOUR;
  static int lprint = TRUE;

  switch ( ISEC1_TimeUnit )
    {
    case ISEC1_TABLE4_MINUTE:  timeunit = TUNIT_MINUTE;  break;
    case ISEC1_TABLE4_QUARTER: timeunit = TUNIT_QUARTER; break;
    case ISEC1_TABLE4_HOUR:    timeunit = TUNIT_HOUR;    break;
    case ISEC1_TABLE4_3HOURS:  timeunit = TUNIT_3HOURS;  break;
    case ISEC1_TABLE4_6HOURS:  timeunit = TUNIT_6HOURS;  break;
    case ISEC1_TABLE4_12HOURS: timeunit = TUNIT_12HOURS; break;
    case ISEC1_TABLE4_DAY:     timeunit = TUNIT_DAY;     break;
    default:
      if ( lprint )
	{
	  Message("GRIB time unit %d unsupported!", ISEC1_TimeUnit);
	  lprint = FALSE;
	}
    }

  return (timeunit);
}
#endif

#if  defined  (HAVE_LIBCGRIBEX)
static
int cgribexTimeIsFC(int *isec1)
{
  int isFC = TRUE;

  if ( ISEC1_TimeRange == 10 && ISEC1_TimePeriod1 == 0 && ISEC1_TimePeriod2 == 0 )
    isFC = FALSE;

  return (isFC);
}
#endif

#if  defined  (HAVE_LIBCGRIBEX)
static
int cgribexGetTsteptype(int timerange)
{
  int tsteptype = 0;
  static int lprint = TRUE;

  switch ( timerange )
    {
    case  0:  tsteptype = TSTEP_INSTANT;  break;
    case  1:  tsteptype = TSTEP_INSTANT2; break;
    case  2:  tsteptype = TSTEP_RANGE;    break;
    case  3:  tsteptype = TSTEP_AVG;      break;
    case  4:  tsteptype = TSTEP_ACCUM;    break;
    case  5:  tsteptype = TSTEP_DIFF;     break;
    case 10:  tsteptype = TSTEP_INSTANT3; break;
    default:
      if ( lprint )
	{
	  Message("GRIB time range %d unsupported!", timerange);
	  lprint = FALSE;
	}
    }

  return (tsteptype);
}
#endif

#if  defined  (HAVE_LIBCGRIBEX)
static
void cgribexGetGrid(stream_t *streamptr, int *isec2, int *isec4, grid_t *grid)
{
  int gridtype;

  gridtype = cgribexGetGridType(isec2);

  if ( streamptr->unreduced && gridtype == GRID_GAUSSIAN_REDUCED )
    {
      gridtype = GRID_GAUSSIAN;
      ISEC2_NumLon = 2*ISEC2_NumLat;
      ISEC4_NumValues = ISEC2_NumLon*ISEC2_NumLat;
    }

  memset(grid, 0, sizeof(grid_t));
  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
      {
	if ( ISEC4_NumValues != ISEC2_NumLon*ISEC2_NumLat )
	  Error("numberOfPoints (%d) and gridSize (%d) differ!",
		ISEC4_NumValues, ISEC2_NumLon*ISEC2_NumLat);
	grid->size  = ISEC4_NumValues;
	grid->xsize = ISEC2_NumLon;
	grid->ysize = ISEC2_NumLat;
        if ( gridtype == GRID_GAUSSIAN ) grid->np = ISEC2_NumPar;
	grid->xinc  = 0;
	grid->yinc  = 0;
	grid->xdef  = 0;
	/* if ( ISEC2_FirstLon != 0 || ISEC2_LastLon != 0 ) */
	  {
	    if ( grid->xsize > 1 )
	      {
		if ( ISEC2_ResFlag && ISEC2_LonIncr > 0 )
		  if ( ISEC2_LonIncr*(grid->xsize-1) > 360000 ) ISEC2_LonIncr = 0;

		if ( ISEC2_ResFlag && ISEC2_LonIncr > 0 )
		  grid->xinc = ISEC2_LonIncr * 0.001;
		else
		  grid->xinc = (ISEC2_LastLon - ISEC2_FirstLon) * 0.001 / (grid->xsize - 1);

		/* correct xinc if necessary */
		if ( ISEC2_FirstLon == 0 && ISEC2_LastLon > 354000 )
		  {
		    double xinc = 360. / grid->xsize;

		    if ( fabs(grid->xinc-xinc) > 0.0 )
		      {
			grid->xinc = xinc;
			if ( CDI_Debug ) Message("set xinc to %g", grid->xinc);
		      }
		  }
	      }
	    grid->xfirst = ISEC2_FirstLon * 0.001;
	    grid->xlast  = ISEC2_LastLon  * 0.001;
	    grid->xdef   = 2;
	  }
	grid->ydef  = 0;
	/* if ( ISEC2_FirstLat != 0 || ISEC2_LastLat != 0 ) */
	  {
	    if ( grid->ysize > 1 )
	      {
		if ( ISEC2_ResFlag && ISEC2_LatIncr > 0 )
		  if ( ISEC2_LatIncr*(grid->ysize-1) > 180000 ) ISEC2_LatIncr = 0;

		if ( ISEC2_ResFlag && ISEC2_LatIncr > 0 )
		  grid->yinc = ISEC2_LatIncr * 0.001;
		else
		  grid->yinc = (ISEC2_LastLat - ISEC2_FirstLat) * 0.001 / (grid->ysize - 1);
	      }
	    grid->yfirst = ISEC2_FirstLat * 0.001;
	    grid->ylast  = ISEC2_LastLat  * 0.001;
	    grid->ydef   = 2;
	  }
	break;
      }
    case GRID_GAUSSIAN_REDUCED:
      {
        grid->np     = ISEC2_NumPar;
	grid->size   = ISEC4_NumValues;
        grid->rowlon = ISEC2_RowLonPtr;
	grid->ysize  = ISEC2_NumLat;
	grid->xinc   = 0;
	grid->yinc   = 0;
	grid->xdef   = 0;
	/* if ( ISEC2_FirstLon != 0 || ISEC2_LastLon != 0 ) */
	  {
	    if ( grid->xsize > 1 )
	      {
		if ( ISEC2_ResFlag && ISEC2_LonIncr > 0 )
		  grid->xinc = ISEC2_LonIncr * 0.001;
		else
		  grid->xinc = (ISEC2_LastLon - ISEC2_FirstLon) * 0.001 / (grid->xsize - 1);
	      }
	    grid->xfirst = ISEC2_FirstLon * 0.001;
	    grid->xlast  = ISEC2_LastLon  * 0.001;
	    grid->xdef   = 2;
	  }
	grid->ydef  = 0;
	/* if ( ISEC2_FirstLat != 0 || ISEC2_LastLat != 0 ) */
	  {
	    if ( grid->ysize > 1 )
	      {
		if ( ISEC2_ResFlag && ISEC2_LatIncr > 0 )
		  grid->yinc = ISEC2_LatIncr * 0.001;
		else
		  grid->yinc = (ISEC2_LastLat - ISEC2_FirstLat) * 0.001 / (grid->ysize - 1);
	      }
	    grid->yfirst = ISEC2_FirstLat * 0.001;
	    grid->ylast  = ISEC2_LastLat  * 0.001;
	    grid->ydef   = 2;
	  }
	break;
      }
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
    case GRID_SPECTRAL:
      {
	grid->size  = ISEC4_NumValues;
	grid->trunc = ISEC2_PentaJ;
	if ( ISEC2_RepMode == 2 )
	  grid->lcomplex = 1;
	else
	  grid->lcomplex = 0;

	break;
      }
    case GRID_GME:
      {
	grid->size  = ISEC4_NumValues;
	grid->nd    = ISEC2_GME_ND;
	grid->ni    = ISEC2_GME_NI;
	grid->ni2   = ISEC2_GME_NI2;
	grid->ni3   = ISEC2_GME_NI3;
	break;
      }
    case GRID_GENERIC:
      {
	grid->size  = ISEC4_NumValues;
	grid->xsize = 0;
	grid->ysize = 0;
	break;
      }
    default:
      {
	Error("Unsupported grid type: %s", gridNamePtr(gridtype));
	break;
      }
    }

  grid->isRotated = FALSE;
  if ( cgribexGetIsRotated(isec2) )
    {
      grid->isRotated = TRUE;
      grid->ypole     = - ISEC2_LatSP * 0.001;
      grid->xpole     =   ISEC2_LonSP * 0.001 - 180;
      grid->angle     = 0;
    }

  grid->xvals = NULL;
  grid->yvals = NULL;
  grid->type  = gridtype;
}
#endif

#if  defined  (HAVE_LIBCGRIBEX)
static
void cgribexAddRecord(int streamID, int param, int *isec1, int *isec2, double *fsec2, double *fsec3,
		      int *isec4, long recsize, off_t position, int datatype, int comptype, int lmv)
{
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

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  tsID    = streamptr->curTsID;
  recID   = recordNewEntry(streamID, tsID);
  record  = &streamptr->tsteps[tsID].records[recID];

  tsteptype = cgribexGetTsteptype(ISEC1_TimeRange);
  numavg    = ISEC1_AvgNum;

  level1  = ISEC1_Level1;
  level2  = ISEC1_Level2;

  /* fprintf(stderr, "param %d %d %d %d\n", param, level1, level2, ISEC1_LevelType); */

  (*record).size     = recsize;
  (*record).position = position;
  (*record).param    = param;
  (*record).ilevel   = level1;
  (*record).ilevel2  = level2;
  (*record).ltype    = ISEC1_LevelType;

  cgribexGetGrid(streamptr, isec2, isec4, &grid);

  gridID = varDefGrid(vlistID, grid, 0);

  zaxistype = grib1ltypeToZaxisType(ISEC1_LevelType);

  if ( zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF )
    {
      int vctsize = ISEC2_NumVCP;
      double *vctptr = &fsec2[10];

      varDefVCT(vctsize, vctptr);
    }

  lbounds = cgribexGetZaxisHasBounds(ISEC1_LevelType);

  if ( datatype > 32 ) datatype = DATATYPE_PACK32;
  if ( datatype <  0 ) datatype = DATATYPE_PACK;

  varAddRecord(recID, param, gridID, zaxistype, lbounds, level1, level2,
	       datatype, &varID, &levelID, tsteptype, numavg, ISEC1_LevelType, NULL, NULL, NULL);

  (*record).varID   = varID;
  (*record).levelID = levelID;

  varDefCompType(varID, comptype);

  if ( lmv ) varDefMissval(varID, FSEC3_MissVal);

  if ( varInqInst(varID) == CDI_UNDEFID )
    {
      int center, subcenter, instID;
      center    = ISEC1_CenterID;
      subcenter = ISEC1_SubCenterID;
      instID    = institutInq(center, subcenter, NULL, NULL);
      if ( instID == CDI_UNDEFID )
	instID = institutDef(center, subcenter, NULL, NULL);
      varDefInst(varID, instID);
    }

  if ( varInqModel(varID) == CDI_UNDEFID )
    {
      int modelID;
      modelID = modelInq(varInqInst(varID), ISEC1_ModelID, NULL);
      if ( modelID == CDI_UNDEFID )
	modelID = modelDef(varInqInst(varID), ISEC1_ModelID, NULL);
      varDefModel(varID, modelID);
    }

  if ( varInqTable(varID) == CDI_UNDEFID )
    {
      int tableID;

      tableID = tableInq(varInqModel(varID), ISEC1_CodeTable, NULL);

      if ( tableID == CDI_UNDEFID )
	tableID = tableDef(varInqModel(varID), ISEC1_CodeTable, NULL);
      varDefTable(varID, tableID);
    }

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;

  if ( CDI_Debug )
    Message("varID = %d  param = %d  zaxistype = %d  gridID = %d  levelID = %d",
	    varID, param, zaxistype, gridID, levelID);
}
#endif

#if  defined  (HAVE_LIBCGRIBEX)
static
void MCH_get_undef(int *isec1, double *undef_pds, double *undef_eps)
{
  /* 2010-01-13: Oliver Fuhrer */
  if ( ISEC1_CenterID == 215 ) {
    if (isec1[34] != 0 && isec1[34] != 255) {
      if (isec1[34] & 2) {
        if (isec1[34] & 1) {
          *undef_pds = -0.99*pow(10.0,-isec1[35]);
        } else {
          *undef_pds = +0.99*pow(10.0,-isec1[35]);
        }
        *undef_eps = pow(10.0,-isec1[35]-1);
      } else {
        if (isec1[34] & 1) {
          *undef_pds = -0.99*pow(10.0,+isec1[35]);
        } else {
          *undef_pds = +0.99*pow(10.0,+isec1[35]);
        }
        *undef_eps = pow(10.0,isec1[35]-1);
      }
    }
  }
}
#endif

#if  defined  (HAVE_LIBCGRIBEX)
static
void cgribexDecodeHeader(int *isec0, int *isec1, int *isec2, double *fsec2,
			 int *isec3, double *fsec3, int *isec4, double *fsec4, 
			 int *gribbuffer, int recsize, int *lmv)
{
  int iret = 0, ipunp = 0, iword = 0;

  gribExDP(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4,
	   ipunp, (int *) gribbuffer, recsize, &iword, "J", &iret);

  *lmv = 0;

  if ( ISEC1_CenterID == 215 && (isec1[34] != 0 && isec1[34] != 255) )
    {
      double undef_pds, undef_eps;

      MCH_get_undef(isec1, &undef_pds, &undef_eps);
      FSEC3_MissVal = undef_pds;
      *lmv = 1;
    }
}
#endif

int cgribexScanTimestep1(int streamID)
{
#if  defined  (HAVE_LIBCGRIBEX)
  int *isec0, *isec1, *isec2, *isec3, *isec4;
  double fsec2[512], fsec3[2], *fsec4 = NULL;
  int lmv = 0;
  off_t recpos = 0;
  unsigned char *gribbuffer = NULL;
  long buffersize = 0;
  int rstatus;
  int fileID;
  int param = 0;
  int level1 = 0, level2 = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0;
  int tsID;
  int varID;
  size_t readsize;
  int nrecords, nrecs, recID;
  int datatype;
  long recsize = 0;
  int warn_time = TRUE;
  int warn_numavg = TRUE;
  int taxisID = -1;
  int rdate = 0, rtime = 0, tunit = 0, fcast = 0;
  taxis_t *taxis;
  int vlistID;
  int comptype;
  long unzipsize;
  compvar_t compVar, compVar0;
  stream_t *streamptr;
  extern int cdiSkipRecords;
  int nskip = cdiSkipRecords;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  streamptr->curTsID = 0;

  isec0 = streamptr->record->sec0;
  isec1 = streamptr->record->sec1;
  isec2 = streamptr->record->sec2;
  isec3 = streamptr->record->sec3;
  isec4 = streamptr->record->sec4;

  tsID  = tstepsNewEntry(streamID);
  taxis = &streamptr->tsteps[tsID].taxis;

  if ( tsID != 0 )
    Error("Internal problem! tstepsNewEntry returns %d", tsID);

  fileID = streamInqFileID(streamID);

  while ( nskip-- > 0 )
    {
      recsize = gribGetSize(fileID);
      if ( recsize == 0 )
	Error("Skipping of %d records failed!", cdiSkipRecords);

      recpos  = fileGetPos(fileID);
      fileSetPos(fileID, recsize, SEEK_CUR);
    }

  nrecs = 0;
  while ( TRUE )
    {
      recsize = gribGetSize(fileID);
      recpos  = fileGetPos(fileID);

      if ( recsize == 0 )
	{
	  if ( nrecs == 0 )
	    Error("No GRIB records found!");

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
	  unzipsize += 100; /* need 0 to 1 bytes for rounding of bds */
	  if ( (long) buffersize < unzipsize )
	    {
	      buffersize = unzipsize;
	      gribbuffer = (unsigned char *) realloc(gribbuffer, buffersize);
	    }
	}

      cgribexDecodeHeader(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4,
			  (int *) gribbuffer, recsize, &lmv);

      param = cdiEncodeParam(ISEC1_Parameter, ISEC1_CodeTable, 255);
      if ( ISEC1_LevelType == 100 ) ISEC1_Level1 *= 100;
      if ( ISEC1_LevelType ==  99 ) ISEC1_LevelType = 100;
      level1   = ISEC1_Level1;
      level2   = ISEC1_Level2;

      gribDateTime(isec1, &vdate, &vtime);

      if ( ISEC4_NumBits > 0 && ISEC4_NumBits <= 32 )
	datatype = ISEC4_NumBits;
      else
        datatype = DATATYPE_PACK;

      if ( nrecs == 0 )
	{
	  datetime0.date = vdate;
	  datetime0.time = vtime;
	  rdate = gribRefDate(isec1);
	  rtime = gribRefTime(isec1);
	  tunit = cgribexGetTimeUnit(isec1);
	  fcast = cgribexTimeIsFC(isec1);
	}
      else
	{
	  datetime.date  = vdate;
	  datetime.time  = vtime;
	  compVar.param  = param;
          compVar.level1 = level1;
          compVar.level2 = level2;
          compVar.ltype  = ISEC1_LevelType;
	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      compVar0.param  = streamptr->tsteps[0].records[recID].param;
	      compVar0.level1 = streamptr->tsteps[0].records[recID].ilevel;
	      compVar0.level2 = streamptr->tsteps[0].records[recID].ilevel2;
	      compVar0.ltype  = streamptr->tsteps[0].records[recID].ltype;

	      if ( memcmp(&compVar0, &compVar, sizeof(compvar_t)) == 0 ) break;
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

      if ( ISEC1_AvgNum )
	{
	  if (  taxis->numavg && warn_numavg && (taxis->numavg != ISEC1_AvgNum) )
	    {
	      Warning("Changing numavg from %d to %d not supported!",
		      taxis->numavg, ISEC1_AvgNum);
	      warn_numavg = FALSE;
	    }
	  else
	    {
	      taxis->numavg = ISEC1_AvgNum;
	    }
	}

      nrecs++;

      if ( CDI_Debug )
	Message("%4d %8d %4d  %8d %8d %6d", nrecs, (int)recpos, param, level1, vdate, vtime);

      cgribexAddRecord(streamID, param, isec1, isec2, fsec2, fsec3,
		       isec4, recsize, recpos, datatype, comptype, lmv);
    }

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
      taxis->unit  = tunit;
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
  Error("CGRIBEX support not compiled in!");
#endif

  return (0);
}


int cgribexScanTimestep2(int streamID)
{
  int rstatus = 0;
#if  defined  (HAVE_LIBCGRIBEX)
  int *isec0, *isec1, *isec2, *isec3, *isec4;
  double fsec2[512], fsec3[2], *fsec4 = NULL;
  int lmv = 0;
  off_t recpos = 0;
  unsigned char *gribbuffer = NULL;
  long buffersize = 0;
  int fileID;
  int param = 0;
  int level1 = 0, level2 = 0, vdate = 0, vtime = 0;
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
  compvar_t compVar, compVar0;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  streamptr->curTsID = 1;

  isec0 = streamptr->record->sec0;
  isec1 = streamptr->record->sec1;
  isec2 = streamptr->record->sec2;
  isec3 = streamptr->record->sec3;
  isec4 = streamptr->record->sec4;

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

      cgribexDecodeHeader(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4,
			  (int *) gribbuffer, recsize, &lmv);

      param = cdiEncodeParam(ISEC1_Parameter, ISEC1_CodeTable, 255);
      if ( ISEC1_LevelType == 100 ) ISEC1_Level1 *= 100;
      if ( ISEC1_LevelType ==  99 ) ISEC1_LevelType = 100;
      level1    = ISEC1_Level1;
      level2    = ISEC1_Level2;

      gribDateTime(isec1, &vdate, &vtime);

      if ( rindex == 0 )
	{
	  if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
	    {
	      taxis->type  = TAXIS_RELATIVE;
	      taxis->rdate = gribRefDate(isec1);
	      taxis->rtime = gribRefTime(isec1);
	    }
	  else
	    {
	      taxis->type  = TAXIS_ABSOLUTE;
	    }
	  taxis->unit  = cgribexGetTimeUnit(isec1);
	  taxis->vdate = vdate;
	  taxis->vtime = vtime;

	  datetime0.date = vdate;
	  datetime0.time = vtime;
	}

      tsteptype = cgribexGetTsteptype(ISEC1_TimeRange);

      if ( ISEC1_AvgNum )
	{
	  if (  taxis->numavg && warn_numavg &&
		(taxis->numavg != ISEC1_AvgNum) )
	    {
	  /*
	      Warning("Changing numavg from %d to %d not supported!",
		      taxis->numavg, ISEC1_AvgNum);
	  */
	      warn_numavg = FALSE;
	    }
	  else
	    {
	      taxis->numavg = ISEC1_AvgNum;
	    }
	}

      datetime.date  = vdate;
      datetime.time  = vtime;
      compVar.param  = param;
      compVar.level1 = level1;
      compVar.level2 = level2;
      compVar.ltype  = ISEC1_LevelType;
      for ( recID = 0; recID < nrecords; recID++ )
	{
	  compVar0.param  = streamptr->tsteps[tsID].records[recID].param;
	  compVar0.level1 = streamptr->tsteps[tsID].records[recID].ilevel;
	  compVar0.level2 = streamptr->tsteps[tsID].records[recID].ilevel2;
	  compVar0.ltype  = streamptr->tsteps[tsID].records[recID].ltype;

	  if ( memcmp(&compVar0, &compVar, sizeof(compvar_t)) == 0 ) break;
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

      if ( memcmp(&compVar0, &compVar, sizeof(compvar_t)) != 0 )
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
      if ( gridInqSize(gridID) == 1 && gridInqType(gridID) == GRID_LONLAT )
	{
	  if ( IS_NOT_EQUAL(gridInqXval(gridID, 0),ISEC2_FirstLon*0.001) ||
	       IS_NOT_EQUAL(gridInqYval(gridID, 0),ISEC2_FirstLat*0.001) )
	    gridChangeType(gridID, GRID_TRAJECTORY);
	}

      if ( tsteptype != vlistInqVarTsteptype(vlistID, varID) )
	vlistDefVarTsteptype(vlistID, varID, tsteptype);

      rindex++;
    }

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


int cgribexScanTimestep(int streamID)
{
  int rstatus = 0;
#if  defined  (HAVE_LIBCGRIBEX)
  int *isec0, *isec1, *isec2, *isec3, *isec4;
  double fsec2[512], fsec3[2], *fsec4 = NULL;
  int lmv = 0;
  long recsize = 0;
  off_t recpos = 0;
  unsigned char *gribbuffer;
  long buffersize = 0;
  int fileID;
  int param = 0;
  int level1 = 0, level2 = 0, vdate = 0, vtime = 0;
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
  compvar_t compVar, compVar0;
  stream_t *streamptr;

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

  isec0 = streamptr->record->sec0;
  isec1 = streamptr->record->sec1;
  isec2 = streamptr->record->sec2;
  isec3 = streamptr->record->sec3;
  isec4 = streamptr->record->sec4;

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
	  if ( recsize > buffersize )
	    {
	      buffersize = recsize;
	      gribbuffer = (unsigned char *) realloc(gribbuffer, buffersize);
	    }

	  if ( rindex >= nrecs ) break;

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

	  cgribexDecodeHeader(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4,
			      (int *) gribbuffer, recsize, &lmv);

	  param = cdiEncodeParam(ISEC1_Parameter, ISEC1_CodeTable, 255);
	  if ( ISEC1_LevelType == 100 ) ISEC1_Level1 *= 100;
	  if ( ISEC1_LevelType ==  99 ) ISEC1_LevelType = 100;
	  level1   = ISEC1_Level1;
	  level2   = ISEC1_Level2;

	  gribDateTime(isec1, &vdate, &vtime);

	  if ( rindex == nrecs ) break;

	  if ( rindex == 0 )
	    {
	      taxisID = vlistInqTaxis(vlistID);
	      if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
		{
		  taxis->type  = TAXIS_RELATIVE;
		  taxis->rdate = gribRefDate(isec1);
		  taxis->rtime = gribRefTime(isec1);
		}
	      else
		{
		  taxis->type  = TAXIS_ABSOLUTE;
		}
	      taxis->unit  = cgribexGetTimeUnit(isec1);
	      taxis->vdate = vdate;
	      taxis->vtime = vtime;

	      datetime0.date = vdate;
	      datetime0.time = vtime;
	    }

	  if ( ISEC1_AvgNum )
	    {
	      if (  taxis->numavg && warn_numavg &&
		   (taxis->numavg != ISEC1_AvgNum) )
		{
	      /*
	          Warning("Changing numavg from %d to %d not supported!",
			  streamptr->tsteps[tsID].taxis.numavg, ISEC1_AvgNum);
	      */
		  warn_numavg = FALSE;
		}
	      else
		{
		  taxis->numavg = ISEC1_AvgNum;
		}
	    }
	  
	  datetime.date  = vdate;
	  datetime.time  = vtime;
	  compVar.param  = param;
          compVar.level1 = level1;
          compVar.level2 = level2;
          compVar.ltype  = ISEC1_LevelType;
	  for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	    {
	      recID   = streamptr->tsteps[1].recIDs[vrecID];
	      compVar0.param  = streamptr->tsteps[tsID].records[recID].param;
	      compVar0.level1 = streamptr->tsteps[tsID].records[recID].ilevel;
	      compVar0.level2 = streamptr->tsteps[tsID].records[recID].ilevel2;
	      compVar0.ltype  = streamptr->tsteps[tsID].records[recID].ltype;

	      if ( memcmp(&compVar0, &compVar, sizeof(compvar_t)) == 0 ) break;
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

	  if ( memcmp(&compVar0, &compVar, sizeof(compvar_t)) != 0 )
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

	  rindex++;
	}

      for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	{
	  recID   = streamptr->tsteps[tsID].recIDs[vrecID];
	  if ( ! streamptr->tsteps[tsID].records[recID].used ) break;
	}

      if ( vrecID < nrecs )
	{
	  char paramstr[32];
	  cdiParamToString(streamptr->tsteps[tsID].records[recID].param, paramstr, sizeof(paramstr));
	  Warning("Param=%s level=%d not found at timestep %d!",
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
#endif

  return (rstatus);
}


int cgribexDecode(unsigned char *gribbuffer, int gribsize, double *data, int gridsize,
		  int unreduced, int *nmiss, int *zip, double missval)
{
  int status = 0;
#if  defined  (HAVE_LIBCGRIBEX)
  int iret = 0, iword = 0;
  int isec0[2], isec1[4096], isec2[4096], isec3[2], isec4[512];
  int izip;
  long unzipsize;
  double fsec2[512], fsec3[2];
  char hoper[2];

  *zip = 0;

  if ( unreduced ) strcpy(hoper, "R");
  else             strcpy(hoper, "D");

  FSEC3_MissVal = missval;

  if ( (izip = gribGetZip(gribsize, gribbuffer, &unzipsize)) > 0 )
    {
      *zip = izip;
      if ( izip == 128 ) /* szip */
	{
	  unsigned char *itmpbuffer = NULL;
	  size_t itmpbuffersize = 0;

	  if ( unzipsize < (long) gribsize )
	    {
	      fprintf(stderr, "Decompressed size smaller than compressed size (in %d; out %ld)!\n",
		      gribsize, unzipsize);
	      return (status);
	    }

	  if ( itmpbuffersize < (size_t) gribsize )
	    {
	      itmpbuffersize = gribsize;
	      itmpbuffer = (unsigned char *) realloc(itmpbuffer, itmpbuffersize);
	    }

	  memcpy(itmpbuffer, gribbuffer, itmpbuffersize);

	  unzipsize += 100; /* need 0 to 1 bytes for rounding of bds */

	  gribsize = gribUnzip(gribbuffer, unzipsize, itmpbuffer, gribsize);

	  if ( gribsize <= 0 )
	    Error("Decompression problem!");
	  
	  free(itmpbuffer);
	}
      else
	{
	  Error("Decompression for %d not implemented!", izip);
	}
    }

  gribExDP(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, data,
	   gridsize, (int *) gribbuffer, gribsize, &iword, hoper, &iret);

  if ( ISEC1_Sec2Or3Flag & 64 )
    *nmiss = ISEC4_NumValues - ISEC4_NumNonMissValues;
  else
    *nmiss = 0;

  if ( ISEC1_CenterID == 215 && (isec1[34] != 0 && isec1[34] != 255) )
    {
      double undef_pds, undef_eps;
      int i;

      MCH_get_undef(isec1, &undef_pds, &undef_eps);

      *nmiss = 0;
      for ( i = 0; i < gridsize; i++ )
        if ( (abs(data[i]-undef_pds) < undef_eps) || IS_EQUAL(data[i],FSEC3_MissVal) ) {
          data[i] = missval;
          (*nmiss)++;
        }
    }
#else
  Error("CGRIBEX support not compiled in!");
#endif

  return (status);
}


#if  defined  (HAVE_LIBCGRIBEX)
static
void cgribexDefInstitut(int *isec1, int vlistID, int varID)
{
  int instID;

  if ( vlistInqInstitut(vlistID) != CDI_UNDEFID )
    instID = vlistInqInstitut(vlistID);
  else
    instID = vlistInqVarInstitut(vlistID, varID);

  if ( instID != CDI_UNDEFID )
    {
      int center, subcenter;
      center    = institutInqCenter(instID);
      subcenter = institutInqSubcenter(instID);
      ISEC1_CenterID    = center;
      ISEC1_SubCenterID = subcenter;
    }
}

static
void cgribexDefModel(int *isec1, int vlistID, int varID)
{
  int modelID;

  if ( vlistInqModel(vlistID) != CDI_UNDEFID )
    modelID = vlistInqModel(vlistID);
  else
    modelID = vlistInqVarModel(vlistID, varID);

  if ( modelID != CDI_UNDEFID )
    ISEC1_ModelID = modelInqGribID(modelID);
}

static
void cgribexDefParam(int *isec1, int param)
{
  int pdis, pcat, pnum;
  static int lwarn = TRUE;

  cdiDecodeParam(param, &pnum, &pcat, &pdis);

  if ( pdis != 255 && lwarn )
    {
      char paramstr[32];
      cdiParamToString(param, paramstr, sizeof(paramstr));
      Warning("Can not convert GRIB2 parameter (%s) to GRIB1!", paramstr);
      lwarn = FALSE;
    }

  if ( pnum < 0 ) pnum = -pnum;

  ISEC1_CodeTable = pcat;
  ISEC1_Parameter = pnum;
}

static
int cgribexDefTimerange(int tsteptype, int factor, int calendar,
			int rdate, int rtime, int vdate, int vtime, int *pip1, int *pip2)
{
  int timerange = -1;
  int year, month, day, hour, minute, second;
  int julday1, secofday1, julday2, secofday2, days, secs;
  int ip, ip1 = 0, ip2 = 0;

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
	case TSTEP_INSTANT:  timerange =  0; ip1 = ip; ip2 = 0;  break;
	case TSTEP_INSTANT2: timerange =  1; ip1 = 0;  ip2 = 0;  break;
	case TSTEP_RANGE:    timerange =  2; ip1 = 0;  ip2 = ip; break;
	case TSTEP_AVG:      timerange =  3; ip1 = 0;  ip2 = ip; break;
	case TSTEP_ACCUM:    timerange =  4; ip1 = 0;  ip2 = ip; break;
	case TSTEP_DIFF:     timerange =  5; ip1 = 0;  ip2 = ip; break;
	case TSTEP_INSTANT3:
	default:             timerange = 10; ip1 = ip/256; ip2 = ip%256;  break;
	}
    }

  *pip1 = ip1;
  *pip2 = ip2;

  return (timerange);
}

static
int cgribexDefDateTime(int *isec1, int timeunit, int date, int time)
{
  int year, month, day, hour, minute, second;
  int century = 0;
  int factor = 1;

  cdiDecodeDate(date, &year, &month, &day);
  cdiDecodeTime(time, &hour, &minute, &second);

  century =  year / 100;

  ISEC1_Year = year - century*100;
  
  if ( year < 0 )
    {
      century = -century;
      ISEC1_Year = -ISEC1_Year;
    }

  if ( ISEC1_Year == 0 )
    {
      century -= 1;
      ISEC1_Year = 100;
    }

  century += 1;
  if ( year < 0 ) century = -century;

  ISEC1_Month  = month;
  ISEC1_Day    = day;
  ISEC1_Hour   = hour;
  ISEC1_Minute = minute;

  ISEC1_Century = century;

  switch (timeunit)
    {
    case TUNIT_MINUTE:  factor =    60; ISEC1_TimeUnit = ISEC1_TABLE4_MINUTE;  break;
    case TUNIT_QUARTER: factor =   900; ISEC1_TimeUnit = ISEC1_TABLE4_QUARTER; break;
    case TUNIT_HOUR:    factor =  3600; ISEC1_TimeUnit = ISEC1_TABLE4_HOUR;    break;
    case TUNIT_3HOURS:  factor = 10800; ISEC1_TimeUnit = ISEC1_TABLE4_3HOURS;  break;
    case TUNIT_6HOURS:  factor = 21600; ISEC1_TimeUnit = ISEC1_TABLE4_6HOURS;  break;
    case TUNIT_12HOURS: factor = 43200; ISEC1_TimeUnit = ISEC1_TABLE4_12HOURS; break;
    case TUNIT_DAY:     factor = 86400; ISEC1_TimeUnit = ISEC1_TABLE4_DAY;     break;
    default:            factor =  3600; ISEC1_TimeUnit = ISEC1_TABLE4_HOUR;    break;
    }

  return (factor);
}

static
void cgribexDefTime(int *isec1, int vdate, int vtime, int tsteptype, int numavg, int taxisID)
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
      int ip1 = 0, ip2 = 0;
      int calendar;

      calendar = taxisInqCalendar(taxisID);
      rdate    = taxisInqRdate(taxisID);
      rtime    = taxisInqRtime(taxisID);

      factor = cgribexDefDateTime(isec1, timeunit, rdate, rtime);

      timerange = cgribexDefTimerange(tsteptype, factor, calendar,
				      rdate, rtime, vdate, vtime, &ip1, &ip2);

      if ( timerange == -1 || timerange == 3 )
	{
	  timetype = TAXIS_ABSOLUTE;
	}
      /*
      else if ( timerange == 10 )
	{
	  if ( ip1 < 0 || ip1 > 0xFFFF ) timetype = TAXIS_ABSOLUTE;
	  if ( ip2 < 0 || ip2 > 0xFFFF ) timetype = TAXIS_ABSOLUTE;
	}
      */
      else
	{
	  if ( ip1 < 0 || ip1 > 0xFF   ) timetype = TAXIS_ABSOLUTE;
	  if ( ip2 < 0 || ip2 > 0xFF   ) timetype = TAXIS_ABSOLUTE;
	}

      ISEC1_TimeRange   = timerange;
      ISEC1_TimePeriod1 = ip1;
      ISEC1_TimePeriod2 = ip2;
    }

  if ( timetype == TAXIS_ABSOLUTE )
    {
      (void) cgribexDefDateTime(isec1, timeunit, vdate, vtime);

      /*
      if ( numavg > 0 )
	ISEC1_TimeRange = 0;
      else
      */
      if ( ISEC1_TimeRange != 3 )
	ISEC1_TimeRange   = 10;

      ISEC1_TimePeriod1 = 0;
      ISEC1_TimePeriod2 = 0;
    }

  ISEC1_AvgNum         = numavg;
  ISEC1_AvgMiss        = 0;
  ISEC1_DecScaleFactor = 0;
}

static
void cgribexDefGrid(int *isec1, int *isec2, int *isec4, int gridID)
{
  int gridtype;
  int lcurvi = FALSE;
  static short lwarn = TRUE;

  memset(isec2, 0, 16*sizeof(int));

  ISEC1_Sec2Or3Flag = 128;
  
  gridtype = gridInqType(gridID);

  ISEC1_GridDefinition = 255;

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
	  Warning("Curvilinear grids are unsupported in GRIB1! Created wrong GDS!");
	}
      gridtype = GRID_LONLAT;
      lcurvi = TRUE;
    }

  ISEC2_Reduced = FALSE;

  ISEC2_ScanFlag = 0;

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

	if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
          ISEC2_GridType = GRIB1_GTYPE_GAUSSIAN;
        else if ( gridtype == GRID_LONLAT && gridIsRotated(gridID) )
	  ISEC2_GridType = GRIB1_GTYPE_LATLON_ROT;
	else
	  ISEC2_GridType = GRIB1_GTYPE_LATLON;

	nlon = gridInqXsize(gridID);
	nlat = gridInqYsize(gridID);

	if ( gridtype == GRID_GAUSSIAN_REDUCED )
	  {
	    ISEC2_Reduced = TRUE;
	    nlon = 0;
	    gridInqRowlon(gridID, ISEC2_RowLonPtr);
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
		if ( lcurvi )
		  xlast  = gridInqXval(gridID, nlon*nlat-1);
		else
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
	    if ( lcurvi )
	      ylast  = gridInqYval(gridID, nlon*nlat-1);
	    else
	      ylast  = gridInqYval(gridID, nlat-1);
	    yinc   = gridInqYinc(gridID);
	    if ( yinc < 0 ) yinc = -yinc;
	  }

	ISEC2_NumLon   = nlon;
	ISEC2_NumLat   = nlat;
	ISEC2_FirstLat = NINT(yfirst*1000);
	ISEC2_LastLat  = NINT(ylast*1000);
	if ( gridtype == GRID_GAUSSIAN_REDUCED )
	  {
	    ISEC2_FirstLon = 0;
	    ISEC2_LastLon  = NINT(1000*(360.-360./(nlat*2)));
	    ISEC2_LonIncr  = NINT(1000*360./(nlat*2));
	  }
	else
	  {
	    ISEC2_FirstLon = NINT(xfirst*1000);
	    ISEC2_LastLon  = NINT(xlast*1000);
	    ISEC2_LonIncr  = NINT(xinc*1000);
	  }
	if ( fabs(xinc*1000 - ISEC2_LonIncr) > FLT_EPSILON )
	  ISEC2_LonIncr = 0;

	if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
          {
            int np = gridInqNP(gridID);
            if ( np == 0 ) np = nlat/2;
            ISEC2_NumPar = np;
          }
	else
	  {
	    ISEC2_LatIncr = NINT(yinc*1000);
	    if ( fabs(yinc*1000 - ISEC2_LatIncr) > FLT_EPSILON )
	      ISEC2_LatIncr = 0;

	    if ( ISEC2_LatIncr < 0 ) ISEC2_LatIncr = -ISEC2_LatIncr;
	  }

	if ( ISEC2_NumLon > 1 && ISEC2_NumLat == 1 )
	  if ( ISEC2_LonIncr != 0 && ISEC2_LatIncr == 0 ) ISEC2_LatIncr = ISEC2_LonIncr;

	if ( ISEC2_NumLon == 1 && ISEC2_NumLat > 1 )
	  if ( ISEC2_LonIncr == 0 && ISEC2_LatIncr != 0 ) ISEC2_LonIncr = ISEC2_LatIncr;

	if ( ISEC2_LatIncr == 0 || ISEC2_LonIncr == 0 )
	  ISEC2_ResFlag = 0;
	else
	  ISEC2_ResFlag = 128;

	if ( gridIsRotated(gridID) )
	  {
	    ISEC2_LatSP = - NINT(gridInqYpole(gridID) * 1000);
	    ISEC2_LonSP =   NINT((gridInqXpole(gridID) + 180) * 1000);
	  }

	/* East -> West */
	if ( ISEC2_LastLon < ISEC2_FirstLon ) ISEC2_ScanFlag += 128;

	/* South -> North */
	if ( ISEC2_LastLat > ISEC2_FirstLat ) ISEC2_ScanFlag += 64;

	break;
      }
    case GRID_LCC:
      {
	double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
	int xsize, ysize;
	int projflag, scanflag;

	xsize = gridInqXsize(gridID);
	ysize = gridInqYsize(gridID);

	gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
		   &projflag, &scanflag);

	ISEC2_GridType = GRIB1_GTYPE_LCC;
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
    case GRID_SPECTRAL:
      {
	ISEC2_GridType = GRIB1_GTYPE_SPECTRAL;
	ISEC2_PentaJ   = gridInqTrunc(gridID);
	ISEC2_PentaK   = ISEC2_PentaJ;
	ISEC2_PentaM   = ISEC2_PentaJ;
	ISEC2_RepType  = 1;
	isec4[2]       = 128;
	if ( gridInqComplexPacking(gridID) && ISEC2_PentaJ >= 21 )
	  {
	    ISEC2_RepMode  = 2;
	    isec4[3]       = 64;
	    isec4[16]      = 0;
	    isec4[17]      = 20;
	    isec4[18]      = 20;
	    isec4[19]      = 20;
	  }
	else
	  {
	    ISEC2_RepMode  = 1;
	    isec4[3]       = 0;
	  }
	break;
      }
    case GRID_GME:
      {
	ISEC2_GridType   = GRIB1_GTYPE_GME;
	ISEC2_GME_ND     = gridInqGMEnd(gridID);
	ISEC2_GME_NI     = gridInqGMEni(gridID);
	ISEC2_GME_NI2    = gridInqGMEni2(gridID);
	ISEC2_GME_NI3    = gridInqGMEni3(gridID);
	ISEC2_GME_AFlag  = 0;
	ISEC2_GME_LatPP  = 90000;
	ISEC2_GME_LonPP  = 0;
	ISEC2_GME_LonMPL = 0;
	ISEC2_GME_BFlag  = 0;
	break;
      }
    default:
      {
	Warning("The CGRIBEX library can not store fields on the used grid!");
	Error("Unsupported grid type: %s", gridNamePtr(gridtype));
      }
    }
}

static
void cgribexDefLevel(int *isec1, int *isec2, double *fsec2, int zaxisID, int levelID)
{
  double level;
  int ilevel, zaxistype, ltype;
  static int warning = 1;
  static int vct_warning = 1;

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

  ISEC2_NumVCP = 0;

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:
      {
	ISEC1_LevelType = GRIB1_LTYPE_SURFACE;
	ISEC1_Level1    = (int) zaxisInqLevel(zaxisID, levelID);
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_TOA:
      {
	ISEC1_LevelType = GRIB1_LTYPE_TOA;
	ISEC1_Level1    = 0;
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_SEA_BOTTOM:
      {
	ISEC1_LevelType = GRIB1_LTYPE_SEA_BOTTOM;
	ISEC1_Level1    = 0;
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_ATMOSPHERE:
      {
	ISEC1_LevelType = GRIB1_LTYPE_ATMOSPHERE;
	ISEC1_Level1    = 0;
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_MEANSEA:
      {
	ISEC1_LevelType = GRIB1_LTYPE_MEANSEA;
	ISEC1_Level1    = (int) zaxisInqLevel(zaxisID, levelID);
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
      {
	int vctsize;

	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	  {
	    ISEC1_LevelType = GRIB1_LTYPE_HYBRID_LAYER;
	    ISEC1_Level1    = (int) zaxisInqLbound(zaxisID, levelID);
	    ISEC1_Level2    = (int) zaxisInqUbound(zaxisID, levelID);
	  }
	else
	  {
	    ISEC1_LevelType = GRIB1_LTYPE_HYBRID;
	    ISEC1_Level1    = (int) zaxisInqLevel(zaxisID, levelID);
	    ISEC1_Level2    = 0;
	  }

	vctsize = zaxisInqVctSize(zaxisID);
	if ( vctsize == 0 && warning )
	  {
	    Warning("VCT missing. ( param = %d, zaxisID = %d )",
		    ISEC1_Parameter, zaxisID);
	    warning = 0;
	  }
	if ( vctsize > 255 )
	  {
	    ISEC2_NumVCP = 0;
	    if ( vct_warning )
	      {
		Warning("VCT size of %d is too large (maximum is 255). Set to 0!", vctsize);
		vct_warning = 0;
	      }
	  }
	else
	  {
	    ISEC2_NumVCP = vctsize;
	    zaxisInqVct(zaxisID, &fsec2[10]);
	  }
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

	ilevel = (int) level;
	if ( level < 32768 && (level < 100 || modf(level/100, &dum) > 0) )
	  {
	    ISEC1_LevelType = GRIB1_LTYPE_99;
	    ISEC1_Level1    = ilevel;
	    ISEC1_Level2    = 0;
	  }
	else
	  {
	    ISEC1_LevelType = GRIB1_LTYPE_ISOBARIC;
	    ISEC1_Level1    = ilevel/100;
	    ISEC1_Level2    = 0;
	  }
	break;
      }
    case ZAXIS_HEIGHT:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	ISEC1_LevelType = GRIB1_LTYPE_HEIGHT;
	ISEC1_Level1    = ilevel;
	ISEC1_Level2    = 0;

	break;
      }
    case ZAXIS_ALTITUDE:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	ISEC1_LevelType = GRIB1_LTYPE_ALTITUDE;
	ISEC1_Level1    = ilevel;
	ISEC1_Level2    = 0;

	break;
      }
    case ZAXIS_SIGMA:
      {
	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	  {
	    ISEC1_LevelType = GRIB1_LTYPE_SIGMA_LAYER;
	    ISEC1_Level1    = (int) zaxisInqLbound(zaxisID, levelID);
	    ISEC1_Level2    = (int) zaxisInqUbound(zaxisID, levelID);
	  }
	else
	  {
            level = zaxisInqLevel(zaxisID, levelID);

            ilevel = (int) level;
            ISEC1_LevelType = GRIB1_LTYPE_SIGMA;
            ISEC1_Level1    = ilevel;
            ISEC1_Level2    = 0;
          }

	break;
      }
    case ZAXIS_DEPTH_BELOW_LAND:
      {
	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	  {
	    ISEC1_LevelType = GRIB1_LTYPE_LANDDEPTH_LAYER;
	    ISEC1_Level1    = (int) zaxisInqLbound(zaxisID, levelID);
	    ISEC1_Level2    = (int) zaxisInqUbound(zaxisID, levelID);
	  }
	else
	  {
	    level = zaxisInqLevel(zaxisID, levelID);

	    ilevel = (int) level;
	    ISEC1_LevelType = GRIB1_LTYPE_LANDDEPTH;
	    ISEC1_Level1    = ilevel;
	    ISEC1_Level2    = 0;
	  }

	break;
      }
    case ZAXIS_DEPTH_BELOW_SEA:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	ISEC1_LevelType = GRIB1_LTYPE_SEADEPTH;
	ISEC1_Level1    = ilevel;
	ISEC1_Level2    = 0;

	break;
      }
    case ZAXIS_ISENTROPIC:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	ISEC1_LevelType = GRIB1_LTYPE_ISENTROPIC;
	ISEC1_Level1    = ilevel;
	ISEC1_Level2    = 0;

	break;
      }
    case ZAXIS_GENERIC:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	ISEC1_LevelType = ltype;
	ISEC1_Level1    = ilevel;
	ISEC1_Level2    = 0;

	break;
      }
    default:
      {
	Error("Unsupported zaxis type: %s", zaxisNamePtr(zaxistype));
	break;
      }
    }
}

static
void cgribexDefMask(int *isec3)
{
}

static
void cgribexDefaultSec0(int *isec0)
{
  ISEC0_GRIB_Len     = 0;
  ISEC0_GRIB_Version = 0;
}

static
void cgribexDefaultSec1(int *isec1)
{
  ISEC1_CenterID    = 0;
  ISEC1_SubCenterID = 0;
  ISEC1_LocalFLag   = 0;
}

static
void cgribexDefaultSec4(int *isec4)
{
  long i;

  for ( i = 2; i <= 10; ++i ) isec4[i] = 0;
}
#endif


size_t cgribexEncode(int varID, int levelID, int vlistID, int gridID, int zaxisID,
		     int vdate, int vtime, int tsteptype, int numavg, 
		     long datasize, const double *data, int nmiss, unsigned char *gribbuffer, size_t gribbuffersize)
{
  size_t nbytes = 0;
#if  defined  (HAVE_LIBCGRIBEX)
  long gribsize;
  int iret = 0, iword = 0;
  int isec0[2], isec1[4096], isec2[4096], isec3[2], isec4[512];
  double fsec2[512], fsec3[2];
  int datatype;
  int param;

  memset(isec1, 0, 32*sizeof(int));
  fsec2[0] = 0; fsec2[1] = 0;

  gribsize = gribbuffersize / sizeof(int);
  param    = vlistInqVarParam(vlistID, varID);

  cgribexDefaultSec0(isec0);
  cgribexDefaultSec1(isec1);
  cgribexDefaultSec4(isec4);

  cgribexDefInstitut(isec1, vlistID, varID);
  cgribexDefModel(isec1, vlistID, varID);

  datatype = vlistInqVarDatatype(vlistID, varID);

  cgribexDefParam(isec1, param);
  cgribexDefTime(isec1, vdate, vtime, tsteptype, numavg, vlistInqTaxis(vlistID));
  cgribexDefGrid(isec1, isec2, isec4, gridID);
  cgribexDefLevel(isec1, isec2, fsec2, zaxisID, levelID);
  cgribexDefMask(isec3);

  ISEC4_NumValues = gridInqSize(gridID);
  ISEC4_NumBits   = grbBitsPerValue(datatype);

  if ( nmiss > 0 )
    {
      FSEC3_MissVal = vlistInqVarMissval(vlistID, varID);
      ISEC1_Sec2Or3Flag |= 64;
    }

  if ( isec4[2] == 128 && isec4[3] == 64 )
    {
      isec4[16] = (int) (1000*calculate_pfactor(data, ISEC2_PentaJ, isec4[17]));
      if ( isec4[16] < -10000 ) isec4[16] = -10000;
      if ( isec4[16] >  10000 ) isec4[16] =  10000;
    }
  //printf("isec4[16] %d\n", isec4[16]);

  gribExDP(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, (double*) data,
	   datasize, (int *) gribbuffer, gribsize, &iword, "C", &iret);

  if ( iret ) Error("Problem during GRIB encode (errno = %d)!", iret);

  nbytes = iword*sizeof(int);
#else
  Error("CGRIBEX support not compiled in!");
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
