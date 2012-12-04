#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif


#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#include "dmemory.h"

#include "cdi.h"
#include "basetime.h"
#include "gaussgrid.h"
#include "stream_int.h"
#include "stream_cdf.h"
#include "cdf_int.h"
#include "varscan.h"
#include "vlist.h"


#if  defined  (HAVE_LIBNETCDF)
#  include "netcdf.h"
#endif

//#define PROJECTION_TEST

#undef  UNDEFID
#define UNDEFID  CDI_UNDEFID


void cdfDefGlobalAtts(int streamID);
void cdfDefLocalAtts(int streamID);

#define  MAXNAMELEN  256

#define  X_AXIS  1
#define  Y_AXIS  2
#define  Z_AXIS  3
#define  T_AXIS  4

typedef struct {
  int     ncvarid;
  int     dimtype;
  size_t  len;
  char    name[MAXNAMELEN];
}
ncdim_t;

typedef struct {
  int      ignore;
  int      isvar;
  int      islon;
  int      islat;
  int      islev;
  int      warn;
  int      timeID;
  int      param;
  int      code;
  int      tabnum;
  int      bounds;
  int      gridID;
  int      zaxisID;
  int      gridtype;
  int      zaxistype;
  int      xdim;
  int      ydim;
  int      zdim;
  int      xvarid;
  int      yvarid;
  int      zvarid;
  int      tvarid;
  int      ncoordvars;
  int      coordvarids[4];
  int      cellarea;
  int      calendar;
  int      tableID;
  int      truncation;
  int      defmiss;
  int      xtype;
  int      ndims;
  int      gmapid;
  int      positive;
  int      dimids[8];
  int      dimtype[8];
  int      natts;
  int     *atts;
  int      deflate;
  int      lunsigned;
  int      lvalidrange;
  size_t   vlen;
  double  *vdata;
  double   missval;
  double   addoffset;
  double   scalefactor;
  double   validrange[2];
  char     name[MAXNAMELEN];
  char     longname[MAXNAMELEN];
  char     stdname[MAXNAMELEN];
  char     units[MAXNAMELEN];
}
ncvar_t;

static
void strtolower(char *str)
{
  int i, len;

  if ( str )
    {
      len = (int) strlen(str);
      for ( i = 0; i < len; i++ )
        str[i] = tolower((int) str[i]);
    }
}

static
int get_timeunit(int len, char *ptu)
{
  int timeunit = -1;

  if ( len > 2 )
    {
      if      ( memcmp(ptu, "sec",    3) == 0 )          timeunit = TUNIT_SECOND;
      else if ( memcmp(ptu, "minute", 6) == 0 )          timeunit = TUNIT_MINUTE;
      else if ( memcmp(ptu, "hour",   4) == 0 )          timeunit = TUNIT_HOUR;
      else if ( memcmp(ptu, "day",    3) == 0 )          timeunit = TUNIT_DAY;
      else if ( memcmp(ptu, "month",  5) == 0 )          timeunit = TUNIT_MONTH;
      else if ( memcmp(ptu, "calendar_month", 14) == 0 ) timeunit = TUNIT_MONTH;
      else if ( memcmp(ptu, "year",   4) == 0 )          timeunit = TUNIT_YEAR;
    }
  else if ( len == 1 )
    {
      if ( ptu[0] == 's' ) timeunit = TUNIT_SECOND;
    }

  return (timeunit);
}

static
int isTimeUnits(const char *timeunits)
{
  int len, i;
  char *ptu, *tu;
  int timetype = -1;
  int timeunit;
  int status = FALSE;

  len = (int) strlen(timeunits);
  tu = (char *) malloc((len+1)*sizeof(char));
  memcpy(tu, timeunits, (len+1)*sizeof(char));
  ptu = tu;

  for ( i = 0; i < len; i++ ) ptu[i] = tolower((int) ptu[i]);

  timeunit = get_timeunit(len, ptu);
  if ( timeunit != -1 )
    {

      while ( ! isspace(*ptu) && *ptu != 0 ) ptu++;
      if ( *ptu )
        {
          while ( isspace(*ptu) ) ptu++;

          if ( memcmp(ptu, "as", 2) == 0 )
            timetype = TAXIS_ABSOLUTE;
          else if ( memcmp(ptu, "since", 5) == 0 )
            timetype = TAXIS_RELATIVE;

          if ( timetype != -1 ) status = TRUE;
        }
    }

  free(tu);

  return (status);
}

static
int splitBasetime(const char *timeunits, taxis_t *taxis)
{
  int len, i;
  char *ptu, *tu;
  int year, month, day;
  int hour = 0, minute = 0, second = 0;
  int timetype = TAXIS_ABSOLUTE;
  int rdate = -1, rtime = -1;
  int timeunit;

  len = (int) strlen(timeunits);
  tu = (char *) malloc((len+1)*sizeof(char));
  memcpy(tu, timeunits, (len+1)*sizeof(char));
  ptu = tu;

  for ( i = 0; i < len; i++ ) ptu[i] = tolower((int) ptu[i]);

  timeunit = get_timeunit(len, ptu);
  if ( timeunit == -1 )
    {
      Message("Unsupported TIMEUNIT: %s!", timeunits);
      return (1);
    }

  while ( ! isspace(*ptu) && *ptu != 0 ) ptu++;
  if ( *ptu )
    {
      while ( isspace(*ptu) ) ptu++;

      if ( memcmp(ptu, "as", 2) == 0 )
        timetype = TAXIS_ABSOLUTE;
      else if ( memcmp(ptu, "since", 5) == 0 )
        timetype = TAXIS_RELATIVE;

      while ( ! isspace(*ptu) && *ptu != 0 ) ptu++;
      if ( *ptu )
        {
          while ( isspace(*ptu) ) ptu++;

          if ( timetype == TAXIS_ABSOLUTE )
            {
              if ( memcmp(ptu, "%y%m%d.%f", 9) != 0 && timeunit == TUNIT_DAY )
                {
                  Message("Unsupported format %s for TIMEUNIT day!", ptu);
                  timeunit = -1;
                }
              else if ( memcmp(ptu, "%y%m.%f", 7) != 0 && timeunit == TUNIT_MONTH )
                {
                  Message("Unsupported format %s for TIMEUNIT month!", ptu);
                  timeunit = -1;
                }
            }
          else if ( timetype == TAXIS_RELATIVE )
            {
              int v1, v2, v3;
              v1 = atoi(ptu);
              if ( v1 < 0 ) ptu++;
              while ( isdigit((int) *ptu) ) ptu++;
              v2 = atoi(++ptu);
              while ( isdigit((int) *ptu) ) ptu++;
              v3 = atoi(++ptu);
              while ( isdigit((int) *ptu) ) ptu++;

              if ( v3 > 999 && v1 < 32 )
                { year = v3; month = v2; day = v1; }
              else
                { year = v1; month = v2; day = v3; }

              while ( isspace((int) *ptu) ) ptu++;

              if ( *ptu )
                {
                  while ( ! isdigit((int) *ptu) ) ptu++;

                  hour = atoi(ptu);
                  while ( isdigit((int) *ptu) ) ptu++;
                  if ( *ptu == ':' )
                    {
                      ptu++;
                      minute = atoi(ptu);
                      while ( isdigit((int) *ptu) ) ptu++;
                      if ( *ptu == ':' )
                        {
                          ptu++;
                          second = atoi(ptu);
                          /*
                          if ( second != 0 )
                            Message("Seconds not supported in time units!");
                          */
                        }
                    }
                }

              rdate = cdiEncodeDate(year, month, day);
              rtime = cdiEncodeTime(hour, minute, second);
              (*taxis).rdate = rdate;
              (*taxis).rtime = rtime;

              if ( CDI_Debug )
                Message("rdate = %d  rtime = %d", rdate, rtime);
            }
        }
    }

  (*taxis).type = timetype;
  (*taxis).unit = timeunit;

  free(tu);

  if ( CDI_Debug )
    Message("timetype = %d  unit = %d", timetype, timeunit);

  return (0);
}


#if  defined  (HAVE_LIBNETCDF)
static
void cdfGetAttInt(int fileID, int ncvarid, char *attname, int attlen, int *attint)
{
  size_t nc_attlen;
  int *pintatt;

  cdf_inq_attlen(fileID, ncvarid, attname, &nc_attlen);

  if ( (int)nc_attlen > attlen )
    pintatt = (int *) malloc(nc_attlen*sizeof(int));
  else
    pintatt = attint;

  cdf_get_att_int(fileID, ncvarid, attname, pintatt);

  if ( (int)nc_attlen > attlen )
    {
      memcpy(attint, pintatt, attlen*sizeof(int));
      free(pintatt);
    }
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static void cdfGetAttDouble(int fileID, int ncvarid, char *attname, int attlen, double *attdouble)
{
  size_t nc_attlen;
  double *pdoubleatt;

  cdf_inq_attlen(fileID, ncvarid, attname, &nc_attlen);

  if ( (int)nc_attlen > attlen )
    pdoubleatt = (double *) malloc(nc_attlen*sizeof(double));
  else
    pdoubleatt = attdouble;

  cdf_get_att_double(fileID, ncvarid, attname, pdoubleatt);

  if ( (int)nc_attlen > attlen )
    {
      memcpy(attdouble, pdoubleatt, attlen*sizeof(double));
      free(pdoubleatt);
    }
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfGetAttText(int fileID, int ncvarid, char *attname, int attlen, char *atttext)
{
  size_t nc_attlen;
  char attbuf[65636];

  cdf_inq_attlen(fileID, ncvarid, attname, &nc_attlen);

  if ( nc_attlen < sizeof(attbuf) )
    {
      cdf_get_att_text(fileID, ncvarid, attname, attbuf);

      attbuf[nc_attlen++] = 0;

      if ( (int) nc_attlen > attlen ) nc_attlen = attlen;
      memcpy(atttext, attbuf, nc_attlen);
    }
  else
    {
      atttext[0] = 0;
    }
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
int cdfInqDatatype(int xtype, int lunsigned)
{
  int datatype = -1;

  if ( xtype == NC_BYTE && lunsigned ) xtype = NC_UBYTE;

  if      ( xtype == NC_BYTE   )  datatype = DATATYPE_INT8;
  else if ( xtype == NC_UBYTE  )  datatype = DATATYPE_UINT8;
  /* else if ( xtype == NC_CHAR   )  datatype = DATATYPE_UINT8; */
  else if ( xtype == NC_SHORT  )  datatype = DATATYPE_INT16;
  else if ( xtype == NC_INT    )  datatype = DATATYPE_INT32;
  else if ( xtype == NC_FLOAT  )  datatype = DATATYPE_FLT32;
  else if ( xtype == NC_DOUBLE )  datatype = DATATYPE_FLT64;
#if  defined  (HAVE_NETCDF4)
  else if ( xtype == NC_LONG   )  datatype = DATATYPE_INT32;
  else if ( xtype == NC_USHORT )  datatype = DATATYPE_UINT16;
  else if ( xtype == NC_UINT   )  datatype = DATATYPE_UINT32;
  else if ( xtype == NC_INT64  )  datatype = DATATYPE_FLT64;
  else if ( xtype == NC_UINT64 )  datatype = DATATYPE_FLT64;
#endif

  return (datatype);
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
int cdfDefDatatype(int datatype, int filetype)
{
  int xtype;

  if ( datatype == DATATYPE_CPX32 || datatype == DATATYPE_CPX64 )
    Error("CDI/netCDF library does not support complex numbers!");

  if ( filetype == FILETYPE_NC4 )
    {
      if      ( datatype == DATATYPE_INT8   ) xtype = NC_BYTE;
      else if ( datatype == DATATYPE_INT16  ) xtype = NC_SHORT;
      else if ( datatype == DATATYPE_INT32  ) xtype = NC_INT;
#if  defined  (HAVE_NETCDF4)
      else if ( datatype == DATATYPE_UINT8  ) xtype = NC_UBYTE;
      else if ( datatype == DATATYPE_UINT16 ) xtype = NC_USHORT;
      else if ( datatype == DATATYPE_UINT32 ) xtype = NC_UINT;
#else
      else if ( datatype == DATATYPE_UINT8  ) xtype = NC_BYTE;
      else if ( datatype == DATATYPE_UINT16 ) xtype = NC_INT;
      else if ( datatype == DATATYPE_UINT32 ) xtype = NC_INT;
#endif
      else if ( datatype == DATATYPE_FLT64  ) xtype = NC_DOUBLE;
      else                                    xtype = NC_FLOAT;
    }
  else
    {
      if      ( datatype == DATATYPE_INT8   ) xtype = NC_BYTE;
      else if ( datatype == DATATYPE_INT16  ) xtype = NC_SHORT;
      else if ( datatype == DATATYPE_INT32  ) xtype = NC_INT;
      else if ( datatype == DATATYPE_UINT8  ) xtype = NC_BYTE;
      else if ( datatype == DATATYPE_UINT16 ) xtype = NC_INT;
      else if ( datatype == DATATYPE_UINT32 ) xtype = NC_INT;
      else if ( datatype == DATATYPE_FLT64  ) xtype = NC_DOUBLE;
      else                                    xtype = NC_FLOAT;
    }


  return (xtype);
}
#endif


#if  defined  (HAVE_LIBNETCDF)
static
void defineAttributes(int vlistID, int varID, int fileID, int ncvarID)
{
  int natts, iatt;
  int atttype, attlen;
  size_t len;
  char attname[1024];

  vlistInqNatts(vlistID, varID, &natts);

  for ( iatt = 0; iatt < natts; iatt++ )
    {
      vlistInqAtt(vlistID, varID, iatt, attname, &atttype, &attlen);

      if ( atttype == DATATYPE_TXT )
        {
          char *atttxt;
          atttxt = (char *) malloc(attlen*sizeof(char));
          vlistInqAttTxt(vlistID, varID, attname, attlen, atttxt);
          len = attlen;
          cdf_put_att_text(fileID, ncvarID, attname, len, atttxt);
          free(atttxt);
        }
      else if ( atttype == DATATYPE_INT16 || atttype == DATATYPE_INT32 )
        {
          int *attint;
          attint = (int *) malloc(attlen*sizeof(int));
          vlistInqAttInt(vlistID, varID, attname, attlen, &attint[0]);
          len = attlen;
          if ( atttype == DATATYPE_INT16 )
            cdf_put_att_int(fileID, ncvarID, attname, NC_SHORT, len, attint);
          else
            cdf_put_att_int(fileID, ncvarID, attname, NC_INT, len, attint);
          free(attint);
        }
      else if ( atttype == DATATYPE_FLT32 || atttype == DATATYPE_FLT64 )
        {
          double *attflt;
          attflt = (double *) malloc(attlen*sizeof(double));
          vlistInqAttFlt(vlistID, varID, attname, attlen, attflt);
          len = attlen;
          if ( atttype == DATATYPE_FLT32 )
            cdf_put_att_double(fileID, ncvarID, attname, NC_FLOAT, len, attflt);
          else
            cdf_put_att_double(fileID, ncvarID, attname, NC_DOUBLE, len, attflt);
          free(attflt);
        }
    }
}
#endif


int cdfCopyRecord(int streamID2, int streamID1)
{
  double *data;
  int datasize;
  int tsID1, tsID2, recID1, recID2;
  int ivarID, gridID;
  int nmiss;
  int ierr = 0;
  int vlistID1, vlistID2;
  stream_t *streamptr1;
  stream_t *streamptr2;

  streamptr1 = stream_to_pointer(streamID1);
  streamptr2 = stream_to_pointer(streamID2);

  stream_check_ptr(__func__, streamptr1);
  stream_check_ptr(__func__, streamptr2);

  vlistID1 = streamptr1->vlistID;
  vlistID2 = streamptr2->vlistID;

  tsID1 = streamptr1->curTsID;
  tsID2 = streamptr2->curTsID;

  recID1 = streamptr1->tsteps[tsID1].curRecID;
  recID2 = streamptr2->tsteps[tsID2].curRecID;

  ivarID = streamptr1->tsteps[tsID1].records[recID1].varID;

  gridID = vlistInqVarGrid(vlistID1, ivarID);

  datasize = gridInqSize(gridID);
  /* bug fix for constant netCDF fields */
  if ( datasize < 1048576 ) datasize = 1048576;

  data = (double *) malloc(datasize*sizeof(double));

  streamReadRecord(streamID1, data, &nmiss);
  streamWriteRecord(streamID2, data, nmiss);

  free(data);

  return (ierr);
}

/* not used
int cdfInqRecord(int streamID, int *varID, int *levelID)
{
  int tsID, recID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  recID = streamptr->tsteps[0].curRecID++;
  printf("cdfInqRecord recID %d %d\n", recID, streamptr->tsteps[0].curRecID);
  printf("cdfInqRecord tsID %d\n", streamptr->curTsID);

  if ( streamptr->tsteps[0].curRecID >= 
       streamptr->tsteps[0].nrecs )
    {
      streamptr->tsteps[0].curRecID = 0;
    }

  *varID   = streamptr->tsteps[0].records[recID].varID;
  *levelID = streamptr->tsteps[0].records[recID].levelID;

  streamptr->record->varID   = *varID;
  streamptr->record->levelID = *levelID;

  if ( CDI_Debug )
    Message("recID = %d  varID = %d  levelID = %d", recID, *varID, *levelID);
  
  return (recID+1);
}
*/
int cdfDefRecord(int streamID)
{
  int ierr = 0;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( CDI_Debug )
    Message("streamID = %d", streamID);

  stream_check_ptr(__func__, streamptr);

  return (ierr);
}

#if  defined  (HAVE_LIBNETCDF)
static
void cdfWriteGridTraj(int streamID, int gridID)
{
  int tsID, fileID;
  int lonID, latID, gridindex;
  size_t index;
  double xlon, xlat;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  gridindex = vlistGridIndex(vlistID, gridID);
  lonID = streamptr->xdimID[gridindex];
  latID = streamptr->ydimID[gridindex];

  xlon = gridInqXval(gridID, 0);
  xlat = gridInqYval(gridID, 0);
  tsID = streamptr->curTsID;
  index = tsID;

  cdf_put_var1_double(fileID, lonID, &index, &xlon);
  cdf_put_var1_double(fileID, latID, &index, &xlat);
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfReadGridTraj(int streamID, int gridID)
{
  int tsID, fileID;
  int lonID, latID, gridindex;
  size_t index;
  double xlon, xlat;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  gridindex = vlistGridIndex(vlistID, gridID);
  lonID = streamptr->xdimID[gridindex];
  latID = streamptr->ydimID[gridindex];

  tsID = streamptr->curTsID;
  index = tsID;

  cdf_get_var1_double(fileID, lonID, &index, &xlon);
  cdf_get_var1_double(fileID, latID, &index, &xlat);

  gridDefXvals(gridID, &xlon);
  gridDefYvals(gridID, &xlat);
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefVarDeflate(int ncid, int ncvarid, int deflate_level)
{
#if  defined  (HAVE_NETCDF4)
  int retval;
  /* Set chunking, shuffle, and deflate. */
  int shuffle = 1;
  int deflate = 1;

  if ( deflate_level < 1 || deflate_level > 9 ) deflate_level = 1;

  if ((retval = nc_def_var_deflate(ncid, ncvarid, shuffle, deflate, deflate_level)))
    {
      Error("nc_def_var_deflate failed, status = %d", retval);
    }
#else
  static int lwarn = TRUE;

  if ( lwarn )
    {
      lwarn = FALSE;
      Warning("Deflate compression failed, netCDF4 not available!");
    }
#endif
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefVarSzip(int ncid, int ncvarid)
{
#if  defined  (HAVE_NETCDF4) && defined (NC_SZIP_NN_OPTION_MASK)
  int retval;
  /* Set options_mask and bits_per_pixel. */
  int options_mask = NC_SZIP_NN_OPTION_MASK;
  int bits_per_pixel = 16;

  if ((retval = nc_def_var_szip(ncid, ncvarid, options_mask, bits_per_pixel)))
    {
      if ( retval == NC_EINVAL )
        {
          static int lwarn = TRUE;

          if ( lwarn )
            {
              lwarn = FALSE;
              Warning("netCDF4/Szip compression not compiled in!");
            }
        }
      else
        Error("nc_def_var_szip failed, status = %d", retval);
    }
#else
  static int lwarn = TRUE;

  if ( lwarn )
    {
      lwarn = FALSE;
      Warning("netCDF4/Szip compression not available!");
    }
#endif
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefVarMissval(int streamID, int varID, int dtype, int lcheck)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( streamptr->vars[varID].defmiss == FALSE )
    {
      int fileID;
      int ncvarid;
      double missval;
      int vlistID;
      int xtype;

      vlistID = streamInqVlist(streamID);
      fileID  = streamInqFileID(streamID);
      ncvarid = streamptr->vars[varID].ncvarid;
      missval = vlistInqVarMissval(vlistID, varID);
      if ( lcheck && streamptr->ncmode == 2 ) cdf_redef(fileID);

      xtype = cdfDefDatatype(dtype, streamptr->filetype);

      cdf_put_att_double(fileID, ncvarid, "_FillValue", (nc_type) xtype, 1, &missval);

      if ( cdiNcMissingValue == 1 )
        cdf_put_att_double(fileID, ncvarid, "missing_value", (nc_type) xtype, 1, &missval);

      if ( lcheck && streamptr->ncmode == 2 ) cdf_enddef(fileID);

      streamptr->vars[varID].defmiss = TRUE;
    }
}
#endif

void cdfWriteRecord(int streamID, const double *data, int nmiss)
{
#if  defined  (HAVE_LIBNETCDF)
  int varID;
  int levelID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  varID   = streamptr->record->varID;
  levelID = streamptr->record->levelID;

  if ( CDI_Debug )
    Message("streamID = %d  varID = %d", streamID, varID);

  cdfWriteVarSliceDP(streamID, varID, levelID, data, nmiss);
#endif
}


int cdfReadRecord(int streamID, double *data, int *nmiss)
{
  int ierr = 0;
  int levelID, varID, tsID, recID, vrecID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( CDI_Debug ) Message("streamID = %d", streamID);

  tsID    = streamptr->curTsID;
  vrecID  = streamptr->tsteps[tsID].curRecID;
  recID   = streamptr->tsteps[tsID].recIDs[vrecID];
  varID   = streamptr->tsteps[tsID].records[recID].varID;
  levelID = streamptr->tsteps[tsID].records[recID].levelID;

  cdfReadVarSliceDP(streamID, varID, levelID, data, nmiss);

  return (ierr);
}

static
void cdfDefTimeValue(int streamID, int tsID)
{
#if  defined  (HAVE_LIBNETCDF)
  int fileID;
  double timevalue;
  int ncvarid;
  size_t index;
  taxis_t *taxis;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  fileID = streamInqFileID(streamID);

  if ( CDI_Debug )
    Message("streamID = %d, fileID = %d", streamID, fileID);

  taxis = &streamptr->tsteps[tsID].taxis;

  if ( streamptr->ncmode == 1 )
    {
      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  index = tsID;

  timevalue = cdiEncodeTimeval(taxis->vdate, taxis->vtime, &streamptr->tsteps[0].taxis);
  if ( CDI_Debug ) Message("tsID = %d  timevalue = %f", tsID, timevalue);

  ncvarid = streamptr->basetime.ncvarid;
  cdf_put_var1_double(fileID, ncvarid, &index, &timevalue);

  if ( taxis->has_bounds )
    {
      size_t start[2], count[2];

      ncvarid = streamptr->basetime.ncvarboundsid;

      timevalue = cdiEncodeTimeval(taxis->vdate_lb, taxis->vtime_lb, &streamptr->tsteps[0].taxis);
      start[0] = tsID; count[0] = 1; start[1] = 0; count[1] = 1;
      cdf_put_vara_double(fileID, ncvarid, start, count, &timevalue);

      timevalue = cdiEncodeTimeval(taxis->vdate_ub, taxis->vtime_ub, &streamptr->tsteps[0].taxis);
      start[0] = tsID; count[0] = 1; start[1] = 1; count[1] = 1;
      cdf_put_vara_double(fileID, ncvarid, start, count, &timevalue);
    }
  /*
printf("fileID = %d %d %d %f\n", fileID, time_varid, index, timevalue);
  */
#endif
}

static
void cdfDefTime(int streamID)
{
#if  defined  (HAVE_LIBNETCDF)
  int fileID;
  int time_varid;
  int time_bndsid;
  int dims[2];
  int year, month, day, hour, minute, second;
  char unitstr[80];
  char calstr[80];
  char tmpstr[256];
  size_t len;
  taxis_t *taxis;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( streamptr->basetime.ncvarid != UNDEFID ) return;

  fileID = streamInqFileID(streamID);

  if ( streamptr->ncmode == 0 ) streamptr->ncmode = 1;

  if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

  cdf_def_dim(fileID, "time", NC_UNLIMITED, &streamptr->basetime.ncdimid);

  dims[0] = streamptr->basetime.ncdimid;
  cdf_def_var(fileID, "time", NC_DOUBLE, 1, dims, &time_varid);

  streamptr->basetime.ncvarid = time_varid;

  taxis = &streamptr->tsteps[0].taxis;


  strcpy(tmpstr, "time");
  cdf_put_att_text(fileID, time_varid, "standard_name", strlen(tmpstr), tmpstr);

  if ( taxis->has_bounds )
    {
      /* fprintf(stderr, "time has bounds\n"); */

      if ( nc_inq_dimid(fileID, "nb2", &dims[1]) != NC_NOERR )
	cdf_def_dim(fileID, "nb2", 2, &dims[1]);

      cdf_def_var(fileID, "time_bnds", NC_DOUBLE, 2, dims, &time_bndsid);

      streamptr->basetime.ncvarboundsid = time_bndsid;

      cdf_put_att_text(fileID, time_varid, "bounds", 9, "time_bnds");
    }

  unitstr[0] = 0;
  if ( streamptr->tsteps[0].taxis.type == TAXIS_ABSOLUTE )
    {
      if ( streamptr->tsteps[0].taxis.unit == TUNIT_YEAR )
        sprintf(unitstr, "year as %s", "%Y.%f");
      else if ( streamptr->tsteps[0].taxis.unit == TUNIT_MONTH )
        sprintf(unitstr, "month as %s", "%Y%m.%f");
      else
        sprintf(unitstr, "day as %s", "%Y%m%d.%f");
    }
  else
    {
      int rdate, rtime;
      int timeunit;

      timeunit = taxis->unit;
      if ( timeunit == -1 ) timeunit = TUNIT_HOUR;
      rdate    = taxis->rdate;
      rtime    = taxis->rtime;
      if ( rdate == -1 )
        {
          rdate  = taxis->vdate;
          rtime  = taxis->vtime;
        }

      cdiDecodeDate(rdate, &year, &month, &day);
      cdiDecodeTime(rtime, &hour, &minute, &second);

      if ( timeunit == TUNIT_QUARTER ) timeunit = TUNIT_MINUTE;
      if ( timeunit == TUNIT_3HOURS  ||
	   timeunit == TUNIT_6HOURS  ||
	   timeunit == TUNIT_12HOURS ) timeunit = TUNIT_HOUR;

      sprintf(unitstr, "%s since %d-%02d-%02d %02d:%02d:%02d",
              tunitNamePtr(timeunit), year, month, day, hour, minute, second);
    }

  len = strlen(unitstr);
  if ( len )
    cdf_put_att_text(fileID, time_varid, "units", len, unitstr);

  if ( taxis->has_bounds )
    if ( len )
      cdf_put_att_text(fileID, time_bndsid, "units", len, unitstr);

  if ( taxis->calendar != -1 )
    {
      calstr[0] = 0;

      if      ( taxis->calendar == CALENDAR_STANDARD )  strcpy(calstr, "standard");
      else if ( taxis->calendar == CALENDAR_PROLEPTIC ) strcpy(calstr, "proleptic_gregorian");
      else if ( taxis->calendar == CALENDAR_NONE )      strcpy(calstr, "none");
      else if ( taxis->calendar == CALENDAR_360DAYS )   strcpy(calstr, "360_day");
      else if ( taxis->calendar == CALENDAR_365DAYS )   strcpy(calstr, "365_day");
      else if ( taxis->calendar == CALENDAR_366DAYS )   strcpy(calstr, "366_day");

      len = strlen(calstr);
      if ( len )
        {
          cdf_put_att_text(fileID, time_varid, "calendar", len, calstr);

          if ( taxis->has_bounds )
            cdf_put_att_text(fileID, time_bndsid, "calendar", len, calstr);
        }
    }

  if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
#endif
}


void cdfDefTimestep(int streamID, int tsID)
{
  int vlistID;

  vlistID = streamInqVlist(streamID);

  if ( vlistHasTime(vlistID) ) cdfDefTime(streamID);

  cdfDefTimeValue(streamID, tsID);
}

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefComplex(int streamID, int gridID)
{
  char axisname[] = "nc2";
  int index;
  int dimID = UNDEFID;
  int gridID0, gridtype0, gridindex;
  int ngrids;
  int fileID;
  int dimlen;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  ngrids = vlistNgrids(vlistID);

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->xdimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_SPECTRAL || gridtype0 == GRID_FOURIER )
            {
              dimID = streamptr->xdimID[index];
              break;
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      dimlen = 2;

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  gridindex = vlistGridIndex(vlistID, gridID);
  streamptr->xdimID[gridindex] = dimID;
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefSP(int streamID, int gridID)
{
  /*
  char longname[] = "Spherical harmonic coefficient";
  */
  char axisname[5] = "nspX";
  int index, iz = 0;
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int ngrids;
  int fileID;
  int dimlen, dimlen0;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  ngrids = vlistNgrids(vlistID);

  dimlen = gridInqSize(gridID)/2;

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->ydimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_SPECTRAL )
            {
              dimlen0 = gridInqSize(gridID0)/2;
              if ( dimlen == dimlen0 )
                {
                  dimID = streamptr->ydimID[index];
                  break;
                }
              else
                iz++;
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      if ( iz == 0 ) axisname[3] = '\0';
      else           sprintf(&axisname[3], "%1d", iz+1);

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  gridindex = vlistGridIndex(vlistID, gridID);
  streamptr->ydimID[gridindex] = dimID;
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefFC(int streamID, int gridID)
{
  char axisname[5] = "nfcX";
  int index, iz = 0;
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int ngrids;
  int fileID;
  int dimlen, dimlen0;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  ngrids = vlistNgrids(vlistID);

  dimlen = gridInqSize(gridID)/2;

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->ydimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_FOURIER )
            {
              dimlen0 = gridInqSize(gridID0)/2;
              if ( dimlen == dimlen0 )
                {
                  dimID = streamptr->ydimID[index];
                  break;
                }
              else
                iz++;
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      if ( iz == 0 ) axisname[3] = '\0';
      else           sprintf(&axisname[3], "%1d", iz+1);

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  gridindex = vlistGridIndex(vlistID, gridID);
  streamptr->ydimID[gridindex] = dimID;
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefTrajLon(int streamID, int gridID)
{
  char units[256];
  char longname[256];
  char stdname[256];
  char axisname[256];
  int gridtype, gridindex;
  int dimID = UNDEFID;
  int fileID;
  int dimlen;
  size_t len;
  int ncvarid;
  int vlistID;
  int xtype = NC_DOUBLE;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( gridInqPrec(gridID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  gridtype = gridInqType(gridID);
  dimlen = gridInqXsize(gridID);
  if ( dimlen != 1 ) Error("Xsize isn't 1 for %s grid!", gridNamePtr(gridtype));

  gridindex = vlistGridIndex(vlistID, gridID);
  ncvarid = streamptr->xdimID[gridindex];

  gridInqXname(gridID, axisname);
  gridInqXlongname(gridID, longname);
  gridInqXstdname(gridID, stdname);
  gridInqXunits(gridID, units);

  if ( ncvarid == UNDEFID )
    {
      dimID = streamptr->basetime.ncvarid;

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_var(fileID, axisname, (nc_type) xtype, 1, &dimID, &ncvarid);

      if ( (len = strlen(stdname)) )
        cdf_put_att_text(fileID, ncvarid, "standard_name", len, stdname);
      if ( (len = strlen(longname)) )
        cdf_put_att_text(fileID, ncvarid, "long_name", len, longname);
      if ( (len = strlen(units)) )
        cdf_put_att_text(fileID, ncvarid, "units", len, units);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  streamptr->xdimID[gridindex] = ncvarid; /* var ID for trajectory !!! */
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefTrajLat(int streamID, int gridID)
{
  char units[] = "degrees_north";
  char longname[] = "latitude";
  char stdname[] = "latitude";
  char axisname[] = "tlat";
  int gridtype, gridindex;
  int dimID = UNDEFID;
  int fileID;
  int dimlen;
  size_t len;
  int ncvarid;
  int vlistID;
  int xtype = NC_DOUBLE;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( gridInqPrec(gridID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamInqVlist(streamID);
  fileID = streamInqFileID(streamID);

  gridtype = gridInqType(gridID);
  dimlen = gridInqYsize(gridID);
  if ( dimlen != 1 ) Error("Ysize isn't 1 for %s grid!", gridNamePtr(gridtype));

  gridindex = vlistGridIndex(vlistID, gridID);
  ncvarid = streamptr->ydimID[gridindex];

  gridInqYname(gridID, axisname);
  gridInqYlongname(gridID, longname);
  gridInqYstdname(gridID, stdname);
  gridInqYunits(gridID, units);

  if ( ncvarid == UNDEFID )
    {
      dimID = streamptr->basetime.ncvarid;

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_var(fileID, axisname, (nc_type) xtype, 1, &dimID, &ncvarid);

      if ( (len = strlen(stdname)) )
        cdf_put_att_text(fileID, ncvarid, "standard_name", len, stdname);
      if ( (len = strlen(longname)) )
        cdf_put_att_text(fileID, ncvarid, "long_name", len, longname);
      if ( (len = strlen(units)) )
        cdf_put_att_text(fileID, ncvarid, "units", len, units);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  streamptr->ydimID[gridindex] = ncvarid; /* var ID for trajectory !!! */
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
int checkGridName(int type, char *axisname, int fileID, int vlistID, int gridID, int ngrids, int mode)
{
  int iz, index;
  int gridID0;
  int ncdimid;
  char axisname0[256];
  char axisname2[256];
  int checkname;
  int status;

  /* check that the name is not already defined */
  checkname = TRUE;
  iz = 0;

  while ( checkname ) 
    {
      strcpy(axisname2, axisname);
      if ( iz ) sprintf(&axisname2[strlen(axisname2)], "_%d", iz+1);

      //status = nc_inq_varid(fileID, axisname2, &ncvarid);
      if ( type == 'V' ) /* type Var oder Dim */
        status = nc_inq_varid(fileID, axisname2, &ncdimid);
      else
        status = nc_inq_dimid(fileID, axisname2, &ncdimid);

      if ( status != NC_NOERR )
        {
          if ( iz )
            {
              /* check that the name does not exist for other grids */
              for ( index = 0; index < ngrids; index++ )
                {
                  gridID0 = vlistGrid(vlistID, index);
                  if ( gridID != gridID0 )
                    {
                      if ( mode == 'X' ) /* mode X or Y */
                        gridInqXname(gridID0, axisname0);
                      else
                        gridInqYname(gridID0, axisname0);

                      if ( strcmp(axisname0, axisname2) == 0 ) break;
                    }
                }
              if ( index == ngrids ) checkname = FALSE;
            }
          else
            {
              checkname = FALSE;
            }
        }

      if ( checkname ) iz++;

      if ( iz > 99 ) break;
    }

  if ( iz ) sprintf(&axisname[strlen(axisname)], "_%d", iz+1);

  return (iz);
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefXaxis(int streamID, int gridID)
{
  char units[256];
  char longname[256];
  char stdname[256];
  char axisname[256];
  int index;
  /*  int index2; */
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int dimIDs[2];
  int ngrids;
  int fileID;
  int dimlen, dimlen0;
  size_t len;
  int ncvarid = UNDEFID, ncbvarid = UNDEFID;
  int nvertex = 2, nvdimID = UNDEFID;
  int vlistID;
  int xtype = NC_DOUBLE;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( gridInqPrec(gridID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  ngrids = vlistNgrids(vlistID);

  dimlen = gridInqXsize(gridID);
  gridindex = vlistGridIndex(vlistID, gridID);

  gridInqXname(gridID, axisname);
  gridInqXlongname(gridID, longname);
  gridInqXstdname(gridID, stdname);
  gridInqXunits(gridID, units);

  if ( axisname[0] == 0 ) Error("axis name undefined!");

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->xdimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_GAUSSIAN    ||
               gridtype0 == GRID_LONLAT      ||
               gridtype0 == GRID_CURVILINEAR ||
               gridtype0 == GRID_GENERIC )
            {
              dimlen0 = gridInqXsize(gridID0);
              if ( dimlen == dimlen0 )
                if ( IS_EQUAL(gridInqXval(gridID0, 0), gridInqXval(gridID, 0)) &&
                     IS_EQUAL(gridInqXval(gridID0, dimlen-1), gridInqXval(gridID, dimlen-1)) )
                  {
                    dimID = streamptr->xdimID[index];
                    break;
                  }
              /*
              for ( index2 = 0; index2 < index; index2++ )
                if ( streamptr->xdimID[index] == streamptr->xdimID[index2] )
                  break;
              if ( index2 == index ) iz++;
              */
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      int status;
      status = checkGridName('V', axisname, fileID, vlistID, gridID, ngrids, 'X');
      if ( status == 0 )
        status = checkGridName('D', axisname, fileID, vlistID, gridID, ngrids, 'X');

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      if ( gridInqXboundsPtr(gridID) || gridInqYboundsPtr(gridID) )
        {
          if ( nc_inq_dimid(fileID, "nb2", &nvdimID) != NC_NOERR )
            cdf_def_dim(fileID, "nb2", nvertex, &nvdimID);
        }

      if ( gridInqXvalsPtr(gridID) )
        {
          cdf_def_var(fileID, axisname, (nc_type) xtype, 1, &dimID, &ncvarid);

          if ( (len = strlen(stdname)) )
            cdf_put_att_text(fileID, ncvarid, "standard_name", len, stdname);
          if ( (len = strlen(longname)) )
            cdf_put_att_text(fileID, ncvarid, "long_name", len, longname);
          if ( (len = strlen(units)) )
            cdf_put_att_text(fileID, ncvarid, "units", len, units);

          cdf_put_att_text(fileID, ncvarid, "axis", 1, "X");

          if ( gridInqXboundsPtr(gridID) && nvdimID != UNDEFID )
            {
              strcat(axisname, "_bnds");
              dimIDs[0] = dimID;
              dimIDs[1] = nvdimID;
              cdf_def_var(fileID, axisname, (nc_type) xtype, 2, dimIDs, &ncbvarid);
              cdf_put_att_text(fileID, ncvarid, "bounds", strlen(axisname), axisname);
            }
          /*
          if ( gridIsRotated(gridID) )
            {
              double north_pole = gridInqXpole(gridID);
              cdf_put_att_double(fileID, ncvarid, "north_pole", NC_DOUBLE, 1, &north_pole);
            }
          */
        }

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      if ( ncvarid  != UNDEFID ) cdf_put_var_double(fileID, ncvarid, gridInqXvalsPtr(gridID));
      if ( ncbvarid != UNDEFID ) cdf_put_var_double(fileID, ncbvarid, gridInqXboundsPtr(gridID));
    }

  streamptr->xdimID[gridindex] = dimID;
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefYaxis(int streamID, int gridID)
{
  char units[256];
  char longname[256];
  char stdname[256];
  char axisname[256];
  int index;
  /*  int index2; */
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int dimIDs[2];
  int ngrids;
  int fileID;
  int dimlen, dimlen0;
  size_t len;
  int ncvarid = UNDEFID, ncbvarid = UNDEFID;
  int nvertex = 2, nvdimID = UNDEFID;
  int vlistID;
  int xtype = NC_DOUBLE;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( gridInqPrec(gridID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  ngrids = vlistNgrids(vlistID);

  dimlen = gridInqYsize(gridID);
  gridindex = vlistGridIndex(vlistID, gridID);

  gridInqYname(gridID, axisname);
  gridInqYlongname(gridID, longname);
  gridInqYstdname(gridID, stdname);
  gridInqYunits(gridID, units);

  if ( axisname[0] == 0 ) Error("axis name undefined!");

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->ydimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_GAUSSIAN    ||
               gridtype0 == GRID_LONLAT      ||
               gridtype0 == GRID_CURVILINEAR ||
               gridtype0 == GRID_GENERIC )
            {
              dimlen0 = gridInqYsize(gridID0);
              if ( dimlen == dimlen0 )
                if ( IS_EQUAL(gridInqYval(gridID0, 0), gridInqYval(gridID, 0)) &&
                     IS_EQUAL(gridInqYval(gridID0, dimlen-1), gridInqYval(gridID, dimlen-1)) )
                  {
                    dimID = streamptr->ydimID[index];
                    break;
                  }
              /*
              for ( index2 = 0; index2 < index; index2++ )
                if ( streamptr->ydimID[index] == streamptr->ydimID[index2] )
                  break;
              if ( index2 == index ) iz++;
              */
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      int status;
      status = checkGridName('V', axisname, fileID, vlistID, gridID, ngrids, 'Y');
      if ( status == 0 )
        status = checkGridName('D', axisname, fileID, vlistID, gridID, ngrids, 'Y');

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      if ( gridInqXboundsPtr(gridID) || gridInqYboundsPtr(gridID) )
        {
          if ( nc_inq_dimid(fileID, "nb2", &nvdimID) != NC_NOERR )
            cdf_def_dim(fileID, "nb2", nvertex, &nvdimID);
        }

      if ( gridInqYvalsPtr(gridID) )
        {
          cdf_def_var(fileID, axisname, (nc_type) xtype, 1, &dimID, &ncvarid);

          if ( (len = strlen(stdname)) )
            cdf_put_att_text(fileID, ncvarid, "standard_name", len, stdname);
          if ( (len = strlen(longname)) )
            cdf_put_att_text(fileID, ncvarid, "long_name", len, longname);
          if ( (len = strlen(units)) )
            cdf_put_att_text(fileID, ncvarid, "units", len, units);

          cdf_put_att_text(fileID, ncvarid, "axis", 1, "Y");

          if ( gridInqYboundsPtr(gridID) && nvdimID != UNDEFID )
            {
              strcat(axisname, "_bnds");
              dimIDs[0] = dimID;
              dimIDs[1] = nvdimID;
              cdf_def_var(fileID, axisname, (nc_type) xtype, 2, dimIDs, &ncbvarid);
              cdf_put_att_text(fileID, ncvarid, "bounds", strlen(axisname), axisname);
            }
          /*
          if ( gridIsRotated(gridID) )
            {
              double north_pole = gridInqYpole(gridID);
              cdf_put_att_double(fileID, ncvarid, "north_pole", NC_DOUBLE, 1, &north_pole);
            }
          */
        }

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      if ( ncvarid  != UNDEFID ) cdf_put_var_double(fileID, ncvarid, gridInqYvalsPtr(gridID));
      if ( ncbvarid != UNDEFID ) cdf_put_var_double(fileID, ncbvarid, gridInqYboundsPtr(gridID));
    }

  streamptr->ydimID[gridindex] = dimID;
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefCurvilinear(int streamID, int gridID)
{
  char xunits[256];
  char xlongname[256];
  char xstdname[256];
  char yunits[256];
  char ylongname[256];
  char ystdname[256];
  char xaxisname[256];
  char yaxisname[256];
  char xdimname[4] = "x";
  char ydimname[4] = "y";
  int index;
  int gridID0, gridtype0, gridindex;
  int xdimID = UNDEFID;
  int ydimID = UNDEFID;
  int dimIDs[3];
  int ngrids;
  int fileID;
  int xdimlen, ydimlen, dimlen0;
  size_t len;
  int ncxvarid = UNDEFID, ncyvarid = UNDEFID;
  int ncbxvarid = UNDEFID, ncbyvarid = UNDEFID, ncavarid = UNDEFID;
  int nvertex = 4, nvdimID = UNDEFID;
  int vlistID;
  int xtype = NC_DOUBLE;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( gridInqPrec(gridID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  ngrids = vlistNgrids(vlistID);

  xdimlen = gridInqXsize(gridID);
  ydimlen = gridInqYsize(gridID);
  gridindex = vlistGridIndex(vlistID, gridID);

  gridInqXname(gridID, xaxisname);
  gridInqXlongname(gridID, xlongname);
  gridInqXstdname(gridID, xstdname);
  gridInqXunits(gridID, xunits);
  gridInqYname(gridID, yaxisname);
  gridInqYlongname(gridID, ylongname);
  gridInqYstdname(gridID, ystdname);
  gridInqYunits(gridID, yunits);

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->xdimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_GAUSSIAN    ||
               gridtype0 == GRID_LONLAT      ||
               gridtype0 == GRID_CURVILINEAR ||
               gridtype0 == GRID_GENERIC )
            {
              dimlen0 = gridInqXsize(gridID0);
              if ( xdimlen == dimlen0 )
                if ( IS_EQUAL(gridInqXval(gridID0, 0), gridInqXval(gridID, 0)) &&
                     IS_EQUAL(gridInqXval(gridID0, xdimlen-1), gridInqXval(gridID, xdimlen-1)) )
                  {
                    xdimID = streamptr->xdimID[index];
                    break;
                  }
            }
        }
    }

  if ( xdimID == UNDEFID )
    {
      int status;
      status = checkGridName('V', xaxisname, fileID, vlistID, gridID, ngrids, 'X');
      status = checkGridName('V', yaxisname, fileID, vlistID, gridID, ngrids, 'Y');
      status = checkGridName('D', xdimname, fileID, vlistID, gridID, ngrids, 'X');
      status = checkGridName('D', ydimname, fileID, vlistID, gridID, ngrids, 'Y');

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, xdimname, xdimlen, &xdimID);
      cdf_def_dim(fileID, ydimname, ydimlen, &ydimID);

      if ( gridInqXboundsPtr(gridID) || gridInqYboundsPtr(gridID) )
        {
          if ( nc_inq_dimid(fileID, "nv4", &nvdimID) != NC_NOERR )
            cdf_def_dim(fileID, "nv4", nvertex, &nvdimID);
        }

      dimIDs[0] = ydimID;
      dimIDs[1] = xdimID;

      if ( gridInqXvalsPtr(gridID) )
        {
          cdf_def_var(fileID, xaxisname, (nc_type) xtype, 2, dimIDs, &ncxvarid);

          if ( (len = strlen(xstdname)) )
            cdf_put_att_text(fileID, ncxvarid, "standard_name", len, xstdname);
          if ( (len = strlen(xlongname)) )
            cdf_put_att_text(fileID, ncxvarid, "long_name", len, xlongname);
          if ( (len = strlen(xunits)) )
            cdf_put_att_text(fileID, ncxvarid, "units", len, xunits);

          /* attribute for Panoply */
          cdf_put_att_text(fileID, ncxvarid, "_CoordinateAxisType", 3, "Lon");

          streamptr->ncxvarID[gridindex] = ncxvarid;

          if ( gridInqXboundsPtr(gridID) && nvdimID != UNDEFID )
            {
              strcat(xaxisname, "_bnds");
              dimIDs[0] = ydimID;
              dimIDs[1] = xdimID;
              dimIDs[2] = nvdimID;
              cdf_def_var(fileID, xaxisname, (nc_type) xtype, 3, dimIDs, &ncbxvarid);
              cdf_put_att_text(fileID, ncxvarid, "bounds", strlen(xaxisname), xaxisname);
            }
        }

      if ( gridInqYvalsPtr(gridID) )
        {
          cdf_def_var(fileID, yaxisname, (nc_type) xtype, 2, dimIDs, &ncyvarid);

          if ( (len = strlen(ystdname)) )
            cdf_put_att_text(fileID, ncyvarid, "standard_name", len, ystdname);
          if ( (len = strlen(ylongname)) )
            cdf_put_att_text(fileID, ncyvarid, "long_name", len, ylongname);
          if ( (len = strlen(yunits)) )
            cdf_put_att_text(fileID, ncyvarid, "units", len, yunits);

          /* attribute for Panoply */
          cdf_put_att_text(fileID, ncyvarid, "_CoordinateAxisType", 3, "Lat");

          streamptr->ncyvarID[gridindex] = ncyvarid;

          if ( gridInqYboundsPtr(gridID) && nvdimID != UNDEFID )
            {
              strcat(yaxisname, "_bnds");
              dimIDs[0] = ydimID;
              dimIDs[1] = xdimID;
              dimIDs[2] = nvdimID;
              cdf_def_var(fileID, yaxisname, (nc_type) xtype, 3, dimIDs, &ncbyvarid);
              cdf_put_att_text(fileID, ncyvarid, "bounds", strlen(yaxisname), yaxisname);
            }
        }

      if ( gridInqAreaPtr(gridID) )
        {
          char yaxisname[] = "cell_area";
          char units[] = "m2";
          char longname[] = "area of grid cell";
          char stdname[] = "cell_area";

          cdf_def_var(fileID, yaxisname, (nc_type) xtype, 2, dimIDs, &ncavarid);

          cdf_put_att_text(fileID, ncavarid, "standard_name", strlen(stdname), stdname);
          cdf_put_att_text(fileID, ncavarid, "long_name", strlen(longname), longname);
          cdf_put_att_text(fileID, ncavarid, "units", strlen(units), units);

          streamptr->ncavarID[gridindex] = ncavarid;
        }

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      if ( ncxvarid  != UNDEFID ) cdf_put_var_double(fileID, ncxvarid,  gridInqXvalsPtr(gridID));
      if ( ncbxvarid != UNDEFID ) cdf_put_var_double(fileID, ncbxvarid, gridInqXboundsPtr(gridID));
      if ( ncyvarid  != UNDEFID ) cdf_put_var_double(fileID, ncyvarid,  gridInqYvalsPtr(gridID));
      if ( ncbyvarid != UNDEFID ) cdf_put_var_double(fileID, ncbyvarid, gridInqYboundsPtr(gridID));
      if ( ncavarid  != UNDEFID ) cdf_put_var_double(fileID, ncavarid,  gridInqAreaPtr(gridID));
    }

  streamptr->xdimID[gridindex] = xdimID;
  streamptr->ydimID[gridindex] = ydimID;
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefRgrid(int streamID, int gridID)
{
  char axisname[7] = "rgridX";
  int index, iz = 0;
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int ngrids;
  int fileID;
  int dimlen, dimlen0;
  int vlistID;
  int lwarn = TRUE;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  ngrids = vlistNgrids(vlistID);

  dimlen = gridInqSize(gridID);

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->xdimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_GAUSSIAN_REDUCED )
            {
              dimlen0 = gridInqSize(gridID0);

              if ( dimlen == dimlen0 )
                {
                  dimID = streamptr->xdimID[index];
                  break;
                }
              else
                iz++;   
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      if ( lwarn )
        {
          Warning("Creating a netCDF file with data on a gaussian reduced grid.");
          Warning("The further processing of the resulting file is unsupported!");
          lwarn = FALSE;
        }

      if ( iz == 0 ) axisname[5] = '\0';
      else           sprintf(&axisname[5], "%1d", iz+1);

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  gridindex = vlistGridIndex(vlistID, gridID);
  streamptr->xdimID[gridindex] = dimID;
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefGdim(int streamID, int gridID)
{
  char axisname[7] = "gsizeX";
  int index, iz = 0;
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int ngrids;
  int fileID;
  int dimlen, dimlen0;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  ngrids = vlistNgrids(vlistID);

  dimlen = gridInqSize(gridID);

  if ( gridInqYsize(gridID) == 0 )
    for ( index = 0; index < ngrids; index++ )
      {
        if ( streamptr->xdimID[index] != UNDEFID )
          {
            gridID0 = vlistGrid(vlistID, index);
            gridtype0 = gridInqType(gridID0);
            if ( gridtype0 == GRID_GENERIC )
              {
                dimlen0 = gridInqSize(gridID0);
                if ( dimlen == dimlen0 )
                  {
                    dimID = streamptr->xdimID[index];
                    break;
                  }
                else
                  iz++; 
              }
          }
      }

  if ( gridInqXsize(gridID) == 0 )
    for ( index = 0; index < ngrids; index++ )
      {
        if ( streamptr->ydimID[index] != UNDEFID )
          {
            gridID0 = vlistGrid(vlistID, index);
            gridtype0 = gridInqType(gridID0);
            if ( gridtype0 == GRID_GENERIC )
              {
                dimlen0 = gridInqSize(gridID0);
                if ( dimlen == dimlen0 )
                  {
                    dimID = streamptr->ydimID[index];
                    break;
                  }
                else
                  iz++; 
              }
          }
      }

  if ( dimID == UNDEFID )
    {
      if ( iz == 0 ) axisname[5] = '\0';
      else           sprintf(&axisname[5], "%1d", iz+1);

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  gridindex = vlistGridIndex(vlistID, gridID);
  streamptr->xdimID[gridindex] = dimID;
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefUnstructured(int streamID, int gridID)
{
  char axisname[] = "ncells";
  char vertname[] = "nv";
  char xunits[256];
  char xlongname[256];
  char xstdname[256];
  char yunits[256];
  char ylongname[256];
  char ystdname[256];
  char xaxisname[256];
  char yaxisname[256];
  int index;
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int ngrids;
  int fileID;
  int dimlen, dimlen0;
  size_t len;
  int ncxvarid = UNDEFID, ncyvarid = UNDEFID;
  int ncbxvarid = UNDEFID, ncbyvarid = UNDEFID, ncavarid = UNDEFID;
  int nvertex, nvdimID = UNDEFID;
  int vlistID;
  int xtype = NC_DOUBLE;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( gridInqPrec(gridID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  ngrids = vlistNgrids(vlistID);

  dimlen = gridInqSize(gridID);
  gridindex = vlistGridIndex(vlistID, gridID);

  gridInqXname(gridID, xaxisname);
  gridInqXlongname(gridID, xlongname);
  gridInqXstdname(gridID, xstdname);
  gridInqXunits(gridID, xunits);
  gridInqYname(gridID, yaxisname);
  gridInqYlongname(gridID, ylongname);
  gridInqYstdname(gridID, ystdname);
  gridInqYunits(gridID, yunits);

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->xdimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_UNSTRUCTURED )
            {
              dimlen0 = gridInqSize(gridID0);
              if ( dimlen == dimlen0 )
		if ( gridInqNvertex(gridID0) == gridInqNvertex(gridID) && 
		     IS_EQUAL(gridInqXval(gridID0, 0), gridInqXval(gridID, 0)) &&
                     IS_EQUAL(gridInqXval(gridID0, dimlen-1), gridInqXval(gridID, dimlen-1)) )
		  {
		    dimID = streamptr->xdimID[index];
		    break;
		  }
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      int status;
      status = checkGridName('V', xaxisname, fileID, vlistID, gridID, ngrids, 'X');
      status = checkGridName('V', yaxisname, fileID, vlistID, gridID, ngrids, 'Y');
      status = checkGridName('D', axisname, fileID, vlistID, gridID, ngrids, 'X');
      status = checkGridName('D', vertname, fileID, vlistID, gridID, ngrids, 'X');

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      nvertex = gridInqNvertex(gridID);
      if ( nvertex > 0 ) cdf_def_dim(fileID, vertname, nvertex, &nvdimID);

      if ( gridInqXvalsPtr(gridID) )
        {
          cdf_def_var(fileID, xaxisname, (nc_type) xtype, 1, &dimID, &ncxvarid);

          if ( (len = strlen(xstdname)) )
            cdf_put_att_text(fileID, ncxvarid, "standard_name", len, xstdname);
          if ( (len = strlen(xlongname)) )
            cdf_put_att_text(fileID, ncxvarid, "long_name", len, xlongname);
          if ( (len = strlen(xunits)) )
            cdf_put_att_text(fileID, ncxvarid, "units", len, xunits);

          streamptr->ncxvarID[gridindex] = ncxvarid;

          if ( gridInqXboundsPtr(gridID) && nvdimID != UNDEFID )
            {
              int dimIDs[2];
              dimIDs[0] = dimID;
              dimIDs[1] = nvdimID;
              strcat(xaxisname, "_vertices");
              cdf_def_var(fileID, xaxisname, (nc_type) xtype, 2, dimIDs, &ncbxvarid);
              cdf_put_att_text(fileID, ncxvarid, "bounds", strlen(xaxisname), xaxisname);
            }
        }

      if ( gridInqYvalsPtr(gridID) )
        {
          cdf_def_var(fileID, yaxisname, (nc_type) xtype, 1, &dimID, &ncyvarid);

          if ( (len = strlen(ystdname)) )
            cdf_put_att_text(fileID, ncyvarid, "standard_name", len, ystdname);
          if ( (len = strlen(ylongname)) )
            cdf_put_att_text(fileID, ncyvarid, "long_name", len, ylongname);
          if ( (len = strlen(yunits)) )
            cdf_put_att_text(fileID, ncyvarid, "units", len, yunits);

          streamptr->ncyvarID[gridindex] = ncyvarid;

          if ( gridInqYboundsPtr(gridID) && nvdimID != UNDEFID )
            {
              int dimIDs[2];
              dimIDs[0] = dimID;
              dimIDs[1] = nvdimID;
              strcat(yaxisname, "_vertices");
              cdf_def_var(fileID, yaxisname, (nc_type) xtype, 2, dimIDs, &ncbyvarid);
              cdf_put_att_text(fileID, ncyvarid, "bounds", strlen(yaxisname), yaxisname);
            }
        }

      if ( gridInqAreaPtr(gridID) )
        {
          char yaxisname[] = "cell_area";
          char units[] = "m2";
          char longname[] = "area of grid cell";
          char stdname[] = "cell_area";

          cdf_def_var(fileID, yaxisname, (nc_type) xtype, 1, &dimID, &ncavarid);

          cdf_put_att_text(fileID, ncavarid, "standard_name", strlen(stdname), stdname);
          cdf_put_att_text(fileID, ncavarid, "long_name", strlen(longname), longname);
          cdf_put_att_text(fileID, ncavarid, "units", strlen(units), units);

          streamptr->ncavarID[gridindex] = ncavarid;
        }

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      if ( ncxvarid  != UNDEFID ) cdf_put_var_double(fileID, ncxvarid,  gridInqXvalsPtr(gridID));
      if ( ncbxvarid != UNDEFID ) cdf_put_var_double(fileID, ncbxvarid, gridInqXboundsPtr(gridID));
      if ( ncyvarid  != UNDEFID ) cdf_put_var_double(fileID, ncyvarid,  gridInqYvalsPtr(gridID));
      if ( ncbyvarid != UNDEFID ) cdf_put_var_double(fileID, ncbyvarid, gridInqYboundsPtr(gridID));
      if ( ncavarid  != UNDEFID ) cdf_put_var_double(fileID, ncavarid,  gridInqAreaPtr(gridID));
    }

  streamptr->xdimID[gridindex] = dimID;
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefVCT(int streamID, int zaxisID)
{
  int type;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  type = zaxisInqType(zaxisID);
  if ( type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF )
    {
      int i;
      int fileID;
      int ilev = zaxisInqVctSize(zaxisID)/2;
      int mlev = ilev - 1;
      size_t start;
      size_t count = 1;
      int ncdimid, ncdimid2;
      int hyaiid, hybiid, hyamid, hybmid;
      double mval;
      char tmpname[256];

      if ( streamptr->vct.ilev > 0 )
        {
          if ( streamptr->vct.ilev != ilev )
            Error("more than one VCT for each file unsupported!");
          return;
        }

      if ( ilev == 0 )
        {
          Warning("VCT missing");
          return;
        }

      fileID = streamInqFileID(streamID);

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, "nhym", mlev, &ncdimid);
      cdf_def_dim(fileID, "nhyi", ilev, &ncdimid2);

      streamptr->vct.mlev   = mlev;
      streamptr->vct.ilev   = ilev;
      streamptr->vct.mlevID = ncdimid;
      streamptr->vct.ilevID = ncdimid2;

      cdf_def_var(fileID, "hyai", NC_DOUBLE, 1, &ncdimid2, &hyaiid);
      cdf_def_var(fileID, "hybi", NC_DOUBLE, 1, &ncdimid2, &hybiid);
      cdf_def_var(fileID, "hyam", NC_DOUBLE, 1, &ncdimid,  &hyamid);
      cdf_def_var(fileID, "hybm", NC_DOUBLE, 1, &ncdimid,  &hybmid);

      strcpy(tmpname, "hybrid A coefficient at layer interfaces");
      cdf_put_att_text(fileID, hyaiid, "long_name", strlen(tmpname), tmpname);
      strcpy(tmpname, "Pa");
      cdf_put_att_text(fileID, hyaiid, "units", strlen(tmpname), tmpname);
      strcpy(tmpname, "hybrid B coefficient at layer interfaces");
      cdf_put_att_text(fileID, hybiid, "long_name", strlen(tmpname), tmpname);
      strcpy(tmpname, "1");
      cdf_put_att_text(fileID, hybiid, "units", strlen(tmpname), tmpname);
      strcpy(tmpname, "hybrid A coefficient at layer midpoints");
      cdf_put_att_text(fileID, hyamid, "long_name", strlen(tmpname), tmpname);
      strcpy(tmpname, "Pa");
      cdf_put_att_text(fileID, hyamid, "units", strlen(tmpname), tmpname);
      strcpy(tmpname, "hybrid B coefficient at layer midpoints");
      cdf_put_att_text(fileID, hybmid, "long_name", strlen(tmpname), tmpname);
      strcpy(tmpname, "1");
      cdf_put_att_text(fileID, hybmid, "units", strlen(tmpname), tmpname);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      const double *vctptr = zaxisInqVctPtr(zaxisID);

      cdf_put_var_double(fileID, hyaiid, vctptr);
      cdf_put_var_double(fileID, hybiid, vctptr+ilev);

      for ( i = 0; i < mlev; i++ )
        {
          start = i;
          mval = (vctptr[i] + vctptr[i+1]) * 0.5;
          cdf_put_vara_double(fileID, hyamid, &start, &count, &mval);
          mval = (vctptr[ilev+i] + vctptr[ilev+i+1]) * 0.5;
          cdf_put_vara_double(fileID, hybmid, &start, &count, &mval);
        }
    }
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefZaxis(int streamID, int zaxisID)
{
  /*  char zaxisname0[256]; */
  char axisname[256];
  char stdname[256];
  char longname[256];
  char units[256];
  char tmpname[256];
  int index;
  int zaxisID0;
  int dimID = UNDEFID;
  int dimIDs[2];
  int fileID;
  int dimlen;
  size_t len;
  int ncvarid = UNDEFID, ncbvarid = UNDEFID;
  int nvertex = 2, nvdimID = UNDEFID;
  int type;
  int nzaxis;
  int ilevel = 0;
  int vlistID;
  int zaxisindex;
  int xtype = NC_DOUBLE;
  int positive;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( zaxisInqPrec(zaxisID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  zaxisindex = vlistZaxisIndex(vlistID, zaxisID);

  nzaxis = vlistNzaxis(vlistID);

  dimlen = zaxisInqSize(zaxisID);
  type   = zaxisInqType(zaxisID);

  if ( dimlen == 1 && type == ZAXIS_SURFACE     ) return;
  if ( dimlen == 1 && type == ZAXIS_TOA         ) return;
  if ( dimlen == 1 && type == ZAXIS_SEA_BOTTOM  ) return;
  if ( dimlen == 1 && type == ZAXIS_ATMOSPHERE  ) return;
  if ( dimlen == 1 && type == ZAXIS_MEANSEA     ) return;

  zaxisInqName(zaxisID, axisname);
  /*
  for ( index = 0; index < nzaxis; index++ )
    {
      if ( streamptr->zaxisID[index] != UNDEFID )
        {
          zaxisID0 = vlistZaxis(vlistID, index);
          zaxisInqName(zaxisID0, zaxisname0);
          if ( strcmp(zaxisname0, axisname) == 0 ) ilevel++;
        }
    }
  */
  if ( dimID == UNDEFID )
    {
      char axisname0[256];
      char axisname2[256];
      int checkname = FALSE;
      int status;

      /* check that the name is not already defined */
      checkname = TRUE;
      ilevel = 0;

      while ( checkname ) 
        {
          strcpy(axisname2, axisname);
          if ( ilevel ) sprintf(&axisname2[strlen(axisname2)], "_%d", ilevel+1);

          status = nc_inq_varid(fileID, axisname2, &ncvarid);
          if ( status != NC_NOERR )
            {
              if ( ilevel )
                {
                  /* check that the name does not exist for other grids */
                  for ( index = 0; index < nzaxis; index++ )
                    {
                      zaxisID0 = vlistZaxis(vlistID, index);
                      if ( zaxisID != zaxisID0 )
                        {
                          zaxisInqName(zaxisID0, axisname0);
                          if ( strcmp(axisname0, axisname2) == 0 ) break;
                        }
                    }
                  if ( index == nzaxis ) checkname = FALSE;
                }
              else
                {
                  checkname = FALSE;
                }
            }

          if ( checkname ) ilevel++;

          if ( ilevel > 99 ) break;
        }

      if ( ilevel ) sprintf(&axisname[strlen(axisname)], "_%1d", ilevel+1);

      if ( type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF )
        {
          if ( type == ZAXIS_HYBRID )
            {
	      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

	      cdf_def_dim(fileID, axisname, dimlen, &dimID);
	      cdf_def_var(fileID, axisname, (nc_type) xtype, 1, &dimID,  &ncvarid);

	      strcpy(tmpname, "hybrid_sigma_pressure");
	      cdf_put_att_text(fileID, ncvarid, "standard_name", strlen(tmpname), tmpname);
	      strcpy(tmpname, "hybrid level at layer midpoints");
	      cdf_put_att_text(fileID, ncvarid, "long_name", strlen(tmpname), tmpname);
	      strcpy(tmpname, "level");
	      cdf_put_att_text(fileID, ncvarid, "units", strlen(tmpname), tmpname);
	      strcpy(tmpname, "down");
	      cdf_put_att_text(fileID, ncvarid, "positive", strlen(tmpname), tmpname);
	      strcpy(tmpname, "hyam hybm (mlev=hyam+hybm*aps)");
	      cdf_put_att_text(fileID, ncvarid, "formula", strlen(tmpname), tmpname);
	      strcpy(tmpname, "ap: hyam b: hybm ps: aps");
	      cdf_put_att_text(fileID, ncvarid, "formula_terms", strlen(tmpname), tmpname);
	      /*
	      strcpy(tmpname, "ilev");
	      cdf_put_att_text(fileID, ncvarid, "borders", strlen(tmpname), tmpname);
	      */
	      cdf_enddef(fileID);
	      streamptr->ncmode = 2;

	      cdf_put_var_double(fileID, ncvarid, zaxisInqLevelsPtr(zaxisID));
            }

          if ( type == ZAXIS_HYBRID_HALF )
            {
	      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

	      cdf_def_dim(fileID, axisname, dimlen, &dimID);
	      cdf_def_var(fileID, axisname, (nc_type) xtype, 1, &dimID,  &ncvarid);

	      strcpy(tmpname, "hybrid_sigma_pressure");
	      cdf_put_att_text(fileID, ncvarid, "standard_name", strlen(tmpname), tmpname);
	      strcpy(tmpname, "hybrid level at layer interfaces");
	      cdf_put_att_text(fileID, ncvarid, "long_name", strlen(tmpname), tmpname);
	      strcpy(tmpname, "level");
	      cdf_put_att_text(fileID, ncvarid, "units", strlen(tmpname), tmpname);
	      strcpy(tmpname, "down");
	      cdf_put_att_text(fileID, ncvarid, "positive", strlen(tmpname), tmpname);
	      strcpy(tmpname, "hyai hybi (ilev=hyai+hybi*aps)");
	      cdf_put_att_text(fileID, ncvarid, "formula", strlen(tmpname), tmpname);
	      strcpy(tmpname, "ap: hyai b: hybi ps: aps");
	      cdf_put_att_text(fileID, ncvarid, "formula_terms", strlen(tmpname), tmpname);

	      cdf_enddef(fileID);
	      streamptr->ncmode = 2;

	      cdf_put_var_double(fileID, ncvarid, zaxisInqLevelsPtr(zaxisID));
            }

          cdfDefVCT(streamID, zaxisID);

          if ( dimID == UNDEFID )
            {
              if ( type == ZAXIS_HYBRID )
                streamptr->zaxisID[zaxisindex] = streamptr->vct.mlevID;
              else
                streamptr->zaxisID[zaxisindex] = streamptr->vct.ilevID;
            }
        }
      else
        {
          if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

          cdf_def_dim(fileID, axisname, dimlen, &dimID);

          zaxisInqLongname(zaxisID, longname);
          zaxisInqUnits(zaxisID, units);
          zaxisInqStdname(zaxisID, stdname);

          cdf_def_var(fileID, axisname, (nc_type) xtype, 1, &dimID, &ncvarid);

          if ( (len = strlen(stdname)) )
            cdf_put_att_text(fileID, ncvarid, "standard_name", len, stdname);
          if ( (len = strlen(longname)) )
            cdf_put_att_text(fileID, ncvarid, "long_name", len, longname);
          if ( (len = strlen(units)) )
            cdf_put_att_text(fileID, ncvarid, "units", len, units);

	  positive = zaxisInqPositive(zaxisID);
	  if ( positive == 1 )
	    {
	      strcpy(tmpname, "up");
	      cdf_put_att_text(fileID, ncvarid, "positive", strlen(tmpname), tmpname);
	    }
	  else if ( positive == 2 )
	    {
	      strcpy(tmpname, "down");
	      cdf_put_att_text(fileID, ncvarid, "positive", strlen(tmpname), tmpname);
	    }

          cdf_put_att_text(fileID, ncvarid, "axis", 1, "Z");

	  if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
            {
	      if ( nc_inq_dimid(fileID, "nb2", &nvdimID) != NC_NOERR )
		cdf_def_dim(fileID, "nb2", nvertex, &nvdimID);

	      if ( nvdimID != UNDEFID )
		{
		  strcat(axisname, "_bnds");
		  dimIDs[0] = dimID;
		  dimIDs[1] = nvdimID;
		  cdf_def_var(fileID, axisname, (nc_type) xtype, 2, dimIDs, &ncbvarid);
		  cdf_put_att_text(fileID, ncvarid, "bounds", strlen(axisname), axisname);
		}
	    }

          cdf_enddef(fileID);
          streamptr->ncmode = 2;

          cdf_put_var_double(fileID, ncvarid, zaxisInqLevelsPtr(zaxisID));

          if ( ncbvarid != UNDEFID )
	    {
	      int i;
	      double *zbounds, *lbounds, *ubounds;

	      lbounds = (double *) malloc(dimlen*sizeof(double));
	      ubounds = (double *) malloc(dimlen*sizeof(double));
	      zbounds = (double *) malloc(2*dimlen*sizeof(double));

	      zaxisInqLbounds(zaxisID, lbounds);
	      zaxisInqUbounds(zaxisID, ubounds);

	      for ( i = 0; i < dimlen; ++i )
		{
		  zbounds[2*i  ] = lbounds[i];
		  zbounds[2*i+1] = ubounds[i];
		}

	      cdf_put_var_double(fileID, ncbvarid, zbounds);

	      free(zbounds);
	      free(ubounds);
	      free(lbounds);
	    }
        }
    }

  if ( dimID != UNDEFID )
    streamptr->zaxisID[zaxisindex] = dimID;
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefPole(int streamID, int gridID)
{
  int fileID;
  int ncvarid = UNDEFID;
  int ncerr;
  double xpole, ypole, angle;
  char varname[] = "rotated_pole";
  char mapname[] = "rotated_latitude_longitude";

  fileID  = streamInqFileID(streamID);

  ypole = gridInqYpole(gridID);
  xpole = gridInqXpole(gridID);
  angle = gridInqAngle(gridID);

  cdf_redef(fileID);

  ncerr = nc_def_var(fileID, varname, (nc_type) NC_CHAR, 0, NULL, &ncvarid);
  if ( ncerr == NC_NOERR )
    {
      cdf_put_att_text(fileID, ncvarid, "grid_mapping_name", strlen(mapname), mapname);
      cdf_put_att_double(fileID, ncvarid, "grid_north_pole_latitude", NC_DOUBLE, 1, &ypole);
      cdf_put_att_double(fileID, ncvarid, "grid_north_pole_longitude", NC_DOUBLE, 1, &xpole);
      if ( angle > 0 )
        cdf_put_att_double(fileID, ncvarid, "north_pole_grid_longitude", NC_DOUBLE, 1, &angle);
    }

  cdf_enddef(fileID);
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefMapping(int streamID, int gridID)
{
  int fileID;
  int ncvarid = UNDEFID;
  int ncerr;

  if ( gridInqType(gridID) == GRID_SINUSOIDAL )
    {
      char varname[] = "sinusoidal";
      char mapname[] = "sinusoidal";

      fileID  = streamInqFileID(streamID);

      cdf_redef(fileID);

      ncerr = nc_def_var(fileID, varname, (nc_type) NC_CHAR, 0, NULL, &ncvarid);
      if ( ncerr == NC_NOERR )
        {
          cdf_put_att_text(fileID, ncvarid, "grid_mapping_name", strlen(mapname), mapname);
          /*
          cdf_put_att_double(fileID, ncvarid, "grid_north_pole_latitude", NC_DOUBLE, 1, &ypole);
          cdf_put_att_double(fileID, ncvarid, "grid_north_pole_longitude", NC_DOUBLE, 1, &xpole);
          */
        }

      cdf_enddef(fileID);
    }
  else if ( gridInqType(gridID) == GRID_LAEA )
    {
      char varname[] = "laea";
      char mapname[] = "lambert_azimuthal_equal_area";

      fileID  = streamInqFileID(streamID);

      cdf_redef(fileID);

      ncerr = nc_def_var(fileID, varname, (nc_type) NC_CHAR, 0, NULL, &ncvarid);
      if ( ncerr == NC_NOERR )
        {
          double a, lon_0, lat_0;

          gridInqLaea(gridID, &a, &lon_0, &lat_0);

          cdf_put_att_text(fileID, ncvarid, "grid_mapping_name", strlen(mapname), mapname);
          cdf_put_att_double(fileID, ncvarid, "earth_radius", NC_DOUBLE, 1, &a);
          cdf_put_att_double(fileID, ncvarid, "longitude_of_projection_origin", NC_DOUBLE, 1, &lon_0);
          cdf_put_att_double(fileID, ncvarid, "latitude_of_projection_origin", NC_DOUBLE, 1, &lat_0);
        }

      cdf_enddef(fileID);
    }
  else if ( gridInqType(gridID) == GRID_LCC2 )
    {
      char varname[] = "Lambert_Conformal";
      char mapname[] = "lambert_conformal_conic";

      fileID  = streamInqFileID(streamID);

      cdf_redef(fileID);

      ncerr = nc_def_var(fileID, varname, (nc_type) NC_CHAR, 0, NULL, &ncvarid);
      if ( ncerr == NC_NOERR )
        {
          double radius, lon_0, lat_0, lat_1, lat_2;

          gridInqLcc2(gridID, &radius, &lon_0, &lat_0, &lat_1, &lat_2);

          cdf_put_att_text(fileID, ncvarid, "grid_mapping_name", strlen(mapname), mapname);
          if ( radius > 0 )
            cdf_put_att_double(fileID, ncvarid, "earth_radius", NC_DOUBLE, 1, &radius);
          cdf_put_att_double(fileID, ncvarid, "longitude_of_central_meridian", NC_DOUBLE, 1, &lon_0);
          cdf_put_att_double(fileID, ncvarid, "latitude_of_projection_origin", NC_DOUBLE, 1, &lat_0);
          if ( IS_EQUAL(lat_1, lat_2) )
            cdf_put_att_double(fileID, ncvarid, "standard_parallel", NC_DOUBLE, 1, &lat_1);
          else
            {
              double lat_1_2[2];
              lat_1_2[0] = lat_1;
              lat_1_2[1] = lat_2;
              cdf_put_att_double(fileID, ncvarid, "standard_parallel", NC_DOUBLE, 2, lat_1_2);
            }
        }

      cdf_enddef(fileID);
    }
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfDefGrid(int streamID, int gridID)
{
  int gridtype, size;
  int gridindex;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  gridindex = vlistGridIndex(vlistID, gridID);
  if ( streamptr->xdimID[gridindex] != UNDEFID ) return;

  gridtype = gridInqType(gridID);
  size     = gridInqSize(gridID);

  if ( CDI_Debug )
    Message("gridtype = %d  size = %d", gridtype, size);

  if ( gridtype == GRID_GAUSSIAN ||
       gridtype == GRID_LONLAT   ||
       gridtype == GRID_GENERIC )
    {
      if ( gridtype == GRID_GENERIC && size == 1 && 
           gridInqXsize(gridID) == 0 && gridInqYsize(gridID) == 0 )
        {
          /* no grid information */
        }
      else if ( gridtype == GRID_GENERIC && (gridInqXsize(gridID) == 0 || gridInqYsize(gridID) == 0) )
        {
          cdfDefGdim(streamID, gridID);
        }
      else
        {
          if ( gridInqXsize(gridID) > 0 ) cdfDefXaxis(streamID, gridID);
          if ( gridInqYsize(gridID) > 0 ) cdfDefYaxis(streamID, gridID);
        }

      if ( gridIsRotated(gridID) ) cdfDefPole(streamID, gridID);
    }
  else if ( gridtype == GRID_CURVILINEAR )
    {
      cdfDefCurvilinear(streamID, gridID);
    }
  else if ( gridtype == GRID_UNSTRUCTURED )
    {
      cdfDefUnstructured(streamID, gridID);
    }
  else if ( gridtype == GRID_GAUSSIAN_REDUCED )
    {
      cdfDefRgrid(streamID, gridID);
    }
  else if ( gridtype == GRID_SPECTRAL )
    {
      cdfDefComplex(streamID, gridID);
      cdfDefSP(streamID, gridID);
    }
  else if ( gridtype == GRID_FOURIER )
    {
      cdfDefComplex(streamID, gridID);
      cdfDefFC(streamID, gridID);
    }
  else if ( gridtype == GRID_TRAJECTORY )
    {
      cdfDefTrajLon(streamID, gridID);
      cdfDefTrajLat(streamID, gridID);
    }
  else if ( gridtype == GRID_SINUSOIDAL || gridtype == GRID_LAEA || gridtype == GRID_LCC2 )
    {
      cdfDefXaxis(streamID, gridID);
      cdfDefYaxis(streamID, gridID);

      cdfDefMapping(streamID, gridID);
    }
  /*
  else if ( gridtype == GRID_LCC )
    {
      cdfDefLcc(streamID, gridID);
    }
  */
  else
    {
      Error("Unsupported grid type: %s", gridNamePtr(gridtype));
    }
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
int cdfDefVar(int streamID, int varID)
{
  int ncvarid = -1;
  int fileID;
  int xid = UNDEFID, yid = UNDEFID, zid = UNDEFID, tid = UNDEFID;
  size_t xsize = 0, ysize = 0;
  int code, param, gridID, zaxisID;
  int pnum, pcat, pdis;
  char varname[256];
  const char *name = NULL;
  const char *longname = NULL;
  const char *stdname = NULL;
  const char *units = NULL;
  int dims[4];
  int lchunk = FALSE;
  size_t chunks[4] = {0,0,0,0};
  int tableID;
  int ndims = 0;
  int len;
  int timeID;
  int xtype, dtype;
  int gridtype, gridsize;
  int gridindex, zaxisindex;
  int tablenum;
  int vlistID;
  int dimorder[3];
  int ixyz;
  int iax = 0;
  char axis[5];
  stream_t *streamptr;
  int retval;

  streamptr = stream_to_pointer(streamID);

  fileID  = streamInqFileID(streamID);

  if ( CDI_Debug )
    Message("streamID = %d, fileID = %d, varID = %d", streamID, fileID, varID);

  if ( streamptr->vars[varID].ncvarid != UNDEFID )
    return (streamptr->vars[varID].ncvarid);

  vlistID = streamInqVlist(streamID);
  gridID  = vlistInqVarGrid(vlistID, varID);
  zaxisID = vlistInqVarZaxis(vlistID, varID);
  timeID  = vlistInqVarTime(vlistID, varID);
  code    = vlistInqVarCode(vlistID, varID);
  param   = vlistInqVarParam(vlistID, varID);
  cdiDecodeParam(param, &pnum, &pcat, &pdis);

  ixyz    = vlistInqVarXYZ(vlistID, varID);
  if ( ixyz == 0 ) ixyz = 321; // ZYX

  gridsize  = gridInqSize(gridID);
  if ( gridsize > 1 ) lchunk = TRUE;
  gridtype  = gridInqType(gridID);
  gridindex = vlistGridIndex(vlistID, gridID);
  if ( gridtype != GRID_TRAJECTORY )
    {
      xid = streamptr->xdimID[gridindex];
      yid = streamptr->ydimID[gridindex];
      if ( xid != UNDEFID ) cdf_inq_dimlen(fileID, xid, &xsize);
      if ( yid != UNDEFID ) cdf_inq_dimlen(fileID, yid, &ysize);
    }

  zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  zid = streamptr->zaxisID[zaxisindex];

  dimorder[0] = ixyz/100;
  dimorder[1] = (ixyz-dimorder[0]*100)/10;
  dimorder[2] = (ixyz-dimorder[0]*100-dimorder[1]*10);
  if ( dimorder[0] != 3 ) lchunk = FALSE; /* ZYX and ZXY */

  if ( ((dimorder[0]>0)+(dimorder[1]>0)+(dimorder[2]>0)) < ((xid!=UNDEFID)+(yid!=UNDEFID)+(zid!=UNDEFID)) )
    {
      printf("xyz=%d  zid=%d  yid=%d  xid=%d\n", ixyz, zid, yid, xid);
      Error("Internal problem, dimension order missing!");
    }

  tid = streamptr->basetime.ncdimid;

  if ( timeID == TIME_VARIABLE )
    {
      if ( tid == UNDEFID ) Error("Internal problem, time undefined!");
      chunks[ndims] = 1;
      dims[ndims++] = tid;
      axis[iax++] = 'T';
    }
  /*
  if ( zid != UNDEFID ) axis[iax++] = 'Z';
  if ( zid != UNDEFID ) chunks[ndims] = 1;
  if ( zid != UNDEFID ) dims[ndims++] = zid;

  if ( yid != UNDEFID ) chunks[ndims] = ysize;
  if ( yid != UNDEFID ) dims[ndims++] = yid;

  if ( xid != UNDEFID ) chunks[ndims] = xsize;
  if ( xid != UNDEFID ) dims[ndims++] = xid;
  */
  for ( int id = 0; id < 3; ++id )
    {
      if ( dimorder[id] == 3 && zid != UNDEFID )
        {
          axis[iax++] = 'Z';
          chunks[ndims] = 1;
          dims[ndims] = zid;
          ndims++;
        }
      else if ( dimorder[id] == 2 && yid != UNDEFID )
        {
          chunks[ndims] = ysize;
          dims[ndims] = yid;
          ndims++;
        }
      else if ( dimorder[id] == 1 && xid != UNDEFID )
        {
          chunks[ndims] = xsize;
          dims[ndims] = xid;
          ndims++;
        }
    }

  if ( CDI_Debug )
    fprintf(stderr, "chunks %d %d %d %d\n",
            (int)chunks[0], (int)chunks[1], (int)chunks[2], (int)chunks[3]);

  tableID  = vlistInqVarTable(vlistID, varID);

  name     = vlistInqVarNamePtr(vlistID, varID);
  longname = vlistInqVarLongnamePtr(vlistID, varID);
  stdname  = vlistInqVarStdnamePtr(vlistID, varID);
  units    = vlistInqVarUnitsPtr(vlistID, varID);

  if ( name     == NULL )     name = tableInqParNamePtr(tableID, code);
  if ( longname == NULL ) longname = tableInqParLongnamePtr(tableID, code);
  if ( units    == NULL )    units = tableInqParUnitsPtr(tableID, code);
  if ( name )
    {
      int checkname;
      int iz;
      int status;

      sprintf(varname, "%s", name);

      checkname = TRUE;
      iz = 0;

      while ( checkname ) 
        {
          if ( iz ) sprintf(varname, "%s_%d", name, iz+1);

          status = nc_inq_varid(fileID, varname, &ncvarid);
          if ( status != NC_NOERR )
            {
              checkname = FALSE;
            }

          if ( checkname ) iz++;

          if ( iz > 99 )
            Error("Double entry of variable name '%s'!", name);
        }

      if ( strcmp(name, varname) != 0 )
        {
          if ( iz == 1 )
            Warning("Changed double entry of variable name '%s' to '%s'!", name, varname);
          else
            Warning("Changed multiple entry of variable name '%s' to '%s'!", name, varname);
        }

      name = varname;
    }
  else
    {
      int checkname;
      int iz;
      int status;
      char *varname2;

      if ( code < 0 ) code = -code;
      if ( pnum < 0 ) pnum = -pnum;

      if ( pdis == 255 )
	sprintf(varname, "var%d", code);
      else
	sprintf(varname, "param%d.%d.%d", pnum, pcat, pdis);

      varname2 = varname+strlen(varname);

      checkname = TRUE;
      iz = 0;

      while ( checkname )
        {
          if ( iz ) sprintf(varname2, "_%d", iz+1);

          status = nc_inq_varid(fileID, varname, &ncvarid);
          if ( status != NC_NOERR ) checkname = FALSE;

          if ( checkname ) iz++;

          if ( iz > 99 ) break;
        }

      name = varname;
      code = 0;
      pdis = 255;
    }

  /* if ( streamptr->ncmode == 2 ) cdf_redef(fileID); */

  dtype = vlistInqVarDatatype(vlistID, varID);
  xtype = cdfDefDatatype(dtype, streamptr->filetype);

  cdf_def_var(fileID, name, (nc_type) xtype, ndims, dims, &ncvarid);

#if  defined  (HAVE_NETCDF4)
  if ( lchunk && (streamptr->filetype == FILETYPE_NC4 || streamptr->filetype == FILETYPE_NC4C) )
    {
      if ( (retval = nc_def_var_chunking(fileID, ncvarid, 0, chunks)) )
        Error("nc_def_var_chunking failed, status = %d", retval);
    }
#endif

  if ( streamptr->comptype == COMPRESS_ZIP )
    {
      if ( lchunk && (streamptr->filetype == FILETYPE_NC4 || streamptr->filetype == FILETYPE_NC4C) )
        {
          cdfDefVarDeflate(fileID, ncvarid, streamptr->complevel);
        }
      else
        {
          if ( lchunk )
            {
              static int lwarn = TRUE;

              if ( lwarn )
                {
                  lwarn = FALSE;
                  Warning("Deflate compression is only available for netCDF4!");
                }
            }
        }
    }

  if ( streamptr->comptype == COMPRESS_SZIP )
    {
      if ( lchunk && (streamptr->filetype == FILETYPE_NC4 || streamptr->filetype == FILETYPE_NC4C) )
        {
#if defined (NC_SZIP_NN_OPTION_MASK)
          cdfDefVarSzip(fileID, ncvarid);
#else
          static int lwarn = TRUE;

          if ( lwarn )
            {
              lwarn = FALSE;
              Warning("netCDF4/SZIP compression not available!");
            }
#endif
        }
      else
        {
          static int lwarn = TRUE;

          if ( lwarn )
            {
              lwarn = FALSE;
              Warning("SZIP compression is only available for netCDF4!");
            }
        }
    }

  if ( stdname && *stdname )
    cdf_put_att_text(fileID, ncvarid, "standard_name", strlen(stdname), stdname);

  if ( longname && *longname )
    cdf_put_att_text(fileID, ncvarid, "long_name", strlen(longname), longname);

  if ( units && *units )
    cdf_put_att_text(fileID, ncvarid, "units", strlen(units), units);

  if ( code > 0 && pdis == 255 )
    cdf_put_att_int(fileID, ncvarid, "code", NC_INT, 1, &code);

  if ( pdis != 255 )
    {
      char paramstr[32];
      cdiParamToString(param, paramstr, sizeof(paramstr));
      cdf_put_att_text(fileID, ncvarid, "param", strlen(paramstr), paramstr);
    }

  if ( tableID != UNDEFID )
    {
      tablenum = tableInqNum(tableID);
      if ( tablenum > 0 )
        cdf_put_att_int(fileID, ncvarid, "table", NC_INT, 1, &tablenum);
    }

  if ( gridtype != GRID_GENERIC && gridtype != GRID_LONLAT  &&
       gridtype != GRID_CURVILINEAR )
    {
      len = strlen(gridNamePtr(gridtype));
      if ( len > 0 )
        cdf_put_att_text(fileID, ncvarid, "grid_type", len, gridNamePtr(gridtype));
    }

  if ( gridIsRotated(gridID) )
    {
      char mapping[] = "rotated_pole";
      cdf_put_att_text(fileID, ncvarid, "grid_mapping", strlen(mapping), mapping);
    }

  if ( gridtype == GRID_SINUSOIDAL )
    {
      char mapping[] = "sinusoidal";
      cdf_put_att_text(fileID, ncvarid, "grid_mapping", strlen(mapping), mapping);
    }
  else if ( gridtype == GRID_LAEA )
    {
      char mapping[] = "laea";
      cdf_put_att_text(fileID, ncvarid, "grid_mapping", strlen(mapping), mapping);
    }
  else if ( gridtype == GRID_LCC2 )
    {
      char mapping[] = "Lambert_Conformal";
      cdf_put_att_text(fileID, ncvarid, "grid_mapping", strlen(mapping), mapping);
    }
  else if ( gridtype == GRID_TRAJECTORY )
    {
      cdf_put_att_text(fileID, ncvarid, "coordinates", 9, "tlon tlat" );
    }
  else if ( gridtype == GRID_UNSTRUCTURED || gridtype == GRID_CURVILINEAR )
    {
      char coordinates[256] = "";
      char cellarea[256] = "area: ";
      int ncxvarID, ncyvarID, ncavarID;
      int gridindex;
      size_t len;

      gridindex = vlistGridIndex(vlistID, gridID);
      ncxvarID = streamptr->ncxvarID[gridindex];
      ncyvarID = streamptr->ncyvarID[gridindex];
      ncavarID = streamptr->ncavarID[gridindex];
      if ( ncxvarID != CDI_UNDEFID )
        cdf_inq_varname(fileID, ncxvarID, coordinates);
      len = strlen(coordinates);
      if ( ncyvarID != CDI_UNDEFID )
        {
          if ( len ) coordinates[len++] = ' ';
          cdf_inq_varname(fileID, ncyvarID, coordinates+len);
        }
      len = strlen(coordinates);
      if ( len )
        cdf_put_att_text(fileID, ncvarid, "coordinates", len, coordinates);

      if ( ncavarID != CDI_UNDEFID )
        {
          len = strlen(cellarea);
          cdf_inq_varname(fileID, ncavarID, cellarea+len);
          len = strlen(cellarea);
          cdf_put_att_text(fileID, ncvarid, "cell_measures", len, cellarea);
        }
    }
  else if ( gridtype == GRID_SPECTRAL || gridtype == GRID_FOURIER )
    {
      int gridTruncation = gridInqTrunc(gridID);

      axis[iax++] = '-';
      axis[iax++] = '-';
      cdf_put_att_text(fileID, ncvarid, "axis", iax, axis);
      cdf_put_att_int(fileID, ncvarid, "truncation", NC_INT, 1, &gridTruncation);
    }

  /*  if ( xtype == NC_BYTE || xtype == NC_SHORT || xtype == NC_INT ) */
    {
      int laddoffset, lscalefactor;
      double addoffset, scalefactor;
      int astype = NC_DOUBLE;

      addoffset    = vlistInqVarAddoffset(vlistID, varID);
      scalefactor  = vlistInqVarScalefactor(vlistID, varID);
      laddoffset   = IS_NOT_EQUAL(addoffset, 0);
      lscalefactor = IS_NOT_EQUAL(scalefactor, 1);

      if ( laddoffset || lscalefactor )
        {
          if ( IS_EQUAL(addoffset,   (double) ((float) addoffset)) &&
               IS_EQUAL(scalefactor, (double) ((float) scalefactor)) )
            {
              astype = NC_FLOAT;
            }

          if ( xtype == (int) NC_FLOAT ) astype = NC_FLOAT;

          cdf_put_att_double(fileID, ncvarid, "add_offset",   (nc_type) astype, 1, &addoffset);
          cdf_put_att_double(fileID, ncvarid, "scale_factor", (nc_type) astype, 1, &scalefactor);
        }
    }

  if ( dtype == DATATYPE_UINT8 && xtype == NC_BYTE )
    {
      int validrange[2] = {0, 255};
      cdf_put_att_int(fileID, ncvarid, "valid_range", NC_SHORT, 2, validrange);
      cdf_put_att_text(fileID, ncvarid, "_Unsigned", 4, "true");
    }

  streamptr->vars[varID].ncvarid = ncvarid;

  if ( vlistInqVarMissvalUsed(vlistID, varID) )
    cdfDefVarMissval(streamID, varID, vlistInqVarDatatype(vlistID, varID), 0);

  if ( zid == -1 )
    {
      if ( zaxisInqType(zaxisID) == ZAXIS_TOA         || 
           zaxisInqType(zaxisID) == ZAXIS_SEA_BOTTOM  ||
           zaxisInqType(zaxisID) == ZAXIS_ATMOSPHERE )
        {
          zaxisInqName(zaxisID, varname);
          cdf_put_att_text(fileID, ncvarid, "level_type", strlen(varname), varname);
        }
    }

  /* Attributes */
  defineAttributes(vlistID, varID, fileID, ncvarid);

  /* if ( streamptr->ncmode == 2 ) cdf_enddef(fileID); */

  return (ncvarid);
}
#endif

static
void scale_add(long size, double *data, double addoffset, double scalefactor)
{
  long i;
  int laddoffset;
  int lscalefactor;

  laddoffset   = IS_NOT_EQUAL(addoffset, 0);
  lscalefactor = IS_NOT_EQUAL(scalefactor, 1);

  if ( laddoffset || lscalefactor )
    {
      for ( i = 0; i < size; ++i )
        {
          if ( lscalefactor ) data[i] *= scalefactor;
          if ( laddoffset )   data[i] += addoffset;
        }
    }
}


void cdfReadVarDP(int streamID, int varID, double *data, int *nmiss)
{
#if  defined  (HAVE_LIBNETCDF)
  int fileID;
  int gridID;
  int zaxisID;
  int xid = UNDEFID, yid = UNDEFID, zid = UNDEFID;
  int ncvarid;
  int tsID;
  size_t size;
  size_t start[4];
  size_t count[4];
  int ndims = 0;
  int idim;
  int timeID;
  int gridindex, zaxisindex;
  int vlistID;
  int i;
  double missval;
  int laddoffset, lscalefactor;
  double addoffset, scalefactor;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( CDI_Debug )
    Message("streamID = %d  varID = %d", streamID, varID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  tsID = streamptr->curTsID;

  if ( CDI_Debug ) Message("tsID = %d", tsID);

  ncvarid = streamptr->vars[varID].ncvarid;

  gridID  = vlistInqVarGrid(vlistID, varID);
  zaxisID = vlistInqVarZaxis(vlistID, varID);
  timeID  = vlistInqVarTime(vlistID, varID);

  gridindex = vlistGridIndex(vlistID, gridID);
  if ( gridInqType(gridID) == GRID_TRAJECTORY )
    {
      cdfReadGridTraj(streamID, gridID);
    }
  else
    {
      xid = streamptr->xdimID[gridindex];
      yid = streamptr->ydimID[gridindex];
    }

  zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  zid = streamptr->zaxisID[zaxisindex];

  if ( timeID == TIME_VARIABLE )
    {
      start[ndims] = tsID;
      count[ndims] = 1;
      ndims++;
    }
  if ( zid != UNDEFID )
    {
      start[ndims] = 0;
      count[ndims] = zaxisInqSize(zaxisID);
      ndims++;
    }
  if ( yid != UNDEFID )
    {
      start[ndims] = 0;
      count[ndims] = gridInqYsize(gridID);
      ndims++;
    }
  if ( xid != UNDEFID )
    {
      start[ndims] = 0;
      count[ndims] = gridInqXsize(gridID);
      ndims++;
    }

  if ( CDI_Debug )
    for (idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);

  cdf_get_vara_double(fileID, ncvarid, start, count, data);

  *nmiss = 0;
  if ( vlistInqVarMissvalUsed(vlistID, varID) == TRUE  )
    {
      size    = gridInqSize(gridID)*zaxisInqSize(zaxisID);
      missval = vlistInqVarMissval(vlistID, varID);

      for ( i = 0; i < (int) size; i++ )
        if ( DBL_IS_EQUAL(data[i], missval) ) *nmiss += 1;
    }

  addoffset    = vlistInqVarAddoffset(vlistID, varID);
  scalefactor  = vlistInqVarScalefactor(vlistID, varID);
  laddoffset   = IS_NOT_EQUAL(addoffset, 0);
  lscalefactor = IS_NOT_EQUAL(scalefactor, 1);

  if ( laddoffset || lscalefactor )
    {
      size    = gridInqSize(gridID)*zaxisInqSize(zaxisID);
      missval = vlistInqVarMissval(vlistID, varID);

      if ( *nmiss > 0 )
        {
          for ( i = 0; i < (int) size; i++ )
            {
              if ( !DBL_IS_EQUAL(data[i], missval) )
                {
                  if ( lscalefactor ) data[i] *= scalefactor;
                  if ( laddoffset )   data[i] += addoffset;
                }
            }
        }
      else
        {
          for ( i = 0; i < (int) size; i++ )
            {
              if ( lscalefactor ) data[i] *= scalefactor;
              if ( laddoffset )   data[i] += addoffset;
            }
        }
    }
#endif
}


void cdfWriteVarDP(int streamID, int varID, const double *data, int nmiss)
{
#if  defined  (HAVE_LIBNETCDF)
  int fileID;
  int gridID;
  int zaxisID;
  int xid = UNDEFID, yid = UNDEFID, zid = UNDEFID;
  int ncvarid;
  int ntsteps;
  size_t size;
  size_t start[4];
  size_t count[4];
  int ndims = 0;
  int idim;
  int timeID;
  int gridindex, zaxisindex;
  int i;
  int dtype;
  int vlistID;
  double *mdata = NULL;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( CDI_Debug )
    Message("streamID = %d  varID = %d", streamID, varID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  ntsteps = streamptr->ntsteps;
  if ( CDI_Debug )
    Message("ntsteps = %d", ntsteps); 

  if ( vlistHasTime(vlistID) ) cdfDefTime(streamID);

  ncvarid = cdfDefVar(streamID, varID);

  gridID  = vlistInqVarGrid(vlistID, varID);
  zaxisID = vlistInqVarZaxis(vlistID, varID);
  timeID  = vlistInqVarTime(vlistID, varID);

  gridindex = vlistGridIndex(vlistID, gridID);
  if ( gridInqType(gridID) == GRID_TRAJECTORY )
    {
      cdfWriteGridTraj(streamID, gridID);
    }
  else
    {
      xid = streamptr->xdimID[gridindex];
      yid = streamptr->ydimID[gridindex];
    }

  zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  zid = streamptr->zaxisID[zaxisindex];

  if ( timeID == TIME_VARIABLE )
    {
      start[ndims] = ntsteps - 1;
      count[ndims] = 1;
      ndims++;
    }
  if ( zid != UNDEFID )
    {
      start[ndims] = 0;
      count[ndims] = zaxisInqSize(zaxisID);
      ndims++;
    }
  if ( yid != UNDEFID )
    {
      start[ndims] = 0;
      cdf_inq_dimlen(fileID, yid, &size);
      /*      count[ndims] = gridInqYsize(gridID); */
      count[ndims] = size;
      ndims++;
    }
  if ( xid != UNDEFID )
    {
      start[ndims] = 0;
      cdf_inq_dimlen(fileID, xid, &size);
      /*      count[ndims] = gridInqXsize(gridID); */
      count[ndims] = size;
      ndims++;
    }

  if ( CDI_Debug )
    for (idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);

  if ( streamptr->ncmode == 1 )
    {
      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  dtype = vlistInqVarDatatype(vlistID, varID);

  if ( nmiss > 0 ) cdfDefVarMissval(streamID, varID, dtype, 1);

  /*  if ( dtype == DATATYPE_INT8 || dtype == DATATYPE_INT16 || dtype == DATATYPE_INT32 ) */
    {
      int nvals;
      int laddoffset, lscalefactor;
      double addoffset, scalefactor;
      double missval;

      addoffset    = vlistInqVarAddoffset(vlistID, varID);
      scalefactor  = vlistInqVarScalefactor(vlistID, varID);
      laddoffset   = IS_NOT_EQUAL(addoffset, 0);
      lscalefactor = IS_NOT_EQUAL(scalefactor, 1);

      missval = vlistInqVarMissval(vlistID, varID);

      nvals = gridInqSize(gridID)*zaxisInqSize(zaxisID);

      if ( laddoffset || lscalefactor )
        {
          mdata = (double *) malloc(nvals*sizeof(double));
          memcpy(mdata, data, nvals*sizeof(double));

          if ( nmiss > 0 )
            {
              for ( i = 0; i < nvals; i++ )
                {
                  if ( !DBL_IS_EQUAL(data[i], missval) )
                    {
                      if ( laddoffset )   mdata[i] -= addoffset;
                      if ( lscalefactor ) mdata[i] /= scalefactor;
                    }
                }
            }
          else
            {
              for ( i = 0; i < nvals; i++ )
                {
                  if ( laddoffset )   mdata[i] -= addoffset;
                  if ( lscalefactor ) mdata[i] /= scalefactor;
                }
            }
        }

      if ( dtype == DATATYPE_UINT8 || dtype == DATATYPE_INT8 ||
           dtype == DATATYPE_INT16 || dtype == DATATYPE_INT32 )
        {
          if ( mdata )
            {
              for ( i = 0; i < nvals; i++ ) mdata[i] = NINT(mdata[i]);
            }
          else
            {
              mdata = (double *) malloc(nvals*sizeof(double));
              for ( i = 0; i < nvals; i++ ) mdata[i] = NINT(data[i]);
            }
        }
    }

  if ( mdata )
    cdf_put_vara_double(fileID, ncvarid, start, count, mdata);
  else
    cdf_put_vara_double(fileID, ncvarid, start, count, data);

  if ( mdata ) free(mdata);
#endif
}


int cdfReadVarSliceDP(int streamID, int varID, int levelID, double *data, int *nmiss)
{
#if  defined  (HAVE_LIBNETCDF)
  int fileID;
  int gridID;
  int zaxisID;
  int xid = UNDEFID, yid = UNDEFID, zid = UNDEFID;
  int ncvarid;
  int tsID;
  int gridsize, xsize, ysize;
  size_t size;
  size_t start[4];
  size_t count[4];
  int ndims = 0;
  int idim;
  int timeID;
  int gridindex;
  int zaxisindex;
  int vlistID;
  int i, j;
  int dimorder[3];
  int ixyz;
  int swapxy = FALSE;
  int lvalidrange;
  double validrange[2];
  double missval;
  int laddoffset, lscalefactor;
  double addoffset, scalefactor;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( CDI_Debug )
    Message("streamID = %d  varID = %d  levelID = %d", streamID, varID, levelID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  tsID = streamptr->curTsID;
  if ( CDI_Debug )
    Message("tsID = %d", tsID);

  ncvarid = streamptr->vars[varID].ncvarid;

  gridID  = vlistInqVarGrid(vlistID, varID);
  zaxisID = vlistInqVarZaxis(vlistID, varID);
  timeID  = vlistInqVarTime(vlistID, varID);
  ixyz    = vlistInqVarXYZ(vlistID, varID);
  if ( ixyz == 0 ) ixyz = 321; // ZYX

  gridsize = gridInqSize(gridID);
  xsize = gridInqXsize(gridID);
  ysize = gridInqYsize(gridID);

  streamptr->numvals += gridsize;

  gridindex = vlistGridIndex(vlistID, gridID);
  if ( gridInqType(gridID) == GRID_TRAJECTORY )
    {
      cdfReadGridTraj(streamID, gridID);
    }
  else if ( gridInqType(gridID) == GRID_UNSTRUCTURED )
    {
      xid = streamptr->xdimID[gridindex];
    }
  else
    {
      xid = streamptr->xdimID[gridindex];
      yid = streamptr->ydimID[gridindex];
    }

  zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  zid = streamptr->zaxisID[zaxisindex];
  /*
  printf("2 %d %d %d %s\n", streamID, zaxisindex, streamptr->zaxisID[zaxisindex], vlistInqVarNamePtr(vlistID, varID));
  */
  dimorder[0] = ixyz/100;
  dimorder[1] = (ixyz-dimorder[0]*100)/10;
  dimorder[2] = (ixyz-dimorder[0]*100-dimorder[1]*10);

  if ( (dimorder[2] == 2 || dimorder[0] == 1) && xid != UNDEFID && yid != UNDEFID ) swapxy = TRUE;
  /*
  printf("swapxy %d\n", swapxy);
  printf("ixyz %d\n", ixyz);
  printf("dimorder: %d %d %d\n", dimorder[0], dimorder[1], dimorder[2]);
  */

  if ( timeID == TIME_VARIABLE )
    {
      start[ndims] = tsID;
      count[ndims] = 1;
      ndims++;
    }

  for ( int id = 0; id < 3; ++id )
    {
      if ( dimorder[id] == 3 && zid != UNDEFID )
        {
          start[ndims] = levelID;
          count[ndims] = 1;
          ndims++;
        }
      else if ( dimorder[id] == 2 && yid != UNDEFID )
        {
          start[ndims] = 0;
          cdf_inq_dimlen(fileID, yid, &size);
          count[ndims] = size;
          ndims++;
        }
      else if ( dimorder[id] == 1 && xid != UNDEFID )
        {
          start[ndims] = 0;
          cdf_inq_dimlen(fileID, xid, &size);
          count[ndims] = size;
          ndims++;
        }
    }

  if ( CDI_Debug )
    for (idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);

  cdf_get_vara_double(fileID, ncvarid, start, count, data);

  if ( swapxy )
    {
      double *tdata;
      tdata = (double *) malloc(gridsize*sizeof(double));
      memcpy(tdata, data, gridsize*sizeof(double));
      for ( j = 0; j < ysize; ++j )
        for ( i = 0; i < xsize; ++i )
          data[j*xsize+i] = tdata[i*ysize+j];
      free(tdata);
    }

  if ( vlistInqVarDatatype(vlistID, varID) == DATATYPE_UINT8 )
    {
      nc_type xtype;
      cdf_inq_vartype(fileID, ncvarid, &xtype);
      if ( xtype == NC_BYTE )
        {
          for ( i = 0; i < gridsize; i++ )
            if ( data[i] < 0 ) data[i] += 256;
        }
    }

  *nmiss = 0;
  if ( vlistInqVarMissvalUsed(vlistID, varID) == TRUE )
    {
      missval = vlistInqVarMissval(vlistID, varID);

      lvalidrange = vlistInqVarValidrange(vlistID, varID, validrange);
      // printf("readvarslice: validrange %d %g %g\n", lvalidrange, validrange[0], validrange[1]);
      if ( lvalidrange )
        for ( i = 0; i < gridsize; i++ )
          {
            if ( IS_NOT_EQUAL(validrange[0], VALIDMISS) && data[i] < validrange[0] ) data[i] = missval;
            if ( IS_NOT_EQUAL(validrange[1], VALIDMISS) && data[i] > validrange[1] ) data[i] = missval;
          }

      // printf("XXX %31.0f %31.0f %31.0f %31.0f\n", missval, (float)data[0]);
      for ( i = 0; i < gridsize; i++ )
        if ( DBL_IS_EQUAL(data[i], missval) ) *nmiss += 1;
    }

  addoffset    = vlistInqVarAddoffset(vlistID, varID);
  scalefactor  = vlistInqVarScalefactor(vlistID, varID);
  laddoffset   = IS_NOT_EQUAL(addoffset, 0);
  lscalefactor = IS_NOT_EQUAL(scalefactor, 1);

  if ( laddoffset || lscalefactor )
    {
      missval = vlistInqVarMissval(vlistID, varID);

      if ( *nmiss > 0 )
        {
          for ( i = 0; i < gridsize; i++ )
            {
              if ( !DBL_IS_EQUAL(data[i], missval) )
                {
                  if ( lscalefactor ) data[i] *= scalefactor;
                  if ( laddoffset )   data[i] += addoffset;
                }
            }
        }
      else
        {
          for ( i = 0; i < gridsize; i++ )
            {
              if ( lscalefactor ) data[i] *= scalefactor;
              if ( laddoffset )   data[i] += addoffset;
            }
        }
    }

#endif
  return (0);
}


int cdfWriteVarSliceDP(int streamID, int varID, int levelID, const double *data, int nmiss)
{
#if  defined  (HAVE_LIBNETCDF)
  int fileID;
  int gridID;
  int zaxisID;
  int xid = UNDEFID, yid = UNDEFID, zid = UNDEFID;
  int ncvarid;
  int ntsteps;
  size_t size, xsize, ysize;
  size_t start[4];
  size_t count[4];
  int nvals;
  int ndims = 0;
  int idim;
  int timeID;
  int gridindex, zaxisindex;
  int i;
  int dimorder[3];
  int ixyz;
  int swapxy = FALSE;
  int dtype;
  int vlistID;
  const double *pdata = data;
  double *mdata = NULL;
  double *sdata = NULL;
  stream_t *streamptr;
  extern int CDF_Debug;

  streamptr = stream_to_pointer(streamID);

  if ( CDI_Debug )
    Message("streamID = %d  varID = %d", streamID, varID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  ntsteps = streamptr->ntsteps;
  if ( CDI_Debug ) Message("ntsteps = %d", ntsteps);

  if ( vlistHasTime(vlistID) ) cdfDefTime(streamID);

  ncvarid = cdfDefVar(streamID, varID);

  gridID  = vlistInqVarGrid(vlistID, varID);
  zaxisID = vlistInqVarZaxis(vlistID, varID);
  timeID  = vlistInqVarTime(vlistID, varID);
  ixyz    = vlistInqVarXYZ(vlistID, varID);
  if ( ixyz == 0 ) ixyz = 321; // ZYX

  gridindex = vlistGridIndex(vlistID, gridID);
  if ( gridInqType(gridID) == GRID_TRAJECTORY )
    {
      cdfWriteGridTraj(streamID, gridID);
    }
  else
    {
      xid = streamptr->xdimID[gridindex];
      yid = streamptr->ydimID[gridindex];
    }

  zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  zid = streamptr->zaxisID[zaxisindex];

  dimorder[0] = ixyz/100;
  dimorder[1] = (ixyz-dimorder[0]*100)/10;
  dimorder[2] = (ixyz-dimorder[0]*100-dimorder[1]*10);

  if ( (dimorder[2] == 2 || dimorder[0] == 1) && xid != UNDEFID && yid != UNDEFID ) swapxy = TRUE;
  /*
  printf("swapxy %d\n", swapxy);
  printf("ixyz %d\n", ixyz);
  printf("dimorder: %d %d %d\n", dimorder[0], dimorder[1], dimorder[2]);
  */

  if ( timeID == TIME_VARIABLE )
    {
      start[ndims] = ntsteps - 1;
      count[ndims] = 1;
      ndims++;
    }

  for ( int id = 0; id < 3; ++id )
    {
      if ( dimorder[id] == 3 && zid != UNDEFID )
        {
          start[ndims] = levelID;
          count[ndims] = 1;
          ndims++;
        }
      else if ( dimorder[id] == 2 && yid != UNDEFID )
        {
          start[ndims] = 0;
          cdf_inq_dimlen(fileID, yid, &ysize);
          count[ndims] = ysize;
          ndims++;
        }
      else if ( dimorder[id] == 1 && xid != UNDEFID )
        {
          start[ndims] = 0;
          cdf_inq_dimlen(fileID, xid, &xsize);
          count[ndims] = xsize;
          ndims++;
        }
    }

  if ( CDI_Debug )
    for (idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);

  dtype = vlistInqVarDatatype(vlistID, varID);

  if ( nmiss > 0 ) cdfDefVarMissval(streamID, varID, dtype, 1);

  nvals = gridInqSize(gridID);

  /*  if ( dtype == DATATYPE_INT8 || dtype == DATATYPE_INT16 || dtype == DATATYPE_INT32 ) */
    {
      int laddoffset, lscalefactor;
      double addoffset, scalefactor;
      double missval;

      addoffset    = vlistInqVarAddoffset(vlistID, varID);
      scalefactor  = vlistInqVarScalefactor(vlistID, varID);
      laddoffset   = IS_NOT_EQUAL(addoffset, 0);
      lscalefactor = IS_NOT_EQUAL(scalefactor, 1);

      missval      = vlistInqVarMissval(vlistID, varID);

      if ( laddoffset || lscalefactor )
        {
          mdata = (double *) malloc(nvals*sizeof(double));
          memcpy(mdata, data, nvals*sizeof(double));

          if ( nmiss > 0 )
            {
              for ( i = 0; i < nvals; i++ )
                {
                  if ( !DBL_IS_EQUAL(mdata[i], missval) )
                    {
                      if ( laddoffset )   mdata[i] -= addoffset;
                      if ( lscalefactor ) mdata[i] /= scalefactor;
                    }
                }
            }
          else
            {
              for ( i = 0; i < nvals; i++ )
                {
                  if ( laddoffset )   mdata[i] -= addoffset;
                  if ( lscalefactor ) mdata[i] /= scalefactor;
                }
            }
        }

      if ( dtype == DATATYPE_UINT8 || dtype == DATATYPE_INT8 ||
           dtype == DATATYPE_INT16 || dtype == DATATYPE_INT32 )
        {
          if ( mdata )
            {
              for ( i = 0; i < nvals; i++ ) mdata[i] = NINT(mdata[i]);
            }
          else
            {
              mdata = (double *) malloc(nvals*sizeof(double));
              for ( i = 0; i < nvals; i++ ) mdata[i] = NINT(data[i]);
            }

          if ( dtype == DATATYPE_UINT8 )
            {
              nc_type xtype;
              cdf_inq_vartype(fileID, ncvarid, &xtype);
              if ( xtype == NC_BYTE )
                {
                  for ( i = 0; i < nvals; ++i )
                    if ( mdata[i] > 127 ) mdata[i] -= 256;
                }
            }
        }

      if ( CDF_Debug )
        {
          double fmin, fmax;
          fmin =  1.0e200;
          fmax = -1.0e200;
          for ( i = 0; i < nvals; ++i )
            {
              if ( !DBL_IS_EQUAL(data[i], missval) )
                {
                  if ( data[i] < fmin ) fmin = data[i];
                  if ( data[i] > fmax ) fmax = data[i];
                }
            }
          Message("nvals = %d, nmiss = %d, missval = %g, minval = %g, maxval = %g",
                  nvals, nmiss, missval, fmin, fmax);
        }

      if ( mdata ) pdata = mdata;
    }

  if ( swapxy )
    {
      sdata = (double *) malloc(nvals*sizeof(double));
      for ( int j = 0; j < (int)ysize; ++j )
        for ( int i = 0; i < (int)xsize; ++i )
           sdata[i*ysize+j] = pdata[j*xsize+i];
      pdata = sdata;
    }

  cdf_put_vara_double(fileID, ncvarid, start, count, pdata);

  if ( mdata ) free(mdata);
  if ( sdata ) free(sdata);
#endif
  return (0);
}

void cdfCreateRecords(int streamID, int tsID)
{
#if  defined  (HAVE_LIBNETCDF)
  int varID, levelID, recID, vrecID, zaxisID;
  int nvars, nlev, nrecs, nvrecs;
  record_t *records = NULL;
  int *recIDs = NULL;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID  = streamInqVlist(streamID);

  if ( tsID < 0 || (tsID >= streamptr->ntsteps && tsID > 0) ) return;

  if ( streamptr->tsteps[tsID].nallrecs > 0 ) return;

  if ( tsID == 0 )
    {
      nvars = vlistNvars(vlistID);
      nrecs = vlistNrecs(vlistID);

      streamptr->nrecs += nrecs;

      if ( nrecs > 0 ) records = (record_t *) malloc(nrecs*sizeof(record_t));
      streamptr->tsteps[tsID].records    = records;
      streamptr->tsteps[tsID].nrecs      = nrecs;
      streamptr->tsteps[tsID].nallrecs   = nrecs;
      streamptr->tsteps[tsID].recordSize = nrecs;
      streamptr->tsteps[tsID].curRecID   = UNDEFID;

      nvrecs = nrecs; /* use all records at first timestep */
      if ( nvrecs > 0 ) recIDs = (int *) malloc(nvrecs*sizeof(int));
      streamptr->tsteps[tsID].recIDs     = recIDs;
      for ( recID = 0; recID < nvrecs; recID++ )
        recIDs[recID] = recID;

      recID = 0;
      for ( varID = 0; varID < nvars; varID++ )
        {
          zaxisID = vlistInqVarZaxis(vlistID, varID);
          nlev    = zaxisInqSize(zaxisID);
          for ( levelID = 0; levelID < nlev; levelID++ )
            {
              recordInitEntry(&records[recID]);
              records[recID].varID   = varID;
              records[recID].levelID = levelID;
              recID++;
            }
        }
    }
  else if ( tsID == 1 )
    {
      nvars = vlistNvars(vlistID);
      nrecs = vlistNrecs(vlistID);

      nvrecs = 0;
      for ( varID = 0; varID < nvars; varID++ )
        {
          if ( vlistInqVarTime(vlistID, varID) == TIME_VARIABLE )
            {
              zaxisID = vlistInqVarZaxis(vlistID, varID);
              nvrecs += zaxisInqSize(zaxisID);
            }
        }

      streamptr->nrecs += nvrecs;

      records = (record_t *) malloc(nrecs*sizeof(record_t));
      streamptr->tsteps[tsID].records    = records;
      streamptr->tsteps[tsID].nrecs      = nvrecs;
      streamptr->tsteps[tsID].nallrecs   = nrecs;
      streamptr->tsteps[tsID].recordSize = nrecs;
      streamptr->tsteps[tsID].curRecID   = UNDEFID;

      memcpy(streamptr->tsteps[tsID].records,
             streamptr->tsteps[0].records,
             nrecs*sizeof(record_t));

      if ( nvrecs )
        {
          recIDs = (int *) malloc(nvrecs*sizeof(int));
          streamptr->tsteps[tsID].recIDs     = recIDs;
          vrecID = 0;
          for ( recID = 0; recID < nrecs; recID++ )
            {
              varID = records[recID].varID;
              if ( vlistInqVarTime(vlistID, varID) == TIME_VARIABLE )
                {
                  recIDs[vrecID++] = recID;
                }
            }
        }
    }
  else
    {
      nvars = vlistNvars(vlistID);
      nrecs = vlistNrecs(vlistID);

      nvrecs = streamptr->tsteps[1].nrecs;

      streamptr->nrecs += nvrecs;

      records = (record_t *) malloc(nrecs*sizeof(record_t));
      streamptr->tsteps[tsID].records    = records;
      streamptr->tsteps[tsID].nrecs      = nvrecs;
      streamptr->tsteps[tsID].nallrecs   = nrecs;
      streamptr->tsteps[tsID].recordSize = nrecs;
      streamptr->tsteps[tsID].curRecID   = UNDEFID;

      memcpy(streamptr->tsteps[tsID].records,
             streamptr->tsteps[0].records,
             nrecs*sizeof(record_t));

      recIDs = (int *) malloc(nvrecs*sizeof(int));
      streamptr->tsteps[tsID].recIDs     = recIDs;

      memcpy(streamptr->tsteps[tsID].recIDs,
             streamptr->tsteps[1].recIDs,
             nvrecs*sizeof(int));
    }
#endif
}

static
int cdfTimeDimID(int fileID, int ndims, int nvars)
{
  int dimid = 0;
#if  defined  (HAVE_LIBNETCDF)
  char dimname[80];
  char timeunits[256];
  char attname[256];
  char name[256];
  nc_type xtype;
  int nvdims, nvatts;
  int dimids[9];
  int varid, iatt;

  for ( dimid = 0; dimid < ndims; dimid++ )
    {
      cdf_inq_dimname(fileID, dimid, dimname);
      if ( memcmp(dimname, "time", 4) == 0 ) break;
    }

  if ( dimid == ndims ) dimid = UNDEFID;

  for ( varid = 0; varid < nvars; varid++ )
    {
      if ( dimid != UNDEFID ) break;

      cdf_inq_var(fileID, varid, name, &xtype, &nvdims, dimids, &nvatts);
      if ( nvdims == 1 )
        {
          for ( iatt = 0; iatt < nvatts; iatt++ )
            {
              cdf_inq_attname(fileID, varid, iatt, attname);
              if ( memcmp(attname, "units", 5) == 0 )
                {
                  cdfGetAttText(fileID, varid, "units", sizeof(timeunits), timeunits);
                  strtolower(timeunits);

                  if ( memcmp(timeunits, "sec",    3) == 0 ||
                       memcmp(timeunits, "minute", 6) == 0 ||
                       memcmp(timeunits, "hour",   4) == 0 ||
                       memcmp(timeunits, "day",    3) == 0 ||
                       memcmp(timeunits, "month",  5) == 0 )
                    {
                      dimid = dimids[0];
                      break;
                    }
                }
            }
        }
    }

#endif
  return (dimid);
}

static
void init_ncdims(long ndims, ncdim_t *ncdims)
{
  long ncdimid;

  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      ncdims[ncdimid].ncvarid      = UNDEFID;
      ncdims[ncdimid].dimtype      = UNDEFID;
      ncdims[ncdimid].len          = 0;
      ncdims[ncdimid].name[0]      = 0;
    }
}

static
void init_ncvars(long nvars, ncvar_t *ncvars)
{
  long ncvarid;

  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      ncvars[ncvarid].ignore          = FALSE;
      ncvars[ncvarid].isvar           = UNDEFID;
      ncvars[ncvarid].islon           = FALSE;
      ncvars[ncvarid].islat           = FALSE;
      ncvars[ncvarid].islev           = FALSE;
      ncvars[ncvarid].warn            = FALSE;
      ncvars[ncvarid].timeID          = TIME_CONSTANT;
      ncvars[ncvarid].param           = UNDEFID;
      ncvars[ncvarid].code            = UNDEFID;
      ncvars[ncvarid].tabnum          = 0;
      ncvars[ncvarid].calendar        = FALSE;
      ncvars[ncvarid].bounds          = UNDEFID;
      ncvars[ncvarid].gridID          = UNDEFID;
      ncvars[ncvarid].zaxisID         = UNDEFID;
      ncvars[ncvarid].gridtype        = UNDEFID;
      ncvars[ncvarid].zaxistype       = UNDEFID;
      ncvars[ncvarid].xdim            = UNDEFID;
      ncvars[ncvarid].ydim            = UNDEFID;
      ncvars[ncvarid].zdim            = UNDEFID;
      ncvars[ncvarid].xvarid          = UNDEFID;
      ncvars[ncvarid].yvarid          = UNDEFID;
      ncvars[ncvarid].zvarid          = UNDEFID;
      ncvars[ncvarid].tvarid          = UNDEFID;
      ncvars[ncvarid].ncoordvars      = 0;
      ncvars[ncvarid].coordvarids[0]  = UNDEFID;
      ncvars[ncvarid].coordvarids[1]  = UNDEFID;
      ncvars[ncvarid].coordvarids[2]  = UNDEFID;
      ncvars[ncvarid].coordvarids[3]  = UNDEFID;
      ncvars[ncvarid].cellarea        = UNDEFID;
      ncvars[ncvarid].tableID         = UNDEFID;
      ncvars[ncvarid].xtype           = 0;
      ncvars[ncvarid].ndims           = 0;
      ncvars[ncvarid].gmapid          = UNDEFID;
      ncvars[ncvarid].vlen            = 0;
      ncvars[ncvarid].vdata           = NULL;
      ncvars[ncvarid].truncation      = 0;
      ncvars[ncvarid].positive        = 0;
      ncvars[ncvarid].defmiss         = 0;
      ncvars[ncvarid].missval         = 0;
      ncvars[ncvarid].addoffset       = 0;
      ncvars[ncvarid].scalefactor     = 1;
      ncvars[ncvarid].name[0]         = 0;
      ncvars[ncvarid].longname[0]     = 0;
      ncvars[ncvarid].stdname[0]      = 0;
      ncvars[ncvarid].units[0]        = 0;
      ncvars[ncvarid].natts           = 0;
      ncvars[ncvarid].atts            = NULL;
      ncvars[ncvarid].deflate         = 0;
      ncvars[ncvarid].lunsigned       = 0;
      ncvars[ncvarid].lvalidrange     = 0;
      ncvars[ncvarid].validrange[0]   = VALIDMISS;
      ncvars[ncvarid].validrange[1]   = VALIDMISS;
    }
}

#if  defined  (HAVE_LIBNETCDF)
static
int isLonAxis(const char *units, const char *stdname)
{
  int status = FALSE;

  if ( memcmp(units, "degrees_east", 12) == 0 ||
       memcmp(units, "degree_east", 11)  == 0 ||
       memcmp(units, "degree_E", 8)      == 0 ||
       memcmp(units, "degrees_E", 9)     == 0 ||
       memcmp(units, "degreeE", 7)       == 0 ||
       memcmp(units, "degreesE", 8)      == 0 ||
       memcmp(stdname, "longitude", 9)   == 0 ||
       (memcmp(units, "degree", 6)            == 0 &&
        memcmp(stdname, "grid_longitude", 14) == 0) ||
       (memcmp(units, "radian", 6)            == 0 &&
        memcmp(stdname, "grid_longitude", 14) == 0) )
    {
      status = TRUE;
    }

  return (status);
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
int isLatAxis(const char *units, const char *stdname)
{
  int status = FALSE;

  if ( memcmp(units, "degrees_north", 13) == 0 ||
       memcmp(units, "degree_north", 12)  == 0 ||
       memcmp(units, "degree_N", 8)       == 0 ||
       memcmp(units, "degrees_N", 9)      == 0 ||
       memcmp(units, "degreeN", 7)        == 0 ||
       memcmp(units, "degreesN", 8)       == 0 ||
       memcmp(stdname, "latitude", 8)     == 0 ||
       (memcmp(units, "degree", 6)           == 0 &&
        memcmp(stdname, "grid_latitude", 13) == 0) ||
       (memcmp(units, "radian", 6)           == 0 &&
        memcmp(stdname, "grid_latitude", 13) == 0) )
    {
      status = TRUE;
    }

  return (status);
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
int isDBLAxis(const char *longname)
{
  int status = FALSE;

  if ( strcmp(longname, "depth below land")         == 0 ||
       strcmp(longname, "depth_below_land")         == 0 ||
       strcmp(longname, "levels below the surface") == 0 )
    {
      status = TRUE;
    }

  return (status);
}

static
int isDepthAxis(const char *stdname, const char *longname)
{
  int status = FALSE;

  if ( strcmp(stdname, "depth") == 0 ) status = TRUE;
 
  if ( status == FALSE )
    if ( strcmp(longname, "depth_below_sea") == 0 ||
         strcmp(longname, "depth below sea") == 0 )
      {
        status = TRUE;
      }

  return (status);
}

static
int isHeightAxis(const char *stdname, const char *longname)
{
  int status = FALSE;

  if ( strcmp(stdname, "height") == 0 ) status = TRUE;

  if ( status == FALSE )
    if ( strcmp(longname, "height") == 0 ||
         strcmp(longname, "height above the surface") == 0 )
      {
        status = TRUE;
      }

  return (status);
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
int unitsIsPressure(const char *units)
{
  int status = FALSE;

  if ( memcmp(units, "millibar", 8) == 0 ||
       memcmp(units, "mb", 2)       == 0 ||
       memcmp(units, "hectopas", 8) == 0 ||
       memcmp(units, "hPa", 3)      == 0 ||
       memcmp(units, "Pa", 2)       == 0 )
    {
      status = TRUE;
    }

  return (status);
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
int isGaussGrid(long ysize, double yinc, double *yvals)
{
  int lgauss = FALSE;
  long i;
  double *yv, *yw;

  if ( IS_EQUAL(yinc, 0) && ysize > 2 ) /* check if gaussian */
    {
      yv = (double *) malloc(ysize*sizeof(double));
      yw = (double *) malloc(ysize*sizeof(double));
      gaussaw(yv, yw, ysize);
      free(yw);
      for ( i = 0; i < ysize; i++ )
        yv[i] = asin(yv[i])/M_PI*180.0;

      for ( i = 0; i < ysize; i++ )
        if ( fabs(yv[i] - yvals[i]) >
             ((yv[0] - yv[1])/500) ) break;

      if ( i == ysize ) lgauss = TRUE;

      /* check S->N */
      if ( lgauss == FALSE )
        {
          for ( i = 0; i < ysize; i++ )
            if ( fabs(yv[i] - yvals[ysize-i-1]) >
                 ((yv[0] - yv[1])/500) ) break;

          if ( i == ysize ) lgauss = TRUE;
        }

      free(yv);
    }

  return (lgauss);
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void cdfSetVar(ncvar_t *ncvars, int ncvarid, int isvar)
{
  if ( isvar != TRUE && isvar != FALSE )
    Error("Internal problem! var %s undefined", ncvars[ncvarid].name);

  if ( ncvars[ncvarid].isvar != UNDEFID &&
       ncvars[ncvarid].isvar != isvar   &&
       ncvars[ncvarid].warn  == FALSE )
    {
      if ( ! ncvars[ncvarid].ignore )
        Warning("Inconsistent variable definition for %s!", ncvars[ncvarid].name);

      ncvars[ncvarid].warn = TRUE;
      isvar = FALSE;
    }

  ncvars[ncvarid].isvar = isvar;
}
#endif

static
void cdfSetDim(ncvar_t *ncvars, int ncvarid, int dimid, int dimtype)
{
  if ( ncvars[ncvarid].dimtype[dimid] != UNDEFID &&
       ncvars[ncvarid].dimtype[dimid] != dimtype )
    {
      Warning("Inconsistent dimension definition for %s! dimid = %d;  type = %d;  newtype = %d",
              ncvars[ncvarid].name, dimid, ncvars[ncvarid].dimtype[dimid], dimtype);
    }

  ncvars[ncvarid].dimtype[dimid] = dimtype;
}

static
void printNCvars(ncvar_t *ncvars, int nvars)
{
  char axis[6];
  int ncvarid, i;
  int ndim;
  int iaxis[] = {'t', 'z', 'y', 'x'};

  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      ndim = 0;
      if ( ncvars[ncvarid].isvar )
        {
          axis[ndim++] = 'v';
          for ( i = 0; i < ncvars[ncvarid].ndims; i++ )
            {/*
              if      ( ncvars[ncvarid].tvarid != -1 ) axis[ndim++] = iaxis[0];
              else if ( ncvars[ncvarid].zvarid != -1 ) axis[ndim++] = iaxis[1];
              else if ( ncvars[ncvarid].yvarid != -1 ) axis[ndim++] = iaxis[2];
              else if ( ncvars[ncvarid].xvarid != -1 ) axis[ndim++] = iaxis[3];
              else
             */
              axis[ndim++] = '?';
            }
        }
      else
        {
          axis[ndim++] = 'c';
          if ( ncvars[ncvarid].islev )
            axis[ndim++] = iaxis[1];
          else if ( ncvars[ncvarid].islat )
            axis[ndim++] = iaxis[2];
          else if ( ncvars[ncvarid].islon )
            axis[ndim++] = iaxis[3];
          else
            axis[ndim++] = '?';
        }

      axis[ndim++] = 0;

      printf("%3d %3d %3d  %-4s %s\n", ncvarid, ncvars[ncvarid].isvar, ndim-2, axis, ncvars[ncvarid].name);
    }
}


typedef struct
{
  int      ncvarid;
  char     name[256];
}
varinfo_t;


static
int cmpvarname(const void *s1, const void *s2)
{
  varinfo_t *x = (varinfo_t *) s1;
  varinfo_t *y = (varinfo_t *) s2;

  return (strcmp(x->name, y->name));
}

#if  defined  (HAVE_LIBNETCDF)
static
void cdfScanVarAttributes(int fileID, int nvars, ncvar_t *ncvars, ncdim_t *ncdims,
                          int timedimid, int modelID, int format)
{
  int ncvarid;
  int ncdimid;
  int nvdims, nvatts;
  int *dimidsp;
  nc_type xtype, atttype;
  size_t attlen;
  char name[256];
  char attname[256];
  const int attstringlen = 8192; char attstring[8192];
  int iatt;
  int i;
  int tablenum;

  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      dimidsp = ncvars[ncvarid].dimids;

      cdf_inq_var(fileID, ncvarid, name, &xtype, &nvdims, dimidsp, &nvatts);
      strcpy(ncvars[ncvarid].name, name);

      for ( ncdimid = 0; ncdimid < nvdims; ncdimid++ )
        ncvars[ncvarid].dimtype[ncdimid] = -1;

      ncvars[ncvarid].xtype = xtype;
      ncvars[ncvarid].ndims = nvdims;

#if  defined  (HAVE_NETCDF4)
      if ( format == NC_FORMAT_NETCDF4_CLASSIC || format == NC_FORMAT_NETCDF4 )
        {
          int shuffle, deflate, deflate_level;
          nc_inq_var_deflate(fileID, ncvarid, &shuffle, &deflate, &deflate_level);
          if ( deflate > 0 )
            ncvars[ncvarid].deflate = 1;
        }
#endif

      if ( nvdims > 0 )
        if ( timedimid == dimidsp[0] )
          {
            ncvars[ncvarid].timeID = TIME_VARIABLE;
            cdfSetDim(ncvars, ncvarid, 0, T_AXIS);
          }

      for ( iatt = 0; iatt < nvatts; iatt++ )
        {
          cdf_inq_attname(fileID, ncvarid, iatt, attname);
          cdf_inq_atttype(fileID, ncvarid, attname, &atttype);
          cdf_inq_attlen(fileID, ncvarid, attname, &attlen);

          if ( strcmp(attname, "long_name") == 0 && atttype == NC_CHAR )
            {
              cdfGetAttText(fileID, ncvarid, attname, MAXNAMELEN, ncvars[ncvarid].longname);
            }
          else if ( strcmp(attname, "standard_name") == 0 && atttype == NC_CHAR )
            {
              cdfGetAttText(fileID, ncvarid, attname, MAXNAMELEN, ncvars[ncvarid].stdname);
            }
          else if ( strcmp(attname, "units") == 0 && atttype == NC_CHAR )
            {
              cdfGetAttText(fileID, ncvarid, attname, MAXNAMELEN, ncvars[ncvarid].units);
            }
          else if ( strcmp(attname, "calendar") == 0 )
            {
              ncvars[ncvarid].calendar = TRUE;
            }
          else if ( strcmp(attname, "param") == 0 && atttype == NC_CHAR )
            {
	      char paramstr[32];
	      int pnum = 0, pcat = 255, pdis = 255;
              cdfGetAttText(fileID, ncvarid, attname, sizeof(paramstr), paramstr);
	      sscanf(paramstr, "%d.%d.%d", &pnum, &pcat, &pdis);
	      ncvars[ncvarid].param = cdiEncodeParam(pnum, pcat, pdis);
              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "code") == 0 && atttype != NC_CHAR )
            {
              cdfGetAttInt(fileID, ncvarid, attname, 1, &ncvars[ncvarid].code);
              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "table") == 0 && atttype != NC_CHAR )
            {
              cdfGetAttInt(fileID, ncvarid, attname, 1, &tablenum);
              if ( tablenum > 0 )
                {
                  ncvars[ncvarid].tabnum = tablenum;
                  ncvars[ncvarid].tableID = tableInq(modelID, tablenum, NULL);
                  if ( ncvars[ncvarid].tableID == CDI_UNDEFID )
                    ncvars[ncvarid].tableID = tableDef(modelID, tablenum, NULL);
                }
              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "trunc_type") == 0 && atttype == NC_CHAR )
            {
              cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);
              if ( memcmp(attstring, "Triangular", attlen) == 0 )
                ncvars[ncvarid].gridtype = GRID_SPECTRAL;
            }
          else if ( strcmp(attname, "grid_type") == 0 && atttype == NC_CHAR )
            {
              cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);
              strtolower(attstring);

              if      ( strcmp(attstring, "gaussian reduced") == 0 )
                ncvars[ncvarid].gridtype = GRID_GAUSSIAN_REDUCED;
              else if ( strcmp(attstring, "gaussian") == 0 )
                ncvars[ncvarid].gridtype = GRID_GAUSSIAN;
              else if ( strncmp(attstring, "spectral", 8) == 0 )
                ncvars[ncvarid].gridtype = GRID_SPECTRAL;
              else if ( strncmp(attstring, "fourier", 7) == 0 )
                ncvars[ncvarid].gridtype = GRID_FOURIER;
              else if ( strcmp(attstring, "trajectory") == 0 )
                ncvars[ncvarid].gridtype = GRID_TRAJECTORY;
              else if ( strcmp(attstring, "generic") == 0 )
                ncvars[ncvarid].gridtype = GRID_GENERIC;
              else if ( strcmp(attstring, "cell") == 0 )
                ncvars[ncvarid].gridtype = GRID_UNSTRUCTURED;
              else if ( strcmp(attstring, "unstructured") == 0 )
                ncvars[ncvarid].gridtype = GRID_UNSTRUCTURED;
              else if ( strcmp(attstring, "curvilinear") == 0 )
                ncvars[ncvarid].gridtype = GRID_CURVILINEAR;
              else if ( strcmp(attstring, "sinusoidal") == 0 )
                ;
              else if ( strcmp(attstring, "laea") == 0 )
                ;
              else if ( strcmp(attstring, "lcc2") == 0 )
                ;
              else if ( strcmp(attstring, "linear") == 0 ) // ignore grid type linear
                ;
              else
                {
                  static int warn = TRUE;
                  if ( warn )
                    {
                      warn = FALSE;
                      Warning("netCDF attribute grid_type='%s' unsupported!", attstring);
                    }
                }

              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "level_type") == 0 && atttype == NC_CHAR )
            {
              cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);
              strtolower(attstring);

              if      ( strcmp(attstring, "toa") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_TOA;
              else if ( strcmp(attstring, "seabottom") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_SEA_BOTTOM;
              else if ( strcmp(attstring, "atmosphere") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_ATMOSPHERE;
              else
                { 
                  static int warn = TRUE;
                  if ( warn )
                    {
                      warn = FALSE;
                      Warning("netCDF attribute level_type='%s' unsupported!", attstring);
                    }
                }

              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "trunc_count") == 0 && atttype != NC_CHAR )
            {
              cdfGetAttInt(fileID, ncvarid, attname, 1, &ncvars[ncvarid].truncation);
            }
          else if ( strcmp(attname, "truncation") == 0 && atttype != NC_CHAR )
            {
              cdfGetAttInt(fileID, ncvarid, attname, 1, &ncvars[ncvarid].truncation);
            }
          else if ( strcmp(attname, "add_offset") == 0 && atttype != NC_CHAR )
            {
	      cdfGetAttDouble(fileID, ncvarid, attname, 1, &ncvars[ncvarid].addoffset);
	      /*
		if ( atttype != NC_BYTE && atttype != NC_SHORT && atttype != NC_INT )
		if ( ncvars[ncvarid].addoffset != 0 )
		Warning("attribute add_offset not supported for atttype %d", atttype);
	      */
	      /* (also used for lon/lat) cdfSetVar(ncvars, ncvarid, TRUE); */
            }
          else if ( strcmp(attname, "scale_factor") == 0 && atttype != NC_CHAR )
            {
	      cdfGetAttDouble(fileID, ncvarid, attname, 1, &ncvars[ncvarid].scalefactor);
	      /*
		if ( atttype != NC_BYTE && atttype != NC_SHORT && atttype != NC_INT )
		if ( ncvars[ncvarid].scalefactor != 1 )
		Warning("attribute scale_factor not supported for atttype %d", atttype);
	      */
	      /* (also used for lon/lat) cdfSetVar(ncvars, ncvarid, TRUE); */
            }
          else if ( strcmp(attname, "bounds") == 0 && atttype == NC_CHAR )
            {
              int status, ncboundsid;

              cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);

              status = nc_inq_varid(fileID, attstring, &ncboundsid);

              if ( status == NC_NOERR )
                {
                  ncvars[ncvarid].bounds = ncboundsid;
                  cdfSetVar(ncvars, ncvars[ncvarid].bounds, FALSE);
                  cdfSetVar(ncvars, ncvarid, FALSE);
                }
              else
                Warning("%s - %s", nc_strerror(status), attstring);
            }
          else if ( strcmp(attname, "cell_measures") == 0 && atttype == NC_CHAR )
            {
              char *pstring, *cell_measures = NULL, *cell_var = NULL;

              cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);
              pstring = attstring;

              while ( isspace((int) *pstring) ) pstring++;
              cell_measures = pstring;
              while ( isalnum((int) *pstring) ) pstring++;
              *pstring++ = 0;
              while ( isspace((int) *pstring) ) pstring++;
              cell_var = pstring;
              while ( ! isspace((int) *pstring) && *pstring != 0 ) pstring++;
              *pstring++ = 0;
              /*
              printf("cell_measures >%s<\n", cell_measures);
              printf("cell_var >%s<\n", cell_var);
              */
              if ( memcmp(cell_measures, "area", 4) == 0 )
                {
                  int status;
                  int nc_cell_id;

                  status = nc_inq_varid(fileID, cell_var, &nc_cell_id);
                  if ( status == NC_NOERR )
                    {
                      ncvars[ncvarid].cellarea = nc_cell_id;
                      /* ncvars[nc_cell_id].isvar = UNDEFID; */
                      cdfSetVar(ncvars, nc_cell_id, FALSE);
                    }
                  else
                    Warning("%s - %s", nc_strerror(status), cell_var);
                }
              else
                {
                  Warning("%s has unexpected contents: %s", attname, cell_measures);
                }
              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          /*
          else if ( strcmp(attname, "coordinates") == 0 )
            {
              char *pstring, *xvarname = NULL, *yvarname = NULL;

              cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);
              pstring = attstring;

              while ( isspace((int) *pstring) ) pstring++;
              xvarname = pstring;
              while ( isgraph((int) *pstring) ) pstring++;
              *pstring++ = 0;
              while ( isspace((int) *pstring) ) pstring++;
              yvarname = pstring;
              while ( isgraph((int) *pstring) ) pstring++;
              *pstring++ = 0;

              cdf_inq_varid(fileID, xvarname, &ncvars[ncvarid].xvarid);
              cdf_inq_varid(fileID, yvarname, &ncvars[ncvarid].yvarid);

              cdfSetVar(ncvars, ncvars[ncvarid].xvarid, FALSE);
              cdfSetVar(ncvars, ncvars[ncvarid].yvarid, FALSE);
              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          */
          else if ( (strcmp(attname, "associate")  == 0 || strcmp(attname, "coordinates") == 0) &&
		    atttype == NC_CHAR )
            {
              int status;
              char *pstring, *varname = NULL;
              int lstop = FALSE;
              int dimvarid;
              extern int cdiIgnoreAttCoordinates;

              cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);
              pstring = attstring;

              for ( i = 0; i < 4; i++ )
                {
                  while ( isspace((int) *pstring) ) pstring++;
                  if ( *pstring == 0 ) break;
                  varname = pstring;
                  while ( !isspace((int) *pstring) && *pstring != 0 ) pstring++;
                  if ( *pstring == 0 ) lstop = TRUE;
                  *pstring++ = 0;

                  status = nc_inq_varid(fileID, varname, &dimvarid);
                  if ( status == NC_NOERR )
                    {
                      cdfSetVar(ncvars, dimvarid, FALSE);
                      if ( cdiIgnoreAttCoordinates == FALSE )
                        {
                          ncvars[ncvarid].coordvarids[i] = dimvarid;
                          ncvars[ncvarid].ncoordvars++;
                        }
                    }
                  else
                    Warning("%s - %s", nc_strerror(status), varname);

                  if ( lstop ) break;
                }

              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "grid_mapping") == 0 && atttype == NC_CHAR )
            {
              int status;
              int nc_gmap_id;

              cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);

              status = nc_inq_varid(fileID, attstring, &nc_gmap_id);
              if ( status == NC_NOERR )
                {
                  ncvars[ncvarid].gmapid = nc_gmap_id;
                  cdfSetVar(ncvars, ncvars[ncvarid].gmapid, FALSE);
                }
              else
                Warning("%s - %s", nc_strerror(status), attstring);

              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "positive") == 0 && atttype == NC_CHAR )
            {
              cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);
              strtolower(attstring);

              if    ( memcmp(attstring, "down", 4) == 0 ) ncvars[ncvarid].positive = -1;
              else if ( memcmp(attstring, "up", 2) == 0 ) ncvars[ncvarid].positive = 1;

              if ( ncvars[ncvarid].ndims == 1 )
                {
                  cdfSetVar(ncvars, ncvarid, FALSE);
                  cdfSetDim(ncvars, ncvarid, 0, Z_AXIS);
                  ncdims[ncvars[ncvarid].dimids[0]].dimtype = Z_AXIS;
                }
            }
          else if ( (strcmp(attname, "_FillValue") == 0 || strcmp(attname, "missing_value") == 0) &&
		    atttype != NC_CHAR )
            {
	      cdfGetAttDouble(fileID, ncvarid, attname, 1, &ncvars[ncvarid].missval);
	      ncvars[ncvarid].defmiss = TRUE;
	      /* cdfSetVar(ncvars, ncvarid, TRUE); */
            }
          else if ( strcmp(attname, "valid_range") == 0 && attlen == 2 )
            {
              if ( ncvars[ncvarid].lvalidrange == FALSE )
                {
                  cdfGetAttDouble(fileID, ncvarid, attname, 2, ncvars[ncvarid].validrange);
                  ncvars[ncvarid].lvalidrange = TRUE;
                  if ( ((int)ncvars[ncvarid].validrange[0]) == 0 && ((int)ncvars[ncvarid].validrange[1]) == 255 )
                    ncvars[ncvarid].lunsigned = TRUE;
                  /* cdfSetVar(ncvars, ncvarid, TRUE); */
                }
            }
          else if ( strcmp(attname, "valid_min") == 0 && attlen == 1 )
            {
              cdfGetAttDouble(fileID, ncvarid, attname, 1, &(ncvars[ncvarid].validrange)[0]);
              ncvars[ncvarid].lvalidrange = TRUE;
            }
          else if ( strcmp(attname, "valid_max") == 0 && attlen == 1 )
            {
              cdfGetAttDouble(fileID, ncvarid, attname, 1, &(ncvars[ncvarid].validrange)[1]);
              ncvars[ncvarid].lvalidrange = TRUE;
            }
          else if ( strcmp(attname, "_Unsigned") == 0 && atttype == NC_CHAR )
            {
              cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);
              strtolower(attstring);

              if ( memcmp(attstring, "true", 4) == 0 )
                {
                  ncvars[ncvarid].lunsigned = TRUE;
                  /*
                  ncvars[ncvarid].lvalidrange = TRUE;
                  ncvars[ncvarid].validrange[0] = 0;
                  ncvars[ncvarid].validrange[1] = 255;
                  */
                }
	      /* cdfSetVar(ncvars, ncvarid, TRUE); */
            }
          else if ( strcmp(attname, "cdi") == 0 && atttype == NC_CHAR )
            {
	      cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);
	      strtolower(attstring);

	      if ( memcmp(attstring, "ignore", 6) == 0 )
		{
		  ncvars[ncvarid].ignore = TRUE;
		  cdfSetVar(ncvars, ncvarid, FALSE);
		}
            }
          else if ( strcmp(attname, "axis") == 0 && atttype == NC_CHAR )
            {
              cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);
	      attlen = strlen(attstring);

	      if ( (int) attlen > nvdims )
		{
		  if ( nvdims > 0 )
		    Warning("Unexpected axis attribute length for %s, ignored!", name);
		}
	      else
		{
		  strtolower(attstring);
		  for ( i = 0; i < (int)attlen; ++i )
		    {
		      if ( attstring[i] != '-' && attstring[i] != 't' && attstring[i] != 'z' &&
			   attstring[i] != 'y' && attstring[i] != 'x' )
			{
			  Warning("Unexpected character in axis attribute for %s, ignored!", name);
			  break;
			}
		    }

		  if ( i == (int) attlen && (int) attlen == nvdims)
		    {
		      while ( attlen-- )
			{
			  if ( tolower((int) attstring[attlen]) == 't' )
			    {
			      if ( attlen != 0 ) Warning("axis attribute 't' not on first position");
			      cdfSetDim(ncvars, ncvarid, attlen, T_AXIS);
			    }
			  else if ( tolower((int) attstring[attlen]) == 'z' )
			    {
			      ncvars[ncvarid].zdim = dimidsp[attlen];
			      cdfSetDim(ncvars, ncvarid, attlen, Z_AXIS);

			      if ( ncvars[ncvarid].ndims == 1 )
				{
				  cdfSetVar(ncvars, ncvarid, FALSE);
				  ncdims[ncvars[ncvarid].dimids[0]].dimtype = Z_AXIS;
				}
			    }
			  else if ( tolower((int) attstring[attlen]) == 'y' )
			    {
			      ncvars[ncvarid].ydim = dimidsp[attlen];
			      cdfSetDim(ncvars, ncvarid, attlen, Y_AXIS);
			    }
			  else if ( tolower((int) attstring[attlen]) == 'x' )
			    {
			      ncvars[ncvarid].xdim = dimidsp[attlen];
			      cdfSetDim(ncvars, ncvarid, attlen, X_AXIS);
			    }
			}
		    }
		}
	    }
	  else
	    {
	      if ( ncvars[ncvarid].natts == 0 )
		ncvars[ncvarid].atts = (int *) malloc(nvatts*sizeof(int));

	      ncvars[ncvarid].atts[ncvars[ncvarid].natts++] = iatt;
	      /*
	      int attrint;
	      double attrflt;
	      nc_type attrtype;
	      cdf_inq_attlen(fileID, ncvarid, attname, &attlen);
	      cdf_inq_atttype(fileID, ncvarid, attname, &attrtype);
	      if ( attlen == 1 && (attrtype == NC_INT || attrtype == NC_SHORT) )
		{
		  cdfGetAttInt(fileID, ncvarid, attname, 1, &attrint);
		  printf("int: %s.%s = %d\n", ncvars[ncvarid].name, attname, attrint);
		}
	      else if ( attlen == 1 && (attrtype == NC_FLOAT || attrtype == NC_DOUBLE) )
		{
		  cdfGetAttDouble(fileID, ncvarid, attname, 1, &attrflt);
		  printf("flt: %s.%s = %g\n", ncvars[ncvarid].name, attname, attrflt);
		}
	      else if ( attrtype == NC_CHAR )
		{
		  cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);
		  attstring[attlen] = 0;
		  printf("txt: %s.%s = %s\n", ncvars[ncvarid].name, attname, attstring);
		}
	      else
		printf("att: %s.%s = unknown\n", ncvars[ncvarid].name, attname);
	      */
	    }
	}
    }
}
#endif

static
void setDimType(int nvars, ncvar_t *ncvars, ncdim_t *ncdims)
{
  int ndims;
  int ncvarid, ncdimid;
  int i;

  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar == TRUE )
	{
	  int lxdim = 0, lydim = 0, lzdim = 0, ltdim = 0;
	  ndims = ncvars[ncvarid].ndims;
	  for ( i = 0; i < ndims; i++ )
	    {
	      ncdimid = ncvars[ncvarid].dimids[i];
	      if      ( ncdims[ncdimid].dimtype == X_AXIS ) cdfSetDim(ncvars, ncvarid, i, X_AXIS);
	      else if ( ncdims[ncdimid].dimtype == Y_AXIS ) cdfSetDim(ncvars, ncvarid, i, Y_AXIS);
	      else if ( ncdims[ncdimid].dimtype == Z_AXIS ) cdfSetDim(ncvars, ncvarid, i, Z_AXIS);
	      else if ( ncdims[ncdimid].dimtype == T_AXIS ) cdfSetDim(ncvars, ncvarid, i, T_AXIS);
	    }

	  if ( CDI_Debug )
	    {
	      Message("var %d %s", ncvarid, ncvars[ncvarid].name);
	      for ( i = 0; i < ndims; i++ )
		printf("  dim %d type %d  ", i, ncvars[ncvarid].dimtype[i]);
	      printf("\n");
	    }

	  for ( i = 0; i < ndims; i++ )
	    {
	      if      ( ncvars[ncvarid].dimtype[i] == X_AXIS ) lxdim = TRUE;
	      else if ( ncvars[ncvarid].dimtype[i] == Y_AXIS ) lydim = TRUE;
	      else if ( ncvars[ncvarid].dimtype[i] == Z_AXIS ) lzdim = TRUE;
	      else if ( ncvars[ncvarid].dimtype[i] == T_AXIS ) ltdim = TRUE;
	    }

	  for ( i = ndims-1; i >= 0; i-- )
	    {
	      if ( ncvars[ncvarid].dimtype[i] == -1 )
		{
		  /*
		  printf("undef dim: %d %d %d %d %d\n", i, lxdim, lydim, lzdim, ltdim);
		  */
		  if ( lxdim == FALSE )
		    {
		      cdfSetDim(ncvars, ncvarid, i, X_AXIS);
		      lxdim = TRUE;
		    }
		  else if ( lydim == FALSE && ncvars[ncvarid].gridtype != GRID_UNSTRUCTURED )
		    {
		      cdfSetDim(ncvars, ncvarid, i, Y_AXIS);
		      lydim = TRUE;
		    }
		  else if ( lzdim == FALSE )
		    {
		      cdfSetDim(ncvars, ncvarid, i, Z_AXIS);
		      lzdim = TRUE;
		    }
		}
	    }
	}
    }
}

#if  defined  (HAVE_LIBNETCDF)
/* verify coordinate vars - first scan (dimname == varname) */
static
void verify_coordinate_vars_1(int ndims, ncdim_t *ncdims, ncvar_t *ncvars, int timedimid)
{
  int ncdimid, ncvarid;

  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      ncvarid = ncdims[ncdimid].ncvarid;
      if ( ncvarid != -1 )
	{
	  if ( ncvars[ncvarid].dimids[0] == timedimid )
	    {
	      ncdims[ncdimid].dimtype = T_AXIS;
	      continue;
	    }

	  if ( ncvars[ncvarid].longname[0] != 0 && ncvars[ncvarid].longname[1] != 0 )
	    {
	      if ( memcmp(ncvars[ncvarid].longname+1, "ongitude", 8) == 0 )
		{
		  ncvars[ncvarid].islon = TRUE;
		  cdfSetVar(ncvars, ncvarid, FALSE);
		  cdfSetDim(ncvars, ncvarid, 0, X_AXIS);
		  ncdims[ncdimid].dimtype = X_AXIS;
		  continue;
		}
	      else if ( memcmp(ncvars[ncvarid].longname+1, "atitude", 7) == 0 )
		{
		  ncvars[ncvarid].islat = TRUE;
		  cdfSetVar(ncvars, ncvarid, FALSE);
		  cdfSetDim(ncvars, ncvarid, 0, Y_AXIS);
		  ncdims[ncdimid].dimtype = Y_AXIS;
		  continue;
		}
	    }

	  if ( ncvars[ncvarid].units[0] != 0 )
	    {
	      if ( isLonAxis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].islon = TRUE;
		  cdfSetVar(ncvars, ncvarid, FALSE);
		  cdfSetDim(ncvars, ncvarid, 0, X_AXIS);
		  ncdims[ncdimid].dimtype = X_AXIS;
		}
	      else if ( isLatAxis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].islat = TRUE;
		  cdfSetVar(ncvars, ncvarid, FALSE);
		  cdfSetDim(ncvars, ncvarid, 0, Y_AXIS);
		  ncdims[ncdimid].dimtype = Y_AXIS;
		}
	      else if ( unitsIsPressure(ncvars[ncvarid].units) )
		{
		  ncvars[ncvarid].zaxistype = ZAXIS_PRESSURE;
		}
	      else if ( strcmp(ncvars[ncvarid].units, "level") == 0 || strcmp(ncvars[ncvarid].units, "1") == 0 )
		{
		  if      ( strcmp(ncvars[ncvarid].longname, "hybrid level at layer midpoints") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID;
		  else if ( memcmp(ncvars[ncvarid].longname, "hybrid level at midpoints", 25) == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID;
		  else if ( strcmp(ncvars[ncvarid].longname, "hybrid level at layer interfaces") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID_HALF;
		  else if ( memcmp(ncvars[ncvarid].longname, "hybrid level at interfaces", 26) == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID_HALF;
		  else if ( strcmp(ncvars[ncvarid].units, "level") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_GENERIC;
		}
	      else if ( strcmp(ncvars[ncvarid].units, "cm")  == 0 )
		{
		  if ( isDBLAxis(ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_DEPTH_BELOW_LAND;
		}
	      else if ( strcmp(ncvars[ncvarid].units, "m")   == 0 )
		{
		  if ( isDepthAxis(ncvars[ncvarid].stdname, ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_DEPTH_BELOW_SEA;
		  else if ( isHeightAxis(ncvars[ncvarid].stdname, ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_HEIGHT;
		}
	    }

	  if ( ncvars[ncvarid].zaxistype != UNDEFID )
	    {
	      cdfSetVar(ncvars, ncvarid, FALSE);
	      cdfSetDim(ncvars, ncvarid, 0, Z_AXIS);
	      ncdims[ncdimid].dimtype = Z_AXIS;
	    }
	}
    }
}
#endif

#if  defined  (HAVE_LIBNETCDF)
/* verify coordinate vars - second scan (all other variables) */
static
void verify_coordinate_vars_2(int nvars, ncvar_t *ncvars)
{
  int ncvarid;

  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar == 0 )
	{
	  /* not needed anymore for rotated grids */
	  if ( ncvars[ncvarid].longname[0] != 0 && ncvars[ncvarid].longname[1] != 0 )
	    {
	      if ( memcmp(ncvars[ncvarid].longname+1, "ongitude", 8) == 0 )
		{
		  ncvars[ncvarid].islon = TRUE;
		  continue;
		}
	      else if ( memcmp(ncvars[ncvarid].longname+1, "atitude", 7) == 0 )
		{
		  ncvars[ncvarid].islat = TRUE;
		  continue;
		}
	    }

	  if ( ncvars[ncvarid].units[0] != 0 )
	    {
	      if ( isLonAxis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].islon = TRUE;
		  continue;
		}
	      else if ( isLatAxis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].islat = TRUE;
		  continue;
		}
	      else if ( unitsIsPressure(ncvars[ncvarid].units) )
		{
		  ncvars[ncvarid].zaxistype = ZAXIS_PRESSURE;
		  continue;
		}
	      else if ( strcmp(ncvars[ncvarid].units, "level") == 0 || strcmp(ncvars[ncvarid].units, "1") == 0 )
		{
		  if      ( strcmp(ncvars[ncvarid].longname, "hybrid level at layer midpoints") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID;
		  else if ( memcmp(ncvars[ncvarid].longname, "hybrid level at midpoints", 25) == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID;
		  else if ( strcmp(ncvars[ncvarid].longname, "hybrid level at layer interfaces") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID_HALF;
		  else if ( memcmp(ncvars[ncvarid].longname, "hybrid level at interfaces", 26) == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID_HALF;
		  else if ( strcmp(ncvars[ncvarid].units, "level") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_GENERIC;
		  continue;
		}
	      else if ( strcmp(ncvars[ncvarid].units, "cm")  == 0 )
		{
		  if ( isDBLAxis(ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_DEPTH_BELOW_LAND;
		  continue;
		}
	      else if ( strcmp(ncvars[ncvarid].units, "m")   == 0 )
		{
		  if ( isDepthAxis(ncvars[ncvarid].stdname, ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_DEPTH_BELOW_SEA;
		  else if ( isHeightAxis(ncvars[ncvarid].stdname, ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_HEIGHT;
		  continue;
		}
	    }
	}
    }
}
#endif

#if  defined  (HAVE_LIBNETCDF)
static
void copy_numeric_projatts(int gridID, int ncvarID, int ncfileID)
{
  int iatt, nvatts;
  size_t attlen;
  char attname[256];
  nc_type xtype;

  cdf_inq_varnatts(ncfileID, ncvarID, &nvatts);

  for ( iatt = 0; iatt < nvatts; iatt++ )
    {
      cdf_inq_attname(ncfileID, ncvarID, iatt, attname);
      cdf_inq_atttype(ncfileID, ncvarID, attname, &xtype);
      cdf_inq_attlen(ncfileID, ncvarID, attname, &attlen);

      printf("%s %d\n", attname, (int)attlen);
    }

}
#endif

#if  defined  (HAVE_LIBNETCDF)
/* define all input grids */
static
void define_all_grids(stream_t *streamptr, int fileID, int vlistID, ncdim_t *ncdims, int nvars, ncvar_t *ncvars, int timedimid)
{
  int ncvarid, ncvarid2;
  int ndims;
  int nbdims;
  int i;
  int nvatts;
  size_t nvertex;
  grid_t grid;
  grid_t proj;
  int gridindex;
  size_t size = 0, xsize, ysize, np;
  char name[256];
  int iatt;
  int ltwarn = TRUE;
  size_t attlen;
  char attname[256];
  const int attstringlen = 8192; char attstring[8192];
  double datt;

  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar && ncvars[ncvarid].gridID == UNDEFID )
	{
	  int xdimids[2] = {-1,-1}, ydimids[2] = {-1,-1};
	  int xdimid = -1, ydimid = -1;
	  int xvarid = -1, yvarid = -1;
	  int islon = 0, islat = 0;
	  int nxdims = 0, nydims = 0;
	  double xinc = 0, yinc = 0;

	  xsize = 0;
	  ysize = 0;
          np    = 0;

	  ndims = ncvars[ncvarid].ndims;
	  for ( i = 0; i < ndims; i++ )
	    {
	      if ( ncvars[ncvarid].dimtype[i] == X_AXIS && nxdims < 2 )
		{
		  xdimids[nxdims] = ncvars[ncvarid].dimids[i];
		  nxdims++;
		}
	      else if ( ncvars[ncvarid].dimtype[i] == Y_AXIS && nydims < 2 )
		{
		  ydimids[nydims] = ncvars[ncvarid].dimids[i];
		  nydims++;
		}
	    }

	  if ( nxdims == 2 )
	    {
	      xdimid = xdimids[1];
	      ydimid = xdimids[0];
	    }
	  else if ( nydims == 2 )
	    {
	      xdimid = ydimids[1];
	      ydimid = ydimids[0];
	    }
	  else
	    {
	      xdimid = xdimids[0];
	      ydimid = ydimids[0];
	    }

	  if ( ncvars[ncvarid].xvarid != UNDEFID )
	    xvarid = ncvars[ncvarid].xvarid;
	  else if ( xdimid != UNDEFID )
	    xvarid = ncdims[xdimid].ncvarid;

	  if ( ncvars[ncvarid].yvarid != UNDEFID )
	    yvarid = ncvars[ncvarid].yvarid;
	  else if ( ydimid != UNDEFID )
	    yvarid = ncdims[ydimid].ncvarid;

	  /*
	  if ( xdimid != UNDEFID )
	    xvarid = ncdims[xdimid].ncvarid;
	  if ( xvarid == UNDEFID && ncvars[ncvarid].xvarid != UNDEFID )
	    xvarid = ncvars[ncvarid].xvarid;

	  if ( ydimid != UNDEFID )
	    yvarid = ncdims[ydimid].ncvarid;
	  if ( yvarid == UNDEFID && ncvars[ncvarid].yvarid != UNDEFID )
	    yvarid = ncvars[ncvarid].yvarid;
	  */

	  if ( xdimid != UNDEFID ) xsize = ncdims[xdimid].len;
	  if ( ydimid != UNDEFID ) ysize = ncdims[ydimid].len;

	  if ( ydimid == UNDEFID && yvarid != UNDEFID )
	    {
	      if ( ncvars[yvarid].ndims == 1 )
		{
		  ydimid = ncvars[yvarid].dimids[0];
		  ysize  = ncdims[ydimid].len;
		}
	    }

	  if ( ncvars[ncvarid].gridtype == UNDEFID || ncvars[ncvarid].gridtype == GRID_GENERIC )
	    if ( ydimid == xdimid ) ncvars[ncvarid].gridtype = GRID_UNSTRUCTURED;

	  grid_init(&grid);
	  grid_init(&proj);

	  grid.prec  = DATATYPE_FLT64;
	  grid.trunc = ncvars[ncvarid].truncation;

	  if ( ncvars[ncvarid].gridtype == GRID_TRAJECTORY )
	    {
	      if ( ncvars[ncvarid].xvarid == UNDEFID )
		Error("Longitude coordinate undefined for %s!", name);
	      if ( ncvars[ncvarid].yvarid == UNDEFID )
		Error("Latitude coordinate undefined for %s!", name);
	    }
	  else
	    {
	      size_t start[3], count[3];
	      int ltgrid = FALSE;

	      if ( xvarid != UNDEFID && yvarid != UNDEFID )
		{
		  if ( ncvars[xvarid].ndims != ncvars[yvarid].ndims )
		    {
		      Warning("Inconsistent grid structure for variable %s!",
			      ncvars[ncvarid].name);
		      ncvars[ncvarid].xvarid = UNDEFID;
		      ncvars[ncvarid].yvarid = UNDEFID;
		      xvarid = UNDEFID;
		      yvarid = UNDEFID;
		    }

		  if ( ncvars[xvarid].ndims > 2 || ncvars[yvarid].ndims > 2 )
		    {
		      if ( ncvars[xvarid].ndims == 3 && ncvars[xvarid].dimids[0] == timedimid &&
			   ncvars[yvarid].ndims == 3 && ncvars[yvarid].dimids[0] == timedimid )
			{
			  if ( ltwarn )
			    Warning("Time varying grids unsupported, using grid at time step 1!");
			  ltgrid = TRUE;
			  ltwarn = FALSE;
			  start[0] = start[1] = start[2] = 0;
			  count[0] = 1; count[1] = ysize; count[2] = xsize;
			}
		      else
			{
			  Warning("Unsupported grid structure for variable %s (grid dims > 2)!", ncvars[ncvarid].name);
			  ncvars[ncvarid].xvarid = UNDEFID;
			  ncvars[ncvarid].yvarid = UNDEFID;
			  xvarid = UNDEFID;
			  yvarid = UNDEFID;
			}
		    }
		}

	      if ( xvarid != UNDEFID )
		{
		  islon = ncvars[xvarid].islon;
		  ndims = ncvars[xvarid].ndims;
		  if ( ndims == 2 || ndims == 3 )
		    {
		      ncvars[ncvarid].gridtype = GRID_CURVILINEAR;
		      size = xsize*ysize;
		      /* Check size of 2 dimensional coordinate variables */
		      {
			int dimid;
			size_t dimsize1, dimsize2;
			dimid = ncvars[xvarid].dimids[ndims-2];
			dimsize1 = ncdims[dimid].len;
			dimid = ncvars[xvarid].dimids[ndims-1];
			dimsize2 = ncdims[dimid].len;
			if ( dimsize1*dimsize2 != size )
			  {
			    Warning("Unsupported array structure, skipped variable %s!", ncvars[ncvarid].name);
			    ncvars[ncvarid].isvar = -1;
			    continue;
			  }
		      }
		    }
		  else
		    {
		      size = xsize;
		      /* Check size of 1 dimensional coordinate variables */
		      {
			int dimid;
			size_t dimsize;
			dimid = ncvars[xvarid].dimids[0];
			dimsize = ncdims[dimid].len;
			if ( dimsize != size )
			  {
			    Warning("Unsupported array structure, skipped variable %s!", ncvars[ncvarid].name);
			    ncvars[ncvarid].isvar = -1;
			    continue;
			  }
		      }
		    }

		  if ( ncvars[xvarid].xtype == NC_FLOAT ) grid.prec = DATATYPE_FLT32;
		  grid.xvals = (double *) malloc(size*sizeof(double));

		  if ( ltgrid )
		    cdf_get_vara_double(fileID, xvarid, start, count, grid.xvals);
		  else
		    cdf_get_var_double(fileID, xvarid, grid.xvals);

                  scale_add(size, grid.xvals, ncvars[xvarid].addoffset, ncvars[xvarid].scalefactor);

		  strcpy(grid.xname, ncvars[xvarid].name);
		  strcpy(grid.xlongname, ncvars[xvarid].longname);
		  strcpy(grid.xunits, ncvars[xvarid].units);
		  /* don't change the name !!! */
		  /*
		  if ( (len = strlen(grid.xname)) > 2 )
		    if ( grid.xname[len-2] == '_' && isdigit((int) grid.xname[len-1]) )
		      grid.xname[len-2] = 0;
		  */
		  if ( islon && xsize > 1 )
		    {
		      xinc = fabs(grid.xvals[0] - grid.xvals[1]);
		      for ( i = 2; i < (int) xsize; i++ )
			if ( (fabs(grid.xvals[i-1] - grid.xvals[i]) - xinc) > (xinc/1000) ) break;

		      if ( i < (int) xsize ) xinc = 0;
		    }
		}

	      if ( yvarid != UNDEFID )
		{
		  islat = ncvars[yvarid].islat;
		  ndims = ncvars[yvarid].ndims;
		  if ( ndims == 2 || ndims == 3 )
		    {
		      ncvars[ncvarid].gridtype = GRID_CURVILINEAR;
		      size = xsize*ysize;
		      /* Check size of 2 dimensional coordinate variables */
		      {
			int dimid;
			size_t dimsize1, dimsize2;
			dimid = ncvars[yvarid].dimids[ndims-2];
			dimsize1 = ncdims[dimid].len;
			dimid = ncvars[yvarid].dimids[ndims-1];
			dimsize2 = ncdims[dimid].len;
			if ( dimsize1*dimsize2 != size )
			  {
			    Warning("Unsupported array structure, skipped variable %s!", ncvars[ncvarid].name);
			    ncvars[ncvarid].isvar = -1;
			    continue;
			  }
		      }
		    }
		  else
		    {
		      if ( (int) ysize == 0 ) size = xsize;
		      else                    size = ysize;

		      /* Check size of 1 dimensional coordinate variables */
		      {
			int dimid;
			size_t dimsize;
			dimid = ncvars[yvarid].dimids[0];
			dimsize = ncdims[dimid].len;
			if ( dimsize != size )
			  {
			    Warning("Unsupported array structure, skipped variable %s!", ncvars[ncvarid].name);
			    ncvars[ncvarid].isvar = -1;
			    continue;
			  }
		      }
		    }

		  if ( ncvars[yvarid].xtype == NC_FLOAT ) grid.prec = DATATYPE_FLT32;
		  grid.yvals = (double *) malloc(size*sizeof(double));

		  if ( ltgrid )
		    cdf_get_vara_double(fileID, yvarid, start, count, grid.yvals);
		  else
		    cdf_get_var_double(fileID, yvarid, grid.yvals);

                  scale_add(size, grid.yvals, ncvars[xvarid].addoffset, ncvars[xvarid].scalefactor);

		  strcpy(grid.yname, ncvars[yvarid].name);
		  strcpy(grid.ylongname, ncvars[yvarid].longname);
		  strcpy(grid.yunits, ncvars[yvarid].units);
		  /* don't change the name !!! */
		  /*
		  if ( (len = strlen(grid.yname)) > 2 )
		    if ( grid.yname[len-2] == '_' && isdigit((int) grid.yname[len-1]) )
		      grid.yname[len-2] = 0;
		  */
		  if ( islon && (int) ysize > 1 )
		    {
		      yinc = fabs(grid.yvals[0] - grid.yvals[1]);
		      for ( i = 2; i < (int) ysize; i++ )
			if ( (fabs(grid.yvals[i-1] - grid.yvals[i]) - yinc) > (yinc/1000) ) break;

		      if ( i < (int) ysize ) yinc = 0;
		    }
		}

	      if      ( (int) ysize == 0 ) size = xsize;
	      else if ( (int) xsize == 0 ) size = ysize;
	      else if ( ncvars[ncvarid].gridtype == GRID_UNSTRUCTURED ) size = xsize; 
	      else                         size = xsize*ysize;
	    }

	  if ( ncvars[ncvarid].gridtype == UNDEFID ||
	       ncvars[ncvarid].gridtype == GRID_GENERIC )
	    {
	      if ( islat && islon )
		{
		  if ( isGaussGrid(ysize, yinc, grid.yvals) )
                    {
                      ncvars[ncvarid].gridtype = GRID_GAUSSIAN;
                      np = ysize/2;
                    }
                  else
		    ncvars[ncvarid].gridtype = GRID_LONLAT;
		}
	      else if ( islat && !islon && xsize == 0 )
		{
		  if ( isGaussGrid(ysize, yinc, grid.yvals) )
                    {
                      ncvars[ncvarid].gridtype = GRID_GAUSSIAN;
                      np = ysize/2;
                    }
                  else
		    ncvars[ncvarid].gridtype = GRID_LONLAT;
		}
	      else if ( islon && !islat && ysize == 0 )
		{
		  ncvars[ncvarid].gridtype = GRID_LONLAT;
		}
	      else
		ncvars[ncvarid].gridtype = GRID_GENERIC;
	    }

	  switch (ncvars[ncvarid].gridtype)
	    {
	    case GRID_GENERIC:
	    case GRID_LONLAT:
	    case GRID_GAUSSIAN:
	    case GRID_UNSTRUCTURED:
	    case GRID_CURVILINEAR:
	      {
		grid.size  = size;
		grid.xsize = xsize;
		grid.ysize = ysize;
                grid.np    = np;
		if ( xvarid != UNDEFID )
		  {
		    grid.xdef  = 1;
		    if ( ncvars[xvarid].bounds != UNDEFID )
		      {
			nbdims = ncvars[ncvars[xvarid].bounds].ndims;
			if ( nbdims == 2 || nbdims == 3 )
			  {
			    nvertex = ncdims[ncvars[ncvars[xvarid].bounds].dimids[nbdims-1]].len;
			    grid.nvertex = (int) nvertex;
			    grid.xbounds = (double *) malloc(nvertex*size*sizeof(double));
			    cdf_get_var_double(fileID, ncvars[xvarid].bounds, grid.xbounds);
			  }
		      }
		  }
		if ( yvarid != UNDEFID )
		  {
		    grid.ydef  = 1;
		    if ( ncvars[yvarid].bounds != UNDEFID )
		      {
			nbdims = ncvars[ncvars[yvarid].bounds].ndims;
			if ( nbdims == 2 || nbdims == 3 )
			  {
			    nvertex = ncdims[ncvars[ncvars[yvarid].bounds].dimids[nbdims-1]].len;
			    /*
			    if ( nvertex != grid.nvertex )
			      Warning("nvertex problem! nvertex x %d, nvertex y %d",
				      grid.nvertex, (int) nvertex);
			    */
			    grid.ybounds = (double *) malloc(nvertex*size*sizeof(double));
			    cdf_get_var_double(fileID, ncvars[yvarid].bounds, grid.ybounds);
			  }
		      }
		  }

		if ( ncvars[ncvarid].cellarea != UNDEFID )
		  {
		    grid.area = (double *) malloc(size*sizeof(double));
		    cdf_get_var_double(fileID, ncvars[ncvarid].cellarea, grid.area);
		  }

		break;
	      }
	    case GRID_SPECTRAL:
	      {
		grid.size = size;
		grid.lcomplex = 1;
		break;
	      }
	    case GRID_FOURIER:
	      {
		grid.size = size;
		break;
	      }
	    case GRID_TRAJECTORY:
	      {
		grid.size = 1;
		break;
	      }
	    }

	  grid.type = ncvars[ncvarid].gridtype;

	  if ( grid.size == 0 )
	    {
	      if ( (ncvars[ncvarid].ndims == 1 && ncvars[ncvarid].dimtype[0] == T_AXIS) ||
		   (ncvars[ncvarid].ndims == 2 && ncvars[ncvarid].dimtype[0] == T_AXIS &&
		    ncvars[ncvarid].dimtype[1] == Z_AXIS) )
		{
		  grid.type  = GRID_GENERIC;
		  grid.size  = 1;
		  grid.xsize = 0;
		  grid.ysize = 0;
		}
	      else
		{
		  Warning("Variable %s has unsupported grid, skipped!", ncvars[ncvarid].name);
		  ncvars[ncvarid].isvar = -1;
		  continue;
		}
	    }

	  if ( ncvars[ncvarid].gmapid >= 0 && ncvars[ncvarid].gridtype != GRID_CURVILINEAR )
	    {
	      cdf_inq_varnatts(fileID, ncvars[ncvarid].gmapid, &nvatts);

	      for ( iatt = 0; iatt < nvatts; iatt++ )
		{
		  cdf_inq_attname(fileID, ncvars[ncvarid].gmapid, iatt, attname);
		  cdf_inq_attlen(fileID, ncvars[ncvarid].gmapid, attname, &attlen);

		  if ( strcmp(attname, "grid_mapping_name") == 0 )
		    {
		      cdfGetAttText(fileID, ncvars[ncvarid].gmapid, attname, attstringlen-1, attstring);
		      strtolower(attstring);

		      if ( strcmp(attstring, "rotated_latitude_longitude") == 0 )
			grid.isRotated = TRUE;
		      else if ( strcmp(attstring, "sinusoidal") == 0 )
			grid.type = GRID_SINUSOIDAL;
		      else if ( strcmp(attstring, "lambert_azimuthal_equal_area") == 0 )
			grid.type = GRID_LAEA;
		      else if ( strcmp(attstring, "lambert_conformal_conic") == 0 )
			grid.type = GRID_LCC2;
		      else if ( strcmp(attstring, "lambert_cylindrical_equal_area") == 0 )
			{
			  proj.type = GRID_PROJECTION;
			  proj.name = strdup(attstring);
			}
		    }
		  else if ( strcmp(attname, "earth_radius") == 0 )
		    {
		      cdfGetAttDouble(fileID, ncvars[ncvarid].gmapid, attname, 1, &datt);
		      grid.laea_a = datt;
		      grid.lcc2_a = datt;
		    }
		  else if ( strcmp(attname, "longitude_of_projection_origin") == 0 )
		    {
		      cdfGetAttDouble(fileID, ncvars[ncvarid].gmapid, attname, 1, &grid.laea_lon_0);
		    }
		  else if ( strcmp(attname, "longitude_of_central_meridian") == 0 )
		    {
		      cdfGetAttDouble(fileID, ncvars[ncvarid].gmapid, attname, 1, &grid.lcc2_lon_0);
		    }
		  else if ( strcmp(attname, "latitude_of_projection_origin") == 0 )
		    {
		      cdfGetAttDouble(fileID, ncvars[ncvarid].gmapid, attname, 1, &datt);
		      grid.laea_lat_0 = datt;
		      grid.lcc2_lat_0 = datt;
		    }
		  else if ( strcmp(attname, "standard_parallel") == 0 )
		    {
		      if ( attlen == 1 )
			{
			  cdfGetAttDouble(fileID, ncvars[ncvarid].gmapid, attname, 1, &datt);
			  grid.lcc2_lat_1 = datt;
			  grid.lcc2_lat_2 = datt;
			}
		      else
			{
			  double datt2[2];
			  cdfGetAttDouble(fileID, ncvars[ncvarid].gmapid, attname, 2, datt2);
			  grid.lcc2_lat_1 = datt2[0];
			  grid.lcc2_lat_2 = datt2[1];
			}
		    }
		  else if ( strcmp(attname, "grid_north_pole_latitude") == 0 )
		    {
		      cdfGetAttDouble(fileID, ncvars[ncvarid].gmapid, attname, 1, &grid.ypole);
		    }
		  else if ( strcmp(attname, "grid_north_pole_longitude") == 0 )
		    {
		      cdfGetAttDouble(fileID, ncvars[ncvarid].gmapid, attname, 1, &grid.xpole);
		    }
		  else if ( strcmp(attname, "north_pole_grid_longitude") == 0 )
		    {
		      cdfGetAttDouble(fileID, ncvars[ncvarid].gmapid, attname, 1, &grid.angle);
		    }
		}
	    }

#if defined (PROJECTION_TEST)
	  if ( proj.type == GRID_PROJECTION )
	    {
	      if ( grid.type == GRID_GENERIC )
		{
		  grid.type = GRID_CURVILINEAR;
		}

	      if ( grid.type == GRID_CURVILINEAR )
		{
		  proj.size  = grid.size;
		  proj.xsize = grid.xsize;
                  proj.ysize = grid.ysize;
		}

	      //  grid.proj = gridGenerate(proj);
	    }
#endif

	  if ( CDI_Debug )
	    {
	      Message("grid: type = %d, size = %d, nx = %d, ny %d",
		      grid.type, grid.size, grid.xsize, grid.ysize);
	      Message("proj: type = %d, size = %d, nx = %d, ny %d",
		      proj.type, proj.size, proj.xsize, proj.ysize);
	    }

#if defined (PROJECTION_TEST)
	  if ( proj.type == GRID_PROJECTION )
	    {
	      ncvars[ncvarid].gridID = varDefGrid(vlistID, proj, 1);
	      copy_numeric_projatts(ncvars[ncvarid].gridID, ncvars[ncvarid].gmapid, fileID);
	    }
	  else
#endif
	    ncvars[ncvarid].gridID = varDefGrid(vlistID, grid, 1);

	  gridindex = vlistGridIndex(vlistID, ncvars[ncvarid].gridID);
	  streamptr->xdimID[gridindex] = xdimid;
	  streamptr->ydimID[gridindex] = ydimid;

	  grid_free(&grid);
	  grid_free(&proj);

	  if ( CDI_Debug )
	    Message("gridID %d %d %s", ncvars[ncvarid].gridID, ncvarid, ncvars[ncvarid].name);

	  for ( ncvarid2 = ncvarid+1; ncvarid2 < nvars; ncvarid2++ )
	    if ( ncvars[ncvarid2].isvar == TRUE && ncvars[ncvarid2].gridID == UNDEFID )
	      {
		int xdimid2 = -1, ydimid2 = -1;
		ndims = ncvars[ncvarid2].ndims;
		for ( i = 0; i < ndims; i++ )
		  {
		    if ( ncvars[ncvarid2].dimtype[i] == X_AXIS )
		      xdimid2 = ncvars[ncvarid2].dimids[i];
		    else if ( ncvars[ncvarid2].dimtype[i] == Y_AXIS )
		      ydimid2 = ncvars[ncvarid2].dimids[i];
		  }

		if ( xdimid == xdimid2 &&
		     (ydimid == ydimid2 || (xdimid == ydimid && ydimid2 == UNDEFID)) )
		  {
		    int same_grid = TRUE;

		    if ( xvarid != -1 && ncvars[ncvarid2].xvarid != UNDEFID &&
			 xvarid != ncvars[ncvarid2].xvarid ) same_grid = FALSE;

		    if ( yvarid != -1 && ncvars[ncvarid2].yvarid != UNDEFID &&
			 yvarid != ncvars[ncvarid2].yvarid ) same_grid = FALSE;

		    if ( same_grid )
		      {
			if ( CDI_Debug )
			  Message("gridID %d %d %s",
				  ncvars[ncvarid].gridID, ncvarid2, ncvars[ncvarid2].name);
			ncvars[ncvarid2].gridID = ncvars[ncvarid].gridID;
		      }
		  }
	      }
	}
    }
}
#endif

#if  defined  (HAVE_LIBNETCDF)
/* define all input zaxes */
static
void define_all_zaxes(stream_t *streamptr, int fileID, int vlistID, ncdim_t *ncdims, int nvars, ncvar_t *ncvars,
		      size_t vctsize, double *vct)
{
  int ncvarid, ncvarid2;
  int i, ilev, ndims;
  int zaxisindex;
  int zprec;
  int nbdims, nvertex, nlevel;
  char *pname, *plongname, *punits;

  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar == TRUE && ncvars[ncvarid].zaxisID == UNDEFID )
	{
	  int with_bounds = FALSE;
	  int zdimid = UNDEFID;
	  int zvarid = UNDEFID;
	  int zsize = 1;
	  double *zvar = NULL;
	  double *lbounds = NULL;
	  double *ubounds = NULL;
	  int zaxisType;

	  ndims = ncvars[ncvarid].ndims;
	  for ( i = 0; i < ndims; i++ )
	    {
	      if ( ncvars[ncvarid].dimtype[i] == Z_AXIS )
		zdimid = ncvars[ncvarid].dimids[i];
	    }

	  if ( zdimid != UNDEFID )
	    {
	      zvarid = ncdims[zdimid].ncvarid;
	      zsize  = ncdims[zdimid].len;
	    }

	  if ( CDI_Debug ) Message("nlevs = %d", zsize);

	  zvar = (double *) malloc(zsize*sizeof(double));

	  zaxisType = UNDEFID;

	  if ( zvarid != UNDEFID ) zaxisType = ncvars[zvarid].zaxistype;

	  if ( zaxisType == UNDEFID )  zaxisType = ZAXIS_GENERIC;

	  zprec = DATATYPE_FLT64;

	  if ( zvarid != UNDEFID )
	    {
	      pname     = ncvars[zvarid].name;
	      plongname = ncvars[zvarid].longname;
	      punits    = ncvars[zvarid].units;
	      if ( ncvars[zvarid].xtype == NC_FLOAT ) zprec = DATATYPE_FLT32;
	      /* don't change the name !!! */
	      /*
	      if ( (len = strlen(pname)) > 2 )
		if ( pname[len-2] == '_' && isdigit((int) pname[len-1]) )
		  pname[len-2] = 0;
	      */
	      cdf_get_var_double(fileID, zvarid, zvar);

	      if ( ncvars[zvarid].bounds != UNDEFID )
		{
		  nbdims = ncvars[ncvars[zvarid].bounds].ndims;
		  if ( nbdims == 2 )
		    {
		      nlevel  = ncdims[ncvars[ncvars[zvarid].bounds].dimids[0]].len;
		      nvertex = ncdims[ncvars[ncvars[zvarid].bounds].dimids[1]].len;
		      if ( nlevel == zsize && nvertex == 2 )
			{
			  double *zbounds;
			  with_bounds = TRUE;
			  zbounds = (double *) malloc(2*nlevel*sizeof(double));
			  lbounds = (double *) malloc(nlevel*sizeof(double));
			  ubounds = (double *) malloc(nlevel*sizeof(double));
			  cdf_get_var_double(fileID, ncvars[zvarid].bounds, zbounds);
			  for ( i = 0; i < nlevel; ++i )
			    {
			      lbounds[i] = zbounds[i*2];
			      ubounds[i] = zbounds[i*2+1];
			    }
			  free(zbounds);
			}
		    }
		}
	    }
	  else
	    {
	      pname     = NULL;
	      plongname = NULL;
	      punits    = NULL;

	      if ( zsize == 1 )
		{
                  if ( ncvars[ncvarid].zaxistype != UNDEFID )
                    zaxisType = ncvars[ncvarid].zaxistype;
                  else
                    zaxisType = ZAXIS_SURFACE;

		  zvar[0] = 0;
		  /*
		  if ( zdimid == UNDEFID )
		    zvar[0] = 9999;
		  else
		    zvar[0] = 0;
		  */
		}
	      else
		{
		  for ( ilev = 0; ilev < (int)zsize; ilev++ ) zvar[ilev] = ilev + 1;
		}
	    }

      	  ncvars[ncvarid].zaxisID = varDefZaxis(vlistID, zaxisType, (int) zsize, zvar, with_bounds, lbounds, ubounds,
						vctsize, vct, pname, plongname, punits, zprec, 1, 0);
	  free(zvar);
	  free(lbounds);
	  free(ubounds);

	  zaxisindex = vlistZaxisIndex(vlistID, ncvars[ncvarid].zaxisID);
	  streamptr->zaxisID[zaxisindex]  = zdimid;

	  if ( CDI_Debug )
	    Message("zaxisID %d %d %s", ncvars[ncvarid].zaxisID, ncvarid, ncvars[ncvarid].name);

	  for ( ncvarid2 = ncvarid+1; ncvarid2 < nvars; ncvarid2++ )
	    if ( ncvars[ncvarid2].isvar == TRUE && ncvars[ncvarid2].zaxisID == UNDEFID )
	      {
		int zdimid2 = -1;
		ndims = ncvars[ncvarid2].ndims;
		for ( i = 0; i < ndims; i++ )
		  {
		    if ( ncvars[ncvarid2].dimtype[i] == Z_AXIS )
		      zdimid2 = ncvars[ncvarid2].dimids[i];
		  }
		if ( zdimid == zdimid2 )
		  {
		    if ( CDI_Debug )
		      Message("zaxisID %d %d %s",
			      ncvars[ncvarid].zaxisID, ncvarid2, ncvars[ncvarid2].name);
		    ncvars[ncvarid2].zaxisID = ncvars[ncvarid].zaxisID;
		  }
	      }
	}
    }
}
#endif


#if  defined  (HAVE_LIBNETCDF)
/* define all input data variables */
static
void define_all_vars(int fileID, int streamID, int vlistID, int instID, int modelID, int tableID, int *varids, ncdim_t *ncdims, int nvars, ncvar_t *ncvars)
{
  int varID1, varID, ncvarid;
  int code;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  if ( streamptr->sortname )
    {
      int index;
      varinfo_t **varInfo;
      varInfo    = (varinfo_t **) malloc(nvars*sizeof(varinfo_t *));
      varInfo[0] = (varinfo_t *)  malloc(nvars*sizeof(varinfo_t));

      for ( index = 1; index < nvars; index++ )
	varInfo[index] = varInfo[0] + index;

      for ( varID = 0; varID < nvars; varID++ )
	{
	  ncvarid = varids[varID];
	  varInfo[varID]->ncvarid = ncvarid;
	  strcpy(varInfo[varID]->name, ncvars[ncvarid].name);
	}
      qsort(varInfo[0], nvars, sizeof(varinfo_t), cmpvarname);
      for ( varID = 0; varID < nvars; varID++ )
	{
	  varids[varID] = varInfo[varID]->ncvarid;
	}
      free(varInfo[0]);
      free(varInfo);
    }

  for ( varID1 = 0; varID1 < nvars; varID1++ )
    {
      int gridID, zaxisID;

      ncvarid = varids[varID1];
      gridID  = ncvars[ncvarid].gridID;
      zaxisID = ncvars[ncvarid].zaxisID;

      varID = streamNewVar(streamID, gridID, zaxisID);
      varID = vlistDefVar(vlistID, gridID, zaxisID, ncvars[ncvarid].timeID);

#if  defined  (HAVE_NETCDF4)
      if ( ncvars[ncvarid].deflate )
	vlistDefVarCompType(vlistID, varID, COMPRESS_ZIP);
#endif

      streamptr->vars[varID1].level   = NULL;
      streamptr->vars[varID1].defmiss = 0;
      streamptr->vars[varID1].nlevs   = zaxisInqSize(ncvars[ncvarid].zaxisID);
      streamptr->vars[varID1].ncvarid = ncvarid;

      vlistDefVarName(vlistID, varID, ncvars[ncvarid].name);
      if ( ncvars[ncvarid].param != UNDEFID ) vlistDefVarParam(vlistID, varID, ncvars[ncvarid].param);
      if ( ncvars[ncvarid].code != UNDEFID )  vlistDefVarCode(vlistID, varID, ncvars[ncvarid].code);
      if ( ncvars[ncvarid].code != UNDEFID )
	{
	  int param;
	  param = cdiEncodeParam(ncvars[ncvarid].code, ncvars[ncvarid].tabnum, 255);
	  vlistDefVarParam(vlistID, varID, param);
	}
      if ( ncvars[ncvarid].longname[0] )      vlistDefVarLongname(vlistID, varID, ncvars[ncvarid].longname);
      if ( ncvars[ncvarid].stdname[0] )       vlistDefVarStdname(vlistID, varID, ncvars[ncvarid].stdname);
      if ( ncvars[ncvarid].units[0] )         vlistDefVarUnits(vlistID, varID, ncvars[ncvarid].units);

      if ( ncvars[ncvarid].lvalidrange )
        vlistDefVarValidrange(vlistID, varID, ncvars[ncvarid].validrange);

      if ( IS_NOT_EQUAL(ncvars[ncvarid].addoffset, 0) )
	vlistDefVarAddoffset(vlistID, varID, ncvars[ncvarid].addoffset);
      if ( IS_NOT_EQUAL(ncvars[ncvarid].scalefactor, 1) )
	vlistDefVarScalefactor(vlistID, varID, ncvars[ncvarid].scalefactor);

      vlistDefVarDatatype(vlistID, varID, cdfInqDatatype(ncvars[ncvarid].xtype, ncvars[ncvarid].lunsigned));

      vlistDefVarInstitut(vlistID, varID, instID);
      vlistDefVarModel(vlistID, varID, modelID);
      if ( ncvars[ncvarid].tableID != UNDEFID )
	vlistDefVarTable(vlistID, varID, ncvars[ncvarid].tableID);

      if ( ncvars[ncvarid].defmiss == TRUE ) vlistDefVarMissval(vlistID, varID, ncvars[ncvarid].missval);

      if ( CDI_Debug )
	Message("varID = %d  gridID = %d  zaxisID = %d", varID,
		vlistInqVarGrid(vlistID, varID), vlistInqVarZaxis(vlistID, varID));

      int gridindex = vlistGridIndex(vlistID, gridID);
      int xdimid = streamptr->xdimID[gridindex];
      int ydimid = streamptr->ydimID[gridindex];

      int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
      int zdimid = streamptr->zaxisID[zaxisindex];

      int ndims = ncvars[ncvarid].ndims;
      int iodim = 0;
      int ixyz = 0;
      int ipow10[4] = {1, 10, 100, 1000};

      if ( ncvars[ncvarid].timeID == TIME_VARIABLE ) iodim++;

      if ( gridInqType(gridID) == GRID_UNSTRUCTURED && ndims-iodim <= 2 && ydimid == xdimid )
        {
          if ( xdimid == ncvars[ncvarid].dimids[ndims-1] )
            {
              ixyz = 321;
            }
          else
            {
              ixyz = 213;
            }
        }
      else
        {
          for ( int idim = iodim; idim < ndims; idim++ )
            {
              if      ( xdimid == ncvars[ncvarid].dimids[idim] )
                ixyz += 1*ipow10[ndims-idim-1];
              else if ( ydimid == ncvars[ncvarid].dimids[idim] )
                ixyz += 2*ipow10[ndims-idim-1];
              else if ( zdimid == ncvars[ncvarid].dimids[idim] )
                ixyz += 3*ipow10[ndims-idim-1];
            }
        }

      vlistDefVarXYZ(vlistID, varID, ixyz);
      /*
      printf("ixyz %d\n", ixyz);
      printf("ndims %d\n", ncvars[ncvarid].ndims);
      for ( int i = 0; i < ncvars[ncvarid].ndims; ++i )
        printf("dimids: %d %d\n", i, ncvars[ncvarid].dimids[i]);
      printf("xdimid, ydimid %d %d\n", xdimid, ydimid);
      */
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      ncvarid = varids[varID];

      if ( ncvars[ncvarid].natts )
	{
	  int nvatts;
	  int attnum;
	  int iatt;
	  nc_type attrtype;
	  size_t attlen;
	  char attname[256];
	  const int attstringlen = 8192; char attstring[8192];

	  nvatts = ncvars[ncvarid].natts;
	  for ( iatt = 0; iatt < nvatts; iatt++ )
	    {
	      attnum = ncvars[ncvarid].atts[iatt];
	      cdf_inq_attname(fileID, ncvarid, attnum, attname);
	      cdf_inq_attlen(fileID, ncvarid, attname, &attlen);
	      cdf_inq_atttype(fileID, ncvarid, attname, &attrtype);
	      if ( attrtype == NC_SHORT || attrtype == NC_INT )
		{
		  int *attint;
		  attint = (int *) malloc(attlen*sizeof(int));
		  cdfGetAttInt(fileID, ncvarid, attname, attlen, attint);
		  if ( attrtype == NC_SHORT )
		    vlistDefAttInt(vlistID, varID, attname, DATATYPE_INT16, (int)attlen, attint);
		  else
		    vlistDefAttInt(vlistID, varID, attname, DATATYPE_INT32, (int)attlen, attint);
		  if ( CDI_Debug )
		    printf("int: %s.%s = %d\n", ncvars[ncvarid].name, attname, attint[0]);
		  free(attint);
		}
	      else if ( attrtype == NC_FLOAT || attrtype == NC_DOUBLE )
		{
		  double *attflt;
		  attflt = (double *) malloc(attlen*sizeof(double));
		  cdfGetAttDouble(fileID, ncvarid, attname, attlen, attflt);
		  if ( attrtype == NC_FLOAT )
		    vlistDefAttFlt(vlistID, varID, attname, DATATYPE_FLT32, (int)attlen, attflt);
		  else
		    vlistDefAttFlt(vlistID, varID, attname, DATATYPE_FLT64, (int)attlen, attflt);
		  if ( CDI_Debug )
		    printf("flt: %s.%s = %g\n", ncvars[ncvarid].name, attname, attflt[0]);
		  free(attflt);
		}
	      else if ( attrtype == NC_CHAR )
		{
		  cdfGetAttText(fileID, ncvarid, attname, attstringlen-1, attstring);
		  attlen = 1 + strlen(attstring);
		  vlistDefAttTxt(vlistID, varID, attname, (int)attlen, attstring);
		  if ( CDI_Debug )
		    printf("txt: %s.%s = %s\n", ncvars[ncvarid].name, attname, attstring);
		}
	      else
		{
		  if ( CDI_Debug )
		    printf("att: %s.%s = unknown\n", ncvars[ncvarid].name, attname);
		}
	    }

	  free(ncvars[ncvarid].atts);
	}
    }

  if ( varids ) free(varids);

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vlistInqVarCode(vlistID, varID) == -varID-1 )
	{
	  const char *pname = vlistInqVarNamePtr(vlistID, varID);
	  size_t len = strlen(pname);
	  if ( len > 3 && isdigit((int) pname[3]) )
	    {
	      if ( memcmp("var", pname, 3) == 0 )
		{
		  vlistDefVarCode(vlistID, varID, atoi(pname+3));
		  vlistDestroyVarName(vlistID, varID);
		}
	    }
	  else if ( len > 4 && isdigit((int) pname[4]) )
	    {
	      if ( memcmp("code", pname, 4) == 0 )
		{
		  vlistDefVarCode(vlistID, varID, atoi(pname+4));
		  vlistDestroyVarName(vlistID, varID);
		}
	    }
	  else if ( len > 5 && isdigit((int) pname[5]) )
	    {
	      if ( memcmp("param", pname, 5) == 0 )
		{
		  int pnum = -1, pcat = 255, pdis = 255;
		  sscanf(pname+5, "%d.%d.%d", &pnum, &pcat, &pdis);
		  vlistDefVarParam(vlistID, varID, cdiEncodeParam(pnum, pcat, pdis));
		  vlistDestroyVarName(vlistID, varID);
		}
	    }
	}
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      instID  = vlistInqVarInstitut(vlistID, varID);
      modelID = vlistInqVarModel(vlistID, varID);
      tableID = vlistInqVarTable(vlistID, varID);
      code    = vlistInqVarCode(vlistID, varID);
      if ( cdiDefaultTableID != UNDEFID )
	{
	  if ( tableInqParNamePtr(cdiDefaultTableID, code) )
	    {
	      vlistDestroyVarName(vlistID, varID);
	      vlistDestroyVarLongname(vlistID, varID);
	      vlistDestroyVarUnits(vlistID, varID);

	      if ( tableID != UNDEFID )
		{
		  vlistDefVarName(vlistID, varID, tableInqParNamePtr(cdiDefaultTableID, code));
		  if ( tableInqParLongnamePtr(cdiDefaultTableID, code) )
		    vlistDefVarLongname(vlistID, varID, tableInqParLongnamePtr(cdiDefaultTableID, code));
		  if ( tableInqParUnitsPtr(cdiDefaultTableID, code) )
		    vlistDefVarUnits(vlistID, varID, tableInqParUnitsPtr(cdiDefaultTableID, code));
		}
	      else
		{
		  tableID = cdiDefaultTableID;
		}
	    }
	  if ( cdiDefaultModelID != UNDEFID )
	    modelID = cdiDefaultModelID;
	  if ( cdiDefaultInstID != UNDEFID )
	    instID = cdiDefaultInstID;
	}
      if ( instID  != UNDEFID ) vlistDefVarInstitut(vlistID, varID, instID);
      if ( modelID != UNDEFID ) vlistDefVarModel(vlistID, varID, modelID);
      if ( tableID != UNDEFID ) vlistDefVarTable(vlistID, varID, tableID);
    }
}
#endif

int cdfInqContents(int streamID)
{
#if  defined  (HAVE_LIBNETCDF)
  int ndims, nvars, ngatts, unlimdimid;
  int ncvarid;
  int ncdimid;
  int fileID;
  nc_type xtype;
  size_t ntsteps;
  int timedimid = -1;
  int *varids;
  int nvarids;
  size_t attlen;
  char attname[256];
  const int attstringlen = 8192; char attstring[8192];
  int iatt, timehasunits = FALSE;
  int time_has_bounds = FALSE;
  size_t len;
  int nc_nvars;
  int nvcth_id = UNDEFID, vcta_id = UNDEFID, vctb_id = UNDEFID;
  size_t vctsize = 0;
  int tableID;
  double *vct = NULL;
  int instID  = UNDEFID;
  int modelID = UNDEFID;
  int taxisID;
  int i;
  int nbdims;
  int calendar = UNDEFID;
  ncdim_t *ncdims;
  ncvar_t *ncvars;
  int vlistID;
  stream_t *streamptr;
  int format = 0;
  int ucla_les = FALSE;

  streamptr = stream_to_pointer(streamID);
  
  stream_check_ptr(__func__, streamptr);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  if ( CDI_Debug )
    Message("streamID = %d, fileID = %d", streamID, fileID);

#if  defined  (HAVE_NETCDF4)
  nc_inq_format(fileID, &format);
#endif

  cdf_inq(fileID, &ndims , &nvars, &ngatts, &unlimdimid);

  /* alloc ncdims */
  if ( ndims > 0 )
    ncdims = (ncdim_t *) malloc(ndims*sizeof(ncdim_t));
  else
    {
      Warning("ndims = %d", ndims);
      return (CDI_EUFSTRUCT);
    }
  
  /* alloc ncvars */
  if ( nvars > 0 )
    ncvars = (ncvar_t *) malloc(nvars*sizeof(ncvar_t));
  else
    {
      Warning("nvars = %d", nvars);
      return (CDI_EUFSTRUCT);
    }

  init_ncdims(ndims, ncdims);
  init_ncvars(nvars, ncvars);

  /* read global attributtes */
  for ( iatt = 0; iatt < ngatts; iatt++ )
    {
      cdf_inq_attname(fileID, NC_GLOBAL, iatt, attname);
      cdf_inq_atttype(fileID, NC_GLOBAL, attname, &xtype);
      cdf_inq_attlen(fileID, NC_GLOBAL, attname, &attlen);

      if ( xtype == NC_CHAR )
	{
	  cdfGetAttText(fileID, NC_GLOBAL, attname, attstringlen-1, attstring);

	  if ( attlen > 0 && attstring[0] != 0 )
	    {
	      if ( strcmp(attname, "history") == 0 )
		{
		  streamptr->historyID = iatt;
		}
	      else if ( strcmp(attname, "institution") == 0 )
		{
		  instID = institutInq(0, 0, NULL, attstring);
		  if ( instID == UNDEFID )
		    instID = institutDef(0, 0, NULL, attstring);
		}
	      else if ( strcmp(attname, "source") == 0 )
		{
		  modelID = modelInq(-1, 0, attstring);
		}
	      else if ( strcmp(attname, "Source") == 0 )
		{
		  if ( strncmp(attstring, "UCLA-LES", 8) == 0 )
		    ucla_les = TRUE;
		}
	      /*
	      else if ( strcmp(attname, "Conventions") == 0 )
		{
		}
	      */
	      else if ( strcmp(attname, "CDI") == 0 )
		{
		}
	      else if ( strcmp(attname, "CDO") == 0 )
		{
		}
	      else
		{
		  vlistDefAttTxt(vlistID, CDI_GLOBAL, attname, (int)attlen, attstring);
		}
	    }
	}
      else if ( xtype == NC_SHORT || xtype == NC_INT )
	{
	  int *attint;
	  attint = (int *) malloc(attlen*sizeof(int));
	  cdfGetAttInt(fileID, NC_GLOBAL, attname, attlen, attint);
	  if ( xtype == NC_SHORT )
	    vlistDefAttInt(vlistID, CDI_GLOBAL, attname, DATATYPE_INT16, (int)attlen, attint);
	  else
	    vlistDefAttInt(vlistID, CDI_GLOBAL, attname, DATATYPE_INT32, (int)attlen, attint);
	  free(attint);
	}
      else if ( xtype == NC_FLOAT || xtype == NC_DOUBLE )
	{
	  double *attflt;
	  attflt = (double *) malloc(attlen*sizeof(double));
	  cdfGetAttDouble(fileID, NC_GLOBAL, attname, attlen, attflt);
	  if ( xtype == NC_FLOAT )
	    vlistDefAttFlt(vlistID, CDI_GLOBAL, attname, DATATYPE_FLT32, (int)attlen, attflt);
	  else
	    vlistDefAttFlt(vlistID, CDI_GLOBAL, attname, DATATYPE_FLT64, (int)attlen, attflt);
	  free(attflt);
	}
    }

  /* find time dim */
  if ( unlimdimid >= 0 )
    timedimid = unlimdimid;
  else
    timedimid = cdfTimeDimID(fileID, ndims, nvars);

  streamptr->basetime.ncdimid = timedimid;

  if ( timedimid != UNDEFID )
    cdf_inq_dimlen(fileID, timedimid, &ntsteps);
  else
    ntsteps = 0;

  streamptr->ntsteps = ntsteps;

  if ( CDI_Debug ) Message("Number of timesteps = %d", streamptr->ntsteps);
  if ( CDI_Debug ) Message("Time dimid = %d", streamptr->basetime.ncdimid);

  /* read ncdims */
  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      cdf_inq_dimlen(fileID, ncdimid, &ncdims[ncdimid].len);
      cdf_inq_dimname(fileID, ncdimid, ncdims[ncdimid].name);
      if ( timedimid == ncdimid )
	ncdims[ncdimid].dimtype = T_AXIS;
    }


  /* scan attributes of all variables */
  cdfScanVarAttributes(fileID, nvars, ncvars, ncdims, timedimid, modelID, format);


  if ( CDI_Debug ) printNCvars(ncvars, nvars);

  /* find coordinate vars */
  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
	{
	  if ( ncvars[ncvarid].ndims == 1 )
	    {
	      if ( timedimid != UNDEFID && timedimid == ncvars[ncvarid].dimids[0] )
		{
		  if ( ncvars[ncvarid].isvar != FALSE ) cdfSetVar(ncvars, ncvarid, TRUE);
		}
	      else
		{
		  if ( ncvars[ncvarid].isvar != TRUE ) cdfSetVar(ncvars, ncvarid, FALSE);
		}
	      // if ( ncvars[ncvarid].isvar != TRUE ) cdfSetVar(ncvars, ncvarid, FALSE);

	      if ( ncdimid == ncvars[ncvarid].dimids[0] && ncdims[ncdimid].ncvarid == UNDEFID )
		if ( strcmp(ncvars[ncvarid].name, ncdims[ncdimid].name) == 0 )
		  {
		    ncdims[ncdimid].ncvarid = ncvarid;
		    ncvars[ncvarid].isvar = FALSE;
		  }
	    }
	}
    }

  /* find time vars */
  if ( timedimid != UNDEFID )
    {
      int ltimevar = FALSE;

      if ( ncdims[timedimid].ncvarid != UNDEFID )
	{
	  streamptr->basetime.ncvarid = ncdims[timedimid].ncvarid;
	  ltimevar = TRUE;
	}

      for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
	if ( ncvarid != streamptr->basetime.ncvarid &&
	     ncvars[ncvarid].ndims == 1 &&
	     timedimid == ncvars[ncvarid].dimids[0] &&
	     ncvars[ncvarid].xtype != NC_CHAR &&
	     isTimeUnits(ncvars[ncvarid].units) )
	  {
	    ncvars[ncvarid].isvar = FALSE;

	    if ( !ltimevar )
	      {
		streamptr->basetime.ncvarid = ncvarid;
		ltimevar = TRUE;
		if ( CDI_Debug ) 
		  fprintf(stderr, "timevar %s\n", ncvars[ncvarid].name);
	      }
	    else
	      {
		if ( CDI_Debug )
		  fprintf(stderr, "skip timevar %s\n", ncvars[ncvarid].name);
	      }
	  }

      if ( ltimevar == FALSE ) /* search for WRF time description */
	{
	  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
	    if ( ncvarid != streamptr->basetime.ncvarid &&
		 ncvars[ncvarid].ndims == 2 &&
		 timedimid == ncvars[ncvarid].dimids[0] &&
		 ncvars[ncvarid].xtype == NC_CHAR &&
		 ncdims[ncvars[ncvarid].dimids[1]].len == 19 )
	      {
		streamptr->basetime.ncvarid = ncvarid;
		streamptr->basetime.lwrf    = TRUE;
		break;
	      }
	}

      /* time varID */
      ncvarid = streamptr->basetime.ncvarid;

      if ( ncvarid == UNDEFID )
	Warning("Variable >time< not found!");
      else if ( streamptr->basetime.lwrf == FALSE )
	{
	  if ( ncvars[ncvarid].units[0] != 0 )
	    timehasunits = TRUE;

	  if ( ncvars[ncvarid].bounds != UNDEFID )
	    {
	      nbdims = ncvars[ncvars[ncvarid].bounds].ndims;
	      if ( nbdims == 2 )
		{
		  len = ncdims[ncvars[ncvars[ncvarid].bounds].dimids[nbdims-1]].len;
		  if ( (int)len == 2 && timedimid == ncvars[ncvars[ncvarid].bounds].dimids[0] )
		    {
		      time_has_bounds = TRUE;
		      streamptr->basetime.ncvarboundsid = ncvars[ncvarid].bounds;
		    }
		}
	    }
	}
    }

  /* check ncvars */
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( timedimid != UNDEFID )
	if ( ncvars[ncvarid].isvar == -1 &&
	     ncvars[ncvarid].ndims > 1   &&
	     timedimid == ncvars[ncvarid].dimids[0] )
	  cdfSetVar(ncvars, ncvarid, TRUE);

      if ( ncvars[ncvarid].isvar == -1 && ncvars[ncvarid].ndims == 0 )
	cdfSetVar(ncvars, ncvarid, FALSE);

      if ( ncvars[ncvarid].isvar == -1 && ncvars[ncvarid].ndims > 1 )
	cdfSetVar(ncvars, ncvarid, TRUE);

      if ( ncvars[ncvarid].isvar == -1 )
	{
	  ncvars[ncvarid].isvar = 0;
	  Warning("Variable %s has unknown type, skipped!", ncvars[ncvarid].name);
	  continue;
	}

      if ( ncvars[ncvarid].ndims > 4 )
	{
	  ncvars[ncvarid].isvar = 0;
	  Warning("%d dimensional variables unsupported. Skip variable %s",
		ncvars[ncvarid].ndims, ncvars[ncvarid].name);
	  continue;
	}

      if ( ncvars[ncvarid].xtype == NC_CHAR )
	{
	  ncvars[ncvarid].isvar = 0;
	  continue;
	}

      if ( cdfInqDatatype(ncvars[ncvarid].xtype, ncvars[ncvarid].lunsigned) == -1 )
	{
	  ncvars[ncvarid].isvar = 0;
	  Warning("Variable %s has an unsupported data type, skipped!", ncvars[ncvarid].name);
	  continue;
	}

      if ( timedimid != UNDEFID && ntsteps == 0 && ncvars[ncvarid].ndims > 0 )
	{
	  if ( timedimid == ncvars[ncvarid].dimids[0] )
	    {
	      ncvars[ncvarid].isvar = 0;
	      Warning("Number of time steps undefined, skipped variable %s!", ncvars[ncvarid].name);
	      continue;
	    }
	}
    }


  /* verify coordinate vars - first scan (dimname == varname) */
  verify_coordinate_vars_1(ndims, ncdims, ncvars, timedimid);

  /* verify coordinate vars - second scan (all other variables) */
  verify_coordinate_vars_2(nvars, ncvars);

  if ( CDI_Debug ) printNCvars(ncvars, nvars);

  if ( ucla_les == TRUE )
    {
      for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
	{
	  ncvarid = ncdims[ncdimid].ncvarid;
	  if ( ncvarid != -1 )
	    {
	      if ( ncdims[ncdimid].dimtype == UNDEFID && ncvars[ncvarid].units[0] == 'm' )
		{
		  if      ( ncvars[ncvarid].name[0] == 'x' ) ncdims[ncdimid].dimtype = X_AXIS;
		  else if ( ncvars[ncvarid].name[0] == 'y' ) ncdims[ncdimid].dimtype = Y_AXIS;
		  else if ( ncvars[ncvarid].name[0] == 'z' ) ncdims[ncdimid].dimtype = Z_AXIS;
		}
	    }
	}
    }
  /*
  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      ncvarid = ncdims[ncdimid].ncvarid;
      if ( ncvarid != -1 )
	{
	  printf("coord var %d %s %s\n", ncvarid, ncvars[ncvarid].name, ncvars[ncvarid].units);
	  if ( ncdims[ncdimid].dimtype == X_AXIS )
	    printf("coord var %d %s is x dim\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncdims[ncdimid].dimtype == Y_AXIS )
	    printf("coord var %d %s is y dim\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncdims[ncdimid].dimtype == Z_AXIS )
	    printf("coord var %d %s is z dim\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncdims[ncdimid].dimtype == T_AXIS )
	    printf("coord var %d %s is t dim\n", ncvarid, ncvars[ncvarid].name);

	  if ( ncvars[ncvarid].islon )
	    printf("coord var %d %s is lon\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncvars[ncvarid].islat )
	    printf("coord var %d %s is lat\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncvars[ncvarid].islev )
	    printf("coord var %d %s is lev\n", ncvarid, ncvars[ncvarid].name);
	}
    }
  */
  /* set dim type */
  setDimType(nvars, ncvars, ncdims);

  /* Set coordinate varids (att: associate)  */
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar == TRUE && ncvars[ncvarid].ncoordvars )
	{
	  /* ndims = ncvars[ncvarid].ndims; */
	  ndims = ncvars[ncvarid].ncoordvars;
	  for ( i = 0; i < ndims; i++ )
	    {
	      if ( ncvars[ncvars[ncvarid].coordvarids[i]].islon )
		ncvars[ncvarid].xvarid = ncvars[ncvarid].coordvarids[i];
	      else if ( ncvars[ncvars[ncvarid].coordvarids[i]].islat )
		ncvars[ncvarid].yvarid = ncvars[ncvarid].coordvarids[i];
	      else if ( ncvars[ncvars[ncvarid].coordvarids[i]].islev )
		ncvars[ncvarid].zvarid = ncvars[ncvarid].coordvarids[i];
	    }
	}
    }
  
  if ( CDI_Debug ) printNCvars(ncvars, nvars);


  /* define all grids */
  define_all_grids(streamptr, fileID, vlistID, ncdims, nvars, ncvars, timedimid);


  /* find VCT */
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].ndims == 1 )
	{
	  if ( memcmp(ncvars[ncvarid].name, "hyai", 4) == 0 )
	    {
	      vcta_id = ncvarid;
	      nvcth_id = ncvars[ncvarid].dimids[0];
	      continue;
	    }
	  if ( memcmp(ncvars[ncvarid].name, "hybi", 4) == 0 )
	    {
	      vctb_id = ncvarid;
	      nvcth_id = ncvars[ncvarid].dimids[0];
	      continue;
	    }
	}
    }

  /* read VCT */
  if ( nvcth_id != UNDEFID && vcta_id != UNDEFID && vctb_id != UNDEFID )
    {
      vctsize = ncdims[nvcth_id].len;
      vctsize *= 2;
      vct = (double *) malloc(vctsize*sizeof(double));
      cdf_get_var_double(fileID, vcta_id, vct);
      cdf_get_var_double(fileID, vctb_id, vct+vctsize/2);
    }


  /* define all zaxes */
  define_all_zaxes(streamptr, fileID, vlistID, ncdims, nvars, ncvars, vctsize, vct);


  if ( vct ) free(vct);

  /* select vars */
  varids = (int *) malloc(nvars*sizeof(int));
  nvarids = 0;
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    if ( ncvars[ncvarid].isvar == TRUE ) varids[nvarids++] = ncvarid;

  nc_nvars = nvars;
  nvars = nvarids;

  if ( CDI_Debug ) Message("time varid = %d", streamptr->basetime.ncvarid);
  if ( CDI_Debug ) Message("ntsteps = %d", streamptr->ntsteps);
  if ( CDI_Debug ) Message("nvars = %d", nvars);

  if ( nvars == 0 ) return (CDI_EUFSTRUCT);


  /* define all data variables */
  define_all_vars(fileID, streamID, vlistID, instID, modelID, tableID, varids, ncdims, nvars, ncvars);


  cdiCreateTimesteps(streamID);

  /* time varID */
  ncvarid = streamptr->basetime.ncvarid;

  if ( timehasunits )
    {
      taxis_t *taxis;
      taxis = &streamptr->tsteps[0].taxis;

      cdfGetAttText(fileID, ncvarid, "units", attstringlen-1, attstring);
      if ( splitBasetime(attstring, taxis) == 1 )
	streamptr->basetime.ncvarid = UNDEFID;
    }

  if ( time_has_bounds )
    streamptr->tsteps[0].taxis.has_bounds = TRUE;

  if ( ncvarid != -1 )
    if ( ncvars[ncvarid].calendar == TRUE )
      {
	cdfGetAttText(fileID, ncvarid, "calendar", attstringlen-1, attstring);
	strtolower(attstring);

	if ( memcmp(attstring, "standard", 8)  == 0 ||
	     memcmp(attstring, "gregorian", 9) == 0 )
	  calendar = CALENDAR_STANDARD;
	else if ( memcmp(attstring, "none", 4) == 0 )
	  calendar = CALENDAR_NONE;
	else if ( memcmp(attstring, "proleptic", 9) == 0 )
	  calendar = CALENDAR_PROLEPTIC;
	else if ( memcmp(attstring, "360", 3) == 0 )
	  calendar = CALENDAR_360DAYS;
	else if ( memcmp(attstring, "365", 3) == 0 ||
		  memcmp(attstring, "noleap", 6)  == 0 )
	  calendar = CALENDAR_365DAYS;
	else if ( memcmp(attstring, "366", 3)  == 0 ||
		  memcmp(attstring, "all_leap", 8) == 0 )
	  calendar = CALENDAR_366DAYS;
	else
	  Warning("calendar >%s< unsupported!", attstring);
      }

  if ( streamptr->tsteps[0].taxis.type == TAXIS_RELATIVE )
    taxisID = taxisCreate(TAXIS_RELATIVE);
  else
    {
      taxisID = taxisCreate(TAXIS_ABSOLUTE);
      if ( !timehasunits )
	{
	  taxisDefTunit(taxisID, TUNIT_DAY);
	  streamptr->tsteps[0].taxis.unit = TUNIT_DAY;
	}
    }

  if ( calendar != UNDEFID )
    {
      taxis_t *taxis;
      taxis = &streamptr->tsteps[0].taxis;

      taxis->calendar = calendar;
      taxisDefCalendar(taxisID, calendar);
    }
  else if ( streamptr->tsteps[0].taxis.type == TAXIS_RELATIVE )
    {
      taxis_t *taxis;

      calendar = CALENDAR_STANDARD;

      taxis = &streamptr->tsteps[0].taxis;

      taxis->calendar = calendar;
      taxisDefCalendar(taxisID, calendar);
    }

  vlistDefTaxis(vlistID, taxisID);

  streamptr->curTsID = 0;
  streamptr->rtsteps = 1;

  (void) cdfInqTimestep(streamID, 0);

  cdfCreateRecords(streamID, 0);

  /* free ncdims */
  free (ncdims);

  /* free ncvars */
  free (ncvars);

#endif

  return (0);
}


int cdfInqTimestep(int streamID, int tsID)
{
  long nrecs = 0;
#if  defined  (HAVE_LIBNETCDF)
  double timevalue;
  int nctimevarid;
  int nctimeboundsid;
  int fileID;
  size_t index;
  taxis_t *taxis;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( CDI_Debug )
    Message("streamID = %d  tsID = %d", streamID, tsID);

  stream_check_ptr(__func__, streamptr);

  if ( tsID < 0 ) Error("unexpected tsID = %d", tsID);

  if ( tsID < streamptr->ntsteps && streamptr->ntsteps > 0 )
    {
      cdfCreateRecords(streamID, tsID);

      taxis = &streamptr->tsteps[tsID].taxis;
      if ( tsID > 0 )
	ptaxisCopy(taxis, &streamptr->tsteps[0].taxis);

      timevalue = tsID;

      nctimevarid = streamptr->basetime.ncvarid;
      if ( nctimevarid != UNDEFID )
	{
	  fileID = streamInqFileID(streamID);
	  index  = tsID;

	  if ( streamptr->basetime.lwrf )
	    {
	      size_t start[2], count[2];
	      char stvalue[32];
	      start[0] = index; start[1] = 0;
	      count[0] = 1; count[1] = 19;
	      stvalue[0] = 0;
	      cdf_get_vara_text(fileID, nctimevarid, start, count, stvalue);
	      stvalue[19] = 0;
	      {
		int year = 1, month = 1, day = 1 , hour = 0, minute = 0, second = 0;
		if ( strlen(stvalue) == 19 )
		  sscanf(stvalue, "%d-%d-%d_%d:%d:%d", &year, &month, &day, &hour, &minute, &second);
		taxis->vdate = cdiEncodeDate(year, month, day);
		taxis->vtime = cdiEncodeTime(hour, minute, second);
		taxis->type = TAXIS_ABSOLUTE;
	      }
	    }
	  else
	    {
	      cdf_get_var1_double(fileID, nctimevarid, &index, &timevalue);
	      cdiDecodeTimeval(timevalue, taxis, &taxis->vdate, &taxis->vtime);
	    }

	  nctimeboundsid = streamptr->basetime.ncvarboundsid;
	  if ( nctimeboundsid != UNDEFID )
	    {
	      size_t start[2], count[2];
	      start[0] = tsID; count[0] = 1; start[1] = 0; count[1] = 1;
	      cdf_get_vara_double(fileID, nctimeboundsid, start, count, &timevalue);

	      cdiDecodeTimeval(timevalue, taxis, &taxis->vdate_lb, &taxis->vtime_lb);

	      start[0] = tsID; count[0] = 1; start[1] = 1; count[1] = 1;
	      cdf_get_vara_double(fileID, nctimeboundsid, start, count, &timevalue);

	      cdiDecodeTimeval(timevalue, taxis, &taxis->vdate_ub, &taxis->vtime_ub);
	    }
	}
    }

  streamptr->curTsID = tsID;
  nrecs = streamptr->tsteps[tsID].nrecs;

#endif
  return ((int) nrecs);
}


void cdfEndDef(int streamID)
{
#if  defined  (HAVE_LIBNETCDF)
  int varID, ncvarid;
  int nvars;
  int fileID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  fileID  = streamInqFileID(streamID);

  cdfDefGlobalAtts(streamID);
  cdfDefLocalAtts(streamID);
  if ( streamptr->accessmode == 0 )
    {
      nvars =  streamptr->nvars;

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      for ( varID = 0; varID < nvars; varID++ )
	ncvarid = cdfDefVar(streamID, varID);

      if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);

      streamptr->accessmode = 1;
    }
#endif
}


void cdfDefInstitut(int streamID)
{
#if  defined  (HAVE_LIBNETCDF)
  int fileID, instID;
  char *longname;
  size_t len;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);
  instID  = vlistInqInstitut(vlistID);

  if ( instID != UNDEFID )
    {
      longname = institutInqLongnamePtr(instID);
      if ( longname )
	{
	  len = strlen(longname);
	  if ( len > 0 )
	    {
	      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);
	      cdf_put_att_text(fileID, NC_GLOBAL, "institution", len, longname);
	      if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
	    }
	}
    }
#endif
}

void cdfDefSource(int streamID)
{
#if  defined  (HAVE_LIBNETCDF)
  int fileID, modelID;
  char *longname;
  size_t len;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);
  modelID = vlistInqModel(vlistID);

  if ( modelID != UNDEFID )
    {
      longname = modelInqNamePtr(modelID);
      if ( longname )
	{
	  len = strlen(longname);
	  if ( len > 0 )
	    {
	      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);
	      cdf_put_att_text(fileID, NC_GLOBAL, "source", len, longname);
	      if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
	    }
	}
    }
#endif
}


void cdfDefGlobalAtts(int streamID)
{
#if  defined  (HAVE_LIBNETCDF)
  int fileID, vlistID;
  vlist_t *vlistptr;
  stream_t *streamptr;
  int natts;

  streamptr = stream_to_pointer(streamID);

  if ( streamptr->globalatts ) return;

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  vlistptr = vlist_to_pointer(vlistID);

  cdfDefSource(streamID);
  cdfDefInstitut(streamID);

  vlistInqNatts(vlistID, CDI_GLOBAL, &natts);

  if ( natts > 0 && streamptr->ncmode == 2 ) cdf_redef(fileID);

  defineAttributes(vlistID, CDI_GLOBAL, fileID, NC_GLOBAL);

  if ( natts > 0 && streamptr->ncmode == 2 ) cdf_enddef(fileID);

  streamptr->globalatts = 1;
#endif
}


void cdfDefLocalAtts(int streamID)
{
#if  defined  (HAVE_LIBNETCDF)
  int varID, instID, fileID;
  char *name;
  size_t len;
  int ncvarid;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  if ( streamptr->localatts ) return;
  if ( vlistInqInstitut(vlistID) != UNDEFID ) return;

  streamptr->localatts = 1;

  if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

  for ( varID = 0; varID < streamptr->nvars; varID++ )
    {
      instID = vlistInqVarInstitut(vlistID, varID);
      if ( instID != UNDEFID )
	{
          ncvarid = streamptr->vars[varID].ncvarid;
  	  name = institutInqNamePtr(instID);
	  if ( name )
	    {
	      len = strlen(name);
	      cdf_put_att_text(fileID, ncvarid, "institution", len, name);
	    }
	}
      }

  if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
#endif
}

void cdfDefHistory(int streamID, int size, char *history)
{
#if  defined  (HAVE_LIBNETCDF)
  int ncid;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  ncid = streamptr->fileID;
  cdf_put_att_text(ncid, NC_GLOBAL, "history", (size_t) size, history);
#endif
}

int cdfInqHistorySize(int streamID)
{
  size_t size = 0;
#if  defined  (HAVE_LIBNETCDF)
  int ncid;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  ncid = streamptr->fileID;
  if ( streamptr->historyID != UNDEFID )
    cdf_inq_attlen(ncid, NC_GLOBAL, "history", &size);

#endif
  return ((int) size);
}


void cdfInqHistoryString(int streamID, char *history)
{
#if  defined  (HAVE_LIBNETCDF)
  int ncid;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  ncid = streamptr->fileID;
  if ( streamptr->historyID != UNDEFID )
    cdf_get_att_text(ncid, NC_GLOBAL, "history", history);

#endif
}


void cdfDefVars(int streamID)
{
#if  defined  (HAVE_LIBNETCDF)
  int index, gridID, zaxisID, vlistID;
  int nvars, ngrids, nzaxis;
  /*
  int  ncvarid;
  */
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  if ( vlistID == UNDEFID )
    Error("Internal problem! vlist undefined for streamID %d", streamID);

  nvars  = vlistNvars(vlistID);
  ngrids = vlistNgrids(vlistID);
  nzaxis = vlistNzaxis(vlistID);
  /*
  if ( vlistHasTime(vlistID) ) cdfDefTime(streamID);
  */
  for ( index = 0; index < ngrids; index++ )
    {
      gridID = vlistGrid(vlistID, index);
      cdfDefGrid(streamID, gridID);
    }

  for ( index = 0; index < nzaxis; index++ )
    {
      zaxisID = vlistZaxis(vlistID, index);
      if ( streamptr->zaxisID[index] == UNDEFID ) cdfDefZaxis(streamID, zaxisID);
    }
  /*
    define time first!!!
  for (varID = 0; varID < nvars; varID++ )
    {
      ncvarid = cdfDefVar(streamID, varID);
    }
  */
#endif
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
