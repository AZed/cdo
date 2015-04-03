#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdarg.h>
#include <ctype.h>

#include "dmemory.h"

#include "cdi.h"
#include "cdi_cksum.h"
#include "cdi_int.h"
#include "gribapi.h"
#ifdef HAVE_LIBNETCDF
#include "stream_cdf.h"
#endif
#include "namespace.h"
#include "serialize.h"
#include "resource_handle.h"
#include "resource_unpack.h"

#if  defined  (HAVE_LIBCGRIBEX)
#include "cgribex.h"
#endif

extern int cdiPioSerialOpenFileMap(int streamID);

int cdiDefaultCalendar = CALENDAR_PROLEPTIC;

int cdiDefaultInstID   = CDI_UNDEFID;
int cdiDefaultModelID  = CDI_UNDEFID;
int cdiDefaultTableID  = CDI_UNDEFID;
//int cdiNcMissingValue  = CDI_UNDEFID;
int cdiNcChunksizehint = CDI_UNDEFID;
int cdiChunkType       = CHUNK_GRID;
int cdiSplitLtype105   = CDI_UNDEFID;

int cdiIgnoreAttCoordinates = FALSE;
int cdiIgnoreValidRange     = FALSE;
int cdiSkipRecords          = 0;
int cdiInventoryMode        = 1;
size_t CDI_netcdf_hdr_pad   = 0UL;

char *cdiPartabPath   = NULL;
int   cdiPartabIntern = 1;

double cdiDefaultMissval = -9.E33;

const char Filetypes[][9] = {
  "UNKNOWN",
  "GRIB",
  "GRIB2",
  "netCDF",
  "netCDF2",
  "netCDF4",
  "netCDF4c",
  "SERVICE",
  "EXTRA",
  "IEG",
  "HDF5",
};

#undef  UNDEFID
#define UNDEFID  CDI_UNDEFID


int CDI_Debug   = 0;    /* If set to 1, debugging           */

static int  STREAM_Debug = 0;   /* If set to 1, debugging */

int cdiGribApiDebug     = 0;
int cdiDefaultLeveltype = -1;
static int cdiDataUnreduced = 0;
static int cdiSortName = 0;
static int cdiHaveMissval = 0;


static int    streamCompareP ( void * streamptr1, void * streamptr2 );
static void   streamDestroyP ( void * streamptr );
static void   streamPrintP   ( void * streamptr, FILE * fp );
static int    streamGetPackSize ( void * streamptr, void *context);
static void   streamPack        ( void * streamptr, void * buff, int size, int * position, void *context );
static int    streamTxCode      ( void );

const resOps streamOps = {
  streamCompareP,
  streamDestroyP,
  streamPrintP,
  streamGetPackSize,
  streamPack,
  streamTxCode
};

long cdiGetenvInt(char *envName)
{
  char *envString;
  long envValue = -1;
  long fact = 1;

  envString = getenv(envName);

  if ( envString )
    {
      int loop, len;

      len = (int) strlen(envString);
      for ( loop = 0; loop < len; loop++ )
	{
	  if ( ! isdigit((int) envString[loop]) )
	    {
	      switch ( tolower((int) envString[loop]) )
		{
		case 'k':  fact = 1024;        break;
		case 'm':  fact = 1048576;     break;
		case 'g':  fact = 1073741824;  break;
		default:
		  fact = 0;
		  Message("Invalid number string in %s: %s", envName, envString);
		  Warning("%s must comprise only digits [0-9].",envName);
		  break;
		}
	      break;
	    }
	}

      if ( fact ) envValue = fact*atol(envString);

      if ( CDI_Debug ) Message("set %s to %ld", envName, envValue);
    }

  return (envValue);
}

static
void cdiSetChunk(const char *chunkAlgo)
{
  //char *pch;
  //size_t len = strlen(chunkAlgo);
  int algo = -1;

  if      ( strcmp("auto",  chunkAlgo)   == 0 ) algo = CHUNK_AUTO;
  else if ( strcmp("grid",  chunkAlgo)   == 0 ) algo = CHUNK_GRID;
  else if ( strcmp("lines", chunkAlgo)   == 0 ) algo = CHUNK_LINES;
  /*
  else if ( (pch = strstr(chunkAlgo,"x")) != 0 )
    {
      int ix, iy;
      ix = atoi(chunkAlgo);
      iy = atoi(pch+1);
      if ( ix > 0 && iy > 0 )
        {
          cdiChunkX = ix;
          cdiChunkY = iy;
          algo = CHUNK_USER;
        }
      else
        Warning("Invalid environment variable CDI_CHUNK_ALGO: %s", chunkAlgo);
    }
  */
  else
    Warning("Invalid environment variable CDI_CHUNK_ALGO: %s", chunkAlgo);

  if ( algo != -1 )
    {
      cdiChunkType = algo;
      if ( CDI_Debug ) Message("set ChunkAlgo to %s", chunkAlgo);
    }
}


void cdiInitialize(void)
{
  static int Init_CDI = FALSE;
  char *envString;
  long value;

  if ( ! Init_CDI )
    {
      Init_CDI = TRUE;

#if  defined  (HAVE_LIBCGRIBEX)
      gribFixZSE(1);   // 1: Fix ZeroShiftError of simple packed spherical harmonics
      gribSetConst(1); // 1: Don't pack constant fields on regular grids
#endif

      value = cdiGetenvInt("CDI_DEBUG");
      if ( value >= 0 ) CDI_Debug = (int) value;

      value = cdiGetenvInt("CDI_GRIBAPI_DEBUG");
      if ( value >= 0 ) cdiGribApiDebug = (int) value;

      value = cdiGetenvInt("CDI_REGULARGRID");
      if ( value >= 0 ) cdiDataUnreduced = (int) value;

      value = cdiGetenvInt("CDI_SORTNAME");
      if ( value >= 0 ) cdiSortName = (int) value;

      value = cdiGetenvInt("CDI_HAVE_MISSVAL");
      if ( value >= 0 ) cdiHaveMissval = (int) value;

      value = cdiGetenvInt("CDI_LEVELTYPE");
      if ( value >= 0 ) cdiDefaultLeveltype = (int) value;

      value = cdiGetenvInt("CDI_NETCDF_HDR_PAD");
      if ( value >= 0 ) CDI_netcdf_hdr_pad = (size_t) value;

      envString = getenv("CDI_MISSVAL");
      if ( envString ) cdiDefaultMissval = atof(envString);
      /*
      envString = getenv("NC_MISSING_VALUE");
      if ( envString ) cdiNcMissingValue = atoi(envString);
      */
      envString = getenv("NC_CHUNKSIZEHINT");
      if ( envString ) cdiNcChunksizehint = atoi(envString);

      envString = getenv("CDI_CHUNK_ALGO");
      if ( envString ) cdiSetChunk(envString);

      envString = getenv("SPLIT_LTYPE_105");
      if ( envString ) cdiSplitLtype105 = atoi(envString);

      envString = getenv("IGNORE_ATT_COORDINATES");
      if ( envString ) cdiIgnoreAttCoordinates = atoi(envString);

      envString = getenv("IGNORE_VALID_RANGE");
      if ( envString ) cdiIgnoreValidRange = atoi(envString);

      envString = getenv("CDI_SKIP_RECORDS");
      if ( envString )
	{
	  cdiSkipRecords = atoi(envString);
	  cdiSkipRecords = cdiSkipRecords > 0 ? cdiSkipRecords : 0;
	}

      envString = getenv("CDI_INVENTORY_MODE");
      if ( envString )
	{
	  if ( strncmp(envString, "time", 4) == 0 )
	    {
	      cdiInventoryMode = 2;
	      if ( CDI_Debug )
		Message("Inventory mode was set to timestep!");
	    }
	}

      envString = getenv("CDI_CALENDAR");
      if ( envString )
	{
	  if      ( strncmp(envString, "standard", 8) == 0 )
	    cdiDefaultCalendar = CALENDAR_STANDARD;
	  else if ( strncmp(envString, "proleptic", 9) == 0 )
	    cdiDefaultCalendar = CALENDAR_PROLEPTIC;
	  else if ( strncmp(envString, "360days", 7) == 0 )
	    cdiDefaultCalendar = CALENDAR_360DAYS;
	  else if ( strncmp(envString, "365days", 7) == 0 )
	    cdiDefaultCalendar = CALENDAR_365DAYS;
	  else if ( strncmp(envString, "366days", 7) == 0 )
	    cdiDefaultCalendar = CALENDAR_366DAYS;
	  else if ( strncmp(envString, "none", 4) == 0 )
	    cdiDefaultCalendar = CALENDAR_NONE;

	  if ( CDI_Debug )
	    Message("Default calendar set to %s!", envString);
	}
#if  defined  (HAVE_LIBCGRIBEX)
      gribSetCalendar(cdiDefaultCalendar);
#endif

      envString = getenv("PARTAB_INTERN");
      if ( envString ) cdiPartabIntern = atoi(envString);

      envString = getenv("PARTAB_PATH");
      if ( envString ) cdiPartabPath = strdup(envString);

      envString = getenv("STREAM_DEBUG");
      if ( envString ) STREAM_Debug = atoi(envString);
    }
}


const char *strfiletype(int filetype)
{
  const char *name;
  int size = (int) (sizeof(Filetypes)/sizeof(char *));

  if ( filetype > 0 && filetype < size )
    name = Filetypes[filetype];
  else
    name = Filetypes[0];

  return (name);
}


stream_t *stream_to_pointer(int idx)
{
  return ( stream_t *) reshGetVal ( idx, &streamOps );
}

static
void streamDefaultValue ( stream_t * streamptr )
{
  int i;

  streamptr->self              = UNDEFID;
  streamptr->accesstype        = UNDEFID;
  streamptr->accessmode        = 0;
  streamptr->filetype          = UNDEFID;
  streamptr->byteorder         = UNDEFID;
  streamptr->fileID            = 0;
  streamptr->dimgroupID        = UNDEFID;
  streamptr->filemode          = 0;
  streamptr->numvals           = 0;
  streamptr->filename          = NULL;
  streamptr->record            = NULL;
  streamptr->varsAllocated     = 0;
  streamptr->nrecs             = 0;
  streamptr->nvars             = 0;
  streamptr->vars              = NULL;
  streamptr->varinit           = 0;
  streamptr->ncmode            = 0;
  streamptr->curTsID           = UNDEFID;
  streamptr->rtsteps           = 0;
  streamptr->ntsteps           = UNDEFID;
  streamptr->numTimestep       = 0;
  streamptr->tsteps            = NULL;
  streamptr->tstepsTableSize   = 0;
  streamptr->tstepsNextID      = 0;
  streamptr->historyID         = UNDEFID;
  streamptr->vlistID           = UNDEFID;
  streamptr->globalatts        = 0;
  streamptr->localatts         = 0;
  streamptr->vct.ilev          = 0;
  streamptr->vct.mlev          = 0;
  streamptr->vct.ilevID        = UNDEFID;
  streamptr->vct.mlevID        = UNDEFID;
  streamptr->unreduced         = cdiDataUnreduced;
  streamptr->sortname          = cdiSortName;
  streamptr->have_missval      = cdiHaveMissval;
  streamptr->comptype          = COMPRESS_NONE;
  streamptr->complevel         = 0;

  basetimeInit(&streamptr->basetime);

  for ( i = 0; i < MAX_GRIDS_PS; i++ ) streamptr->xdimID[i]   = UNDEFID;
  for ( i = 0; i < MAX_GRIDS_PS; i++ ) streamptr->ydimID[i]   = UNDEFID;
  for ( i = 0; i < MAX_ZAXES_PS; i++ ) streamptr->zaxisID[i]  = UNDEFID;
  for ( i = 0; i < MAX_GRIDS_PS; i++ ) streamptr->ncxvarID[i] = UNDEFID;
  for ( i = 0; i < MAX_GRIDS_PS; i++ ) streamptr->ncyvarID[i] = UNDEFID;
  for ( i = 0; i < MAX_GRIDS_PS; i++ ) streamptr->ncavarID[i] = UNDEFID;

  streamptr->curfile           = 0;
  streamptr->nfiles            = 0;
  streamptr->fnames            = NULL;

  streamptr->gribContainers    = NULL;
  streamptr->vlistIDorig       = UNDEFID;
}


stream_t *stream_new_entry(void)
{
  stream_t *streamptr;

  cdiInitialize(); /* ***************** make MT version !!! */

  streamptr = (stream_t *) xmalloc(sizeof(stream_t));
  streamDefaultValue ( streamptr );
  streamptr->self = reshPut (( void * ) streamptr, &streamOps );

  return streamptr;
}


void stream_delete_entry(stream_t *streamptr)
{
  int idx;

  xassert ( streamptr );

  idx = streamptr->self;
  free ( streamptr );
  reshRemove ( idx, &streamOps );

  if ( STREAM_Debug )
    Message("Removed idx %d from stream list", idx);
}


void stream_check_ptr(const char *caller, stream_t *streamptr)
{
  if ( streamptr == NULL )
    Errorc("stream undefined!");
}


int streamSize(void)
{
  return reshCountType ( &streamOps );
}


void cdiDefGlobal(const char *string, int val)
{
  if      ( strcmp(string, "REGULARGRID")      == 0 ) cdiDataUnreduced = val;
  else if ( strcmp(string, "GRIBAPI_DEBUG")    == 0 ) cdiGribApiDebug = val;
  else if ( strcmp(string, "SORTNAME")         == 0 ) cdiSortName = val;
  else if ( strcmp(string, "HAVE_MISSVAL")     == 0 ) cdiHaveMissval = val;
  else if ( strcmp(string, "NC_CHUNKSIZEHINT") == 0 ) cdiNcChunksizehint = val;
  else if ( strcmp(string, "NETCDF_HDR_PAD")   == 0 ) CDI_netcdf_hdr_pad = (size_t) val;
  else Warning("Unsupported global key: %s", string);
}


void cdiDefMissval(double missval)
{
  cdiInitialize();

  cdiDefaultMissval = missval;
}


double cdiInqMissval(void)
{
  cdiInitialize();

  return (cdiDefaultMissval);
}


void vlist_check_contents(int vlistID)
{
  int index, nzaxis, zaxisID;

  nzaxis = vlistNzaxis(vlistID);

  for ( index = 0; index < nzaxis; index++ )
    {
      zaxisID = vlistZaxis(vlistID, index);
      if ( zaxisInqType(zaxisID) == ZAXIS_GENERIC )
	cdiCheckZaxis(zaxisID);
    }
}


int streamInqFileID(int streamID)
{
  stream_t *streamptr;

  streamptr = ( stream_t *) reshGetVal ( streamID, &streamOps );

  return (streamptr->fileID);
}

/* not used anymore */
/*
void streamDefineTaxis(int streamID)
{
  stream_t *streamptr;

  streamptr = ( stream_t *) reshGetVal ( streamID, &streamOps );

  if ( streamptr->tsteps == NULL )
    {
      int varID, nvars;
      int vlistID;

      vlistID = streamptr->vlistID;

      nvars = vlistNvars(vlistID);
      for ( varID = 0; varID < nvars; varID++ )
	if ( vlistInqVarTsteptype(vlistID, varID) == TSTEP_CONSTANT ) break;

      if ( varID == nvars )
	{
	  int taxisID;

	  taxisID = vlistInqTaxis(vlistID);
	  if ( taxisID == CDI_UNDEFID )
	    {
	      taxisID = taxisCreate(TAXIS_ABSOLUTE);
	      vlistDefTaxis(vlistID, taxisID);
	    }

	  (void) streamDefTimestep(streamID, 0);
	}
      else
	Error("time axis undefined");
    }
}
*/

void streamDefDimgroupID(int streamID, int dimgroupID)
{
  stream_t *streamptr;

  streamptr = ( stream_t *) reshGetVal ( streamID, &streamOps );

  streamptr->dimgroupID = dimgroupID;
}


int streamInqDimgroupID(int streamID)
{
  stream_t *streamptr;

  streamptr = ( stream_t *) reshGetVal ( streamID, &streamOps );

  return (streamptr->dimgroupID);
}


void cdiDefAccesstype(int streamID, int type)
{
  stream_t *streamptr;

  streamptr = ( stream_t *) reshGetVal ( streamID, &streamOps );

  if ( streamptr->accesstype == UNDEFID )
    {
      streamptr->accesstype = type;
    }
  else
    {
      if ( streamptr->accesstype != type )
	{
	  if ( streamptr->accesstype == TYPE_REC )
	    Error("Changing access type from REC to VAR not allowed!");
	  else
	    Error("Changing access type from VAR to REC not allowed!");
	}
    }
}


int cdiInqAccesstype(int streamID)
{
  stream_t *streamptr;

  streamptr = ( stream_t *) reshGetVal ( streamID, &streamOps );

  return (streamptr->accesstype);
} 


int streamInqNvars ( int streamID )
{
  stream_t * streamptr;
  streamptr = ( stream_t * ) reshGetVal ( streamID, &streamOps );
  return ( streamptr->nvars );
}


int  streamCompareP ( void * streamptr1, void * streamptr2 )
{
  stream_t * s1 = ( stream_t * ) streamptr1;
  stream_t * s2 = ( stream_t * ) streamptr2;
  int differ = -1;
  int equal  = 0;
  int len;

  xassert ( s1 );
  xassert ( s2 );

  if ( s1->filetype  != s2->filetype  ) return differ;
  if (  namespaceAdaptKey2 ( s1->vlistIDorig ) !=
	namespaceAdaptKey2 ( s2->vlistIDorig )) return differ;
  if ( s1->byteorder != s2->byteorder ) return differ;
  if ( s1->comptype  != s2->comptype  ) return differ;
  if ( s1->complevel != s2->complevel ) return differ;

  if ( s1->filename )
    {
      len = strlen ( s1->filename ) + 1;
      if ( memcmp ( s1->filename, s2->filename, len ))
	return differ;
    }
  else if ( s2->filename )
    return differ;

  return equal;
}


void streamDestroyP ( void * streamptr )
{
  int id;
  stream_t * sp = ( stream_t * ) streamptr;

  xassert ( sp );

  id = sp->self;
  streamClose ( id );
}


void streamPrintP   ( void * streamptr, FILE * fp )
{
  stream_t * sp = ( stream_t * ) streamptr;

  if ( !sp ) return;

  fprintf ( fp, "#\n");
  fprintf ( fp, "# streamID %d\n", sp->self);
  fprintf ( fp, "#\n");
  fprintf ( fp, "self          = %d\n", sp->self );
  fprintf ( fp, "accesstype    = %d\n", sp->accesstype );
  fprintf ( fp, "accessmode    = %d\n", sp->accessmode );
  fprintf ( fp, "filetype      = %d\n", sp->filetype );
  fprintf ( fp, "byteorder     = %d\n", sp->byteorder );
  fprintf ( fp, "fileID        = %d\n", sp->fileID );
  fprintf ( fp, "dimgroupID    = %d\n", sp->dimgroupID );
  fprintf ( fp, "filemode      = %d\n", sp->filemode );
  fprintf ( fp, "//off_t numvals;\n" );
  fprintf ( fp, "filename      = %s\n", sp->filename );
  fprintf ( fp, "//Record   *record;\n" );
  fprintf ( fp, "nrecs         = %d\n", sp->nrecs );
  fprintf ( fp, "nvars         = %d\n", sp->nvars );
  fprintf ( fp, "varlocked     = %d\n", sp->varlocked );
  fprintf ( fp, "//svarinfo_t *vars;\n" );
  fprintf ( fp, "varsAllocated = %d\n", sp->varsAllocated );
  fprintf ( fp, "varinit       = %d\n", sp->varinit );
  fprintf ( fp, "curTsID       = %d\n", sp->curTsID );
  fprintf ( fp, "rtsteps       = %d\n", sp->rtsteps );
  fprintf ( fp, "//long ntsteps;\n" );
  fprintf ( fp, "numTimestep   = %d\n", sp->numTimestep );
  fprintf ( fp, "//  tsteps_t   *tsteps;\n" );
  fprintf ( fp, "tstepsTableSize= %d\n", sp->tstepsTableSize );
  fprintf ( fp, "tstepsNextID  = %d\n", sp->tstepsNextID );
  fprintf ( fp, "//basetime_t  basetime;\n" );
  fprintf ( fp, "ncmode        = %d\n", sp->ncmode );
  fprintf ( fp, "vlistID       = %d\n", sp->vlistID );
  fprintf ( fp, "//  int       xdimID[MAX_GRIDS_PS];\n" );
  fprintf ( fp, "//  int       ydimID[MAX_GRIDS_PS];\n" );
  fprintf ( fp, "//  int       zaxisID[MAX_ZAXES_PS];\n" );
  fprintf ( fp, "//  int       ncxvarID[MAX_GRIDS_PS];\n" );
  fprintf ( fp, "//  int       ncyvarID[MAX_GRIDS_PS];\n" );
  fprintf ( fp, "//  int       ncavarID[MAX_GRIDS_PS];\n" );
  fprintf ( fp, "historyID     = %d\n", sp->historyID );
  fprintf ( fp, "globalatts    = %d\n", sp->globalatts );
  fprintf ( fp, "localatts     = %d\n", sp->localatts );
  fprintf ( fp, "//  VCT       vct;\n" );
  fprintf ( fp, "unreduced     = %d\n", sp->unreduced );
  fprintf ( fp, "sortname      = %d\n", sp->sortname );
  fprintf ( fp, "have_missval  = %d\n", sp->have_missval );
  fprintf ( fp, "ztype         = %d\n", sp->comptype );
  fprintf ( fp, "zlevel        = %d\n", sp->complevel );
  fprintf ( fp, "curfile       = %d\n", sp->curfile );
  fprintf ( fp, "nfiles        = %d\n", sp->nfiles );
  fprintf ( fp, "//  char    **fnames;\n" );
  fprintf ( fp, "//  void    **gribContainers;\n" );
  fprintf ( fp, "vlistIDorig   = %d\n", sp->vlistIDorig );
}


void streamGetIndexList ( int nstreams, int * streamIndexList )
{
  reshGetResHListOfType ( nstreams, streamIndexList, &streamOps );
}

void
cdiStreamSetupVlist(stream_t *streamptr, int vlistID, int vlistIDorig)
{
  int nvars = vlistNvars(vlistID);
  streamptr->vlistID = vlistID;
  streamptr->vlistIDorig = vlistIDorig;
  for (int varID = 0; varID < nvars; varID++ )
    {
      int gridID  = vlistInqVarGrid(vlistID, varID);
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      stream_new_var(streamptr, gridID, zaxisID);
      if ( streamptr->have_missval )
        vlistDefVarMissval(vlistID, varID,
                           vlistInqVarMissval(vlistID, varID));
    }

  if (streamptr->filemode == 'w' )
    {
      if ( streamptr->filetype == FILETYPE_NC  ||
           streamptr->filetype == FILETYPE_NC2 ||
           streamptr->filetype == FILETYPE_NC4 ||
           streamptr->filetype == FILETYPE_NC4C )
        {
#ifdef HAVE_LIBNETCDF
          void (*myCdfDefVars)(stream_t *streamptr)
            = (void (*)(stream_t *))
            namespaceSwitchGet(NSSWITCH_CDF_STREAM_SETUP).func;
          myCdfDefVars(streamptr);
#endif
        }
      else if ( streamptr->filetype == FILETYPE_GRB  ||
                streamptr->filetype == FILETYPE_GRB2 )
        {
          gribContainersNew(streamptr);
        }
    }
}


static int
streamTxCode(void)
{
  return STREAM;
}


int streamNint = 11 ;


static int
streamGetPackSize(void * voidP, void *context)
{
  stream_t * streamP = ( stream_t * ) voidP;
  int packBufferSize
    = serializeGetSize(streamNint, DATATYPE_INT, context)
    + serializeGetSize(2, DATATYPE_UINT32, context)
    + serializeGetSize((int)strlen(streamP->filename) + 1,
                       DATATYPE_TXT, context)
    + serializeGetSize(1, DATATYPE_FLT64, context);
  return packBufferSize;
}


static void
streamPack(void * streamptr, void * packBuffer, int packBufferSize,
           int * packBufferPos, void *context)
{
  stream_t * streamP = ( stream_t * ) streamptr;
  int intBuffer[streamNint];

  intBuffer[0]  = streamP->self;
  intBuffer[1]  = streamP->filetype;
  intBuffer[2]  = (int)strlen(streamP->filename) + 1;
  intBuffer[3]  = streamP->vlistID;
  intBuffer[4]  = streamP->vlistIDorig;
  intBuffer[5]  = streamP->byteorder;
  intBuffer[6]  = streamP->comptype;
  intBuffer[7]  = streamP->complevel;
  intBuffer[8]  = cdiDataUnreduced;
  intBuffer[9]  = cdiSortName;
  intBuffer[10] = cdiHaveMissval;

  serializePack(intBuffer, streamNint, DATATYPE_INT, packBuffer, packBufferSize, packBufferPos, context);
  uint32_t d = cdiCheckSum(DATATYPE_INT, streamNint, intBuffer);
  serializePack(&d, 1, DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);

  serializePack(&cdiDefaultMissval, 1, DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
  serializePack(streamP->filename, intBuffer[2], DATATYPE_TXT, packBuffer, packBufferSize, packBufferPos, context);
  d = cdiCheckSum(DATATYPE_TXT, intBuffer[2], streamP->filename);
  serializePack(&d, 1, DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
}

struct streamAssoc
streamUnpack(char * unpackBuffer, int unpackBufferSize,
             int * unpackBufferPos, int originNamespace, void *context)
{
  int intBuffer[streamNint], streamID;
  uint32_t d;
  char filename[CDI_MAX_NAME];

  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                  intBuffer, streamNint, DATATYPE_INT, context);
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                  &d, 1, DATATYPE_UINT32, context);
  xassert(cdiCheckSum(DATATYPE_INT, streamNint, intBuffer) == d);

  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                  &cdiDefaultMissval, 1, DATATYPE_FLT64, context);
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                  &filename, intBuffer[2], DATATYPE_TXT, context);
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                  &d, 1, DATATYPE_UINT32, context);
  xassert(d == cdiCheckSum(DATATYPE_TXT, intBuffer[2], filename));
  streamID = streamOpenWrite ( filename, intBuffer[1] );
  xassert ( streamID >= 0 &&
            namespaceAdaptKey ( intBuffer[0], originNamespace ) == streamID );
  streamDefByteorder(streamID, intBuffer[5]);
  streamDefCompType(streamID, intBuffer[6]);
  streamDefCompLevel(streamID, intBuffer[7]);
  cdiDefGlobal("REGULARGRID", intBuffer[8]);
  cdiDefGlobal("SORTNAME", intBuffer[9]);
  cdiDefGlobal("HAVE_MISSVAL", intBuffer[10]);
  struct streamAssoc retval = { streamID, intBuffer[3], intBuffer[4] };
  return retval;
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
