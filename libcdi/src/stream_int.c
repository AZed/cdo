#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdarg.h>
#include <string.h> 
#include <errno.h>
#include <math.h>
#include <ctype.h>

#include "dmemory.h"

#include "cdi.h"
#include "stream_int.h"

#if  defined  (HAVE_LIBCGRIBEX)
#include "cgribex.h"
#endif

int cdiDefaultCalendar = CALENDAR_PROLEPTIC;

int cdiDefaultInstID   = CDI_UNDEFID;
int cdiDefaultModelID  = CDI_UNDEFID;
int cdiDefaultTableID  = CDI_UNDEFID;
int cdiNcMissingValue  = CDI_UNDEFID;
int cdiSplitLtype105   = CDI_UNDEFID;

int cdiIgnoreAttCoordinates = FALSE;
int cdiSkipRecords          = 0;
int cdiInventoryMode        = 1;

char *cdiPartabPath   = NULL;
int   cdiPartabIntern = 1;

double cdiDefaultMissval = -9.E33;


char *Filetypes[] = {
  "UNKNOWN",
  "GRIB",
  "GRIB2",
  "netCDF",
  "netCDF2",
  "netCDF4",
  "SERVICE",
  "EXTRA",
  "IEG",
  "HDF5",
};

#undef  UNDEFID
#define UNDEFID  CDI_UNDEFID


int CDI_Debug   = 0;    /* If set to 1, debugging           */


int cdiDefaultLeveltype = -1;
static int cdiDataUnreduced = 0;
static int cdiSortName = 0;
static int cdiHaveMissval = 0;


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
		}
	      break;
	    }
	}

      if ( fact ) envValue = fact*atol(envString);

      if ( CDI_Debug ) Message("set %s to %ld", envName, envValue);
    }

  return (envValue);
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

      value = cdiGetenvInt("CD_REGULARGRID");
      if ( value >= 0 ) cdiDataUnreduced = (int) value;

      value = cdiGetenvInt("CDI_REGULARGRID");
      if ( value >= 0 ) cdiDataUnreduced = (int) value;

      value = cdiGetenvInt("CDI_SORTNAME");
      if ( value >= 0 ) cdiSortName = (int) value;

      value = cdiGetenvInt("CDI_HAVE_MISSVAL");
      if ( value >= 0 ) cdiHaveMissval = (int) value;

      value = cdiGetenvInt("CD_LEVELTYPE");
      if ( value >= 0 ) cdiDefaultLeveltype = (int) value;

      value = cdiGetenvInt("CDI_LEVELTYPE");
      if ( value >= 0 ) cdiDefaultLeveltype = (int) value;

      envString = getenv("CD_MISSVAL");
      if ( envString ) cdiDefaultMissval = atof(envString);

      envString = getenv("CDI_MISSVAL");
      if ( envString ) cdiDefaultMissval = atof(envString);

      envString = getenv("NC_MISSING_VALUE");
      if ( envString ) cdiNcMissingValue = atoi(envString);

      envString = getenv("SPLIT_LTYPE_105");
      if ( envString ) cdiSplitLtype105 = atoi(envString);

      envString = getenv("IGNORE_ATT_COORDINATES");
      if ( envString ) cdiIgnoreAttCoordinates = atoi(envString);

      envString = getenv("CDI_SKIP_RECORDS");
      if ( envString )
	{
	  cdiSkipRecords = atoi(envString);
	  cdiSkipRecords = cdiSkipRecords > 0 ? cdiSkipRecords : 0;
	}

      envString = getenv("GRIB_INVENTORY_MODE");
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
      gribSetCalendar(cdiDefaultCalendar);

      envString = getenv("PARTAB_INTERN");
      if ( envString ) cdiPartabIntern = atoi(envString);

      envString = getenv("PARTAB_PATH");
      if ( envString ) cdiPartabPath = strdup(envString);
    }
}


char *strfiletype(int filetype)
{
  char *name;
  int size = (int) (sizeof(Filetypes)/sizeof(char *));

  if ( filetype > 0 && filetype < size )
    name = Filetypes[filetype];
  else
    name = Filetypes[0];  

  return (name);
}


static int  STREAM_Debug = 0;   /* If set to 1, debugging */

static int _stream_min = MIN_STREAMS;
static int _stream_max = MAX_STREAMS;

static void stream_initialize(void);

static int _stream_init = FALSE;

#if  defined  (HAVE_LIBPTHREAD)
#  include <pthread.h>

static pthread_once_t  _stream_init_thread = PTHREAD_ONCE_INIT;
static pthread_mutex_t _stream_mutex;

#  define STREAM_LOCK()         pthread_mutex_lock(&_stream_mutex)
#  define STREAM_UNLOCK()       pthread_mutex_unlock(&_stream_mutex)
#  define STREAM_INIT()        \
   if ( _stream_init == FALSE ) pthread_once(&_stream_init_thread, stream_initialize)

#else

#  define STREAM_LOCK()
#  define STREAM_UNLOCK()
#  define STREAM_INIT()        \
   if ( _stream_init == FALSE ) stream_initialize()

#endif


typedef struct _streamPtrToIdx {
  int       idx;
  int       next;
  stream_t *ptr;
} streamPtrToIdx;


static streamPtrToIdx *_streamList  = NULL;
static int             _streamAvail = -1;

static
void stream_list_new(void)
{
  assert(_streamList == NULL);

  _streamList = (streamPtrToIdx *) malloc(_stream_min*sizeof(streamPtrToIdx));
}

static
void stream_list_delete(void)
{
  if ( _streamList ) free(_streamList);
}

static
void stream_init_pointer(void)
{
  int i;
  
  for ( i = 0; i < _stream_min; ++i )
    {
      _streamList[i].idx  = i;
      _streamList[i].next = i + 1;
      _streamList[i].ptr  = NULL;
    }

  _streamList[_stream_min-1].next = -1;

  _streamAvail = 0;
}

static
void stream_list_extend(void)
{
  int nstreams;
  int i;

  assert(_streamList != NULL);

  nstreams = _stream_min + MIN_STREAMS;

  if ( nstreams <= _stream_max)
    {
      _streamList = (streamPtrToIdx *) realloc(_streamList, nstreams*sizeof(streamPtrToIdx));
  
      for ( i = _stream_min; i < nstreams; ++i )
	{
	  _streamList[i].idx  = i;
	  _streamList[i].next = i + 1;
	  _streamList[i].ptr  = NULL;
	}

      _streamAvail = _stream_min;
      _streamList[_stream_min-1].next = _stream_min;
      _stream_min = nstreams;
      _streamList[_stream_min-1].next = -1;
    }
  else
    Warning("Too many open streams (limit is %d)!", _stream_max);
}


stream_t *stream_to_pointer(int idx)
{
  stream_t *streamptr = NULL;

  STREAM_INIT();

  if ( idx >= 0 && idx < _stream_min )
    {
      STREAM_LOCK();

      streamptr = _streamList[idx].ptr;

      STREAM_UNLOCK();
    }
  else
    Error("stream index %d undefined!", idx);

  return (streamptr);
}

/* Create an index from a pointer (add the pointer to the stream list) */
static
int stream_from_pointer(stream_t *ptr)
{
  int idx = -1;

  if ( ptr )
    {
      STREAM_LOCK();

      if ( _streamAvail < 0 ) stream_list_extend();

      if ( _streamAvail >= 0 )
	{
	  streamPtrToIdx *newptr;

	  newptr       = &_streamList[_streamAvail];
	  _streamAvail = newptr->next;
	  newptr->next = -1;
	  idx	       = newptr->idx;
	  newptr->ptr  = ptr;
      
	  if ( STREAM_Debug )
	    Message("Pointer %p has idx %d from stream list", ptr, idx);
	}

      STREAM_UNLOCK();
    }
  else
    Error("Internal problem (pointer %p undefined)", ptr);

  return (idx);
}

static
void stream_init_entry(stream_t *streamptr)
{
  int i;

  streamptr->self              = stream_from_pointer(streamptr);

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
}


stream_t *stream_new_entry(void)
{
  stream_t *streamptr;

  cdiInitialize(); /* ***************** make MT version !!! */

  STREAM_INIT();

  streamptr = (stream_t *) malloc(sizeof(stream_t));

  if ( streamptr ) stream_init_entry(streamptr);

  return (streamptr);
}


void stream_delete_entry(stream_t *streamptr)
{
  int idx;

  idx = streamptr->self;

  STREAM_LOCK();

  free(streamptr);

  _streamList[idx].next = _streamAvail;
  _streamList[idx].ptr  = 0;
  _streamAvail          = idx;

  STREAM_UNLOCK();

  if ( STREAM_Debug )
    Message("Removed idx %d from stream list", idx);
}


static
void stream_initialize(void)
{
  char *env;

#if  defined  (HAVE_LIBPTHREAD)
  /* initialize global API mutex lock */
  pthread_mutex_init(&_stream_mutex, NULL);
#endif

  env = getenv("STREAM_DEBUG");
  if ( env ) STREAM_Debug = atoi(env);

  stream_list_new();
  atexit(stream_list_delete);

  STREAM_LOCK();

  stream_init_pointer();

  STREAM_UNLOCK();

  _stream_init = TRUE;
}


void stream_check_ptr(const char *caller, stream_t *streamptr)
{
  if ( streamptr == NULL )
    Errorc("stream undefined!");
}


int streamSize(void)
{
  int streamsize = 0;
  int i;
  
  STREAM_INIT();

  STREAM_LOCK();

  for ( i = 0; i < _stream_min; i++ )
    if ( _streamList[i].ptr ) streamsize++;

  STREAM_UNLOCK();

  return (streamsize);
}


void cdiDefGlobal(const char *string, int val)
{
  if ( strcmp(string, "REGULARGRID") == 0 )
    {
      cdiDataUnreduced = val;
    }
  else if ( strcmp(string, "SORTNAME") == 0 )
    {
      cdiSortName = val;
    }
  else if ( strcmp(string, "HAVE_MISSVAL") == 0 )
    {
      cdiHaveMissval = val;
    }
  else
    {
      Warning("Unsupported global key: %s", string);
    }
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


void cdiCheckContents(int streamID)
{
  int index, nzaxis, zaxisID;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  vlistID = streamInqVlist(streamID);
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

  streamptr = stream_to_pointer(streamID);

  return (streamptr->fileID);
}


void streamDefineTaxis(int streamID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( streamptr->tsteps == NULL )
    {
      int varID, nvars;
      int vlistID;
  
      vlistID = streamInqVlist(streamID);

      nvars = vlistNvars(vlistID);
      for ( varID = 0; varID < nvars; varID++ )
	if ( vlistInqVarTime(vlistID, varID) == TIME_VARIABLE ) break;

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


void streamDefDimgroupID(int streamID, int dimgroupID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  streamptr->dimgroupID = dimgroupID;
}


int streamInqDimgroupID(int streamID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  return (streamptr->dimgroupID);
}


void cdiDefAccesstype(int streamID, int type)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

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

  streamptr = stream_to_pointer(streamID);

  return (streamptr->accesstype);
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
