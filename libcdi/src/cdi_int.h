#ifndef _CDI_INT_H
#define _CDI_INT_H

#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <sys/types.h>

#ifndef strdupx
#ifndef strdup
char *strdup(const char *s);
#endif
#define strdupx  strdup
/*
#define strdupx(s)			          \
({					      	  \
   const char *__old = (s);			  \
   size_t __len = strlen(__old) + 1;		  \
   char *__new = (char *) malloc(__len);	  \
   (char *) memcpy(__new, __old, __len);	  \
})
*/
#endif

#ifndef  M_PI
#define  M_PI        3.14159265358979323846  /* pi */
#endif


#ifndef  _ERROR_H
#  include "error.h"
#endif
#ifndef _BASETIME_H
#  include "basetime.h"
#endif
#ifndef _TIMEBASE_H
#  include "timebase.h"
#endif
#ifndef  _TAXIS_H
#  include "taxis.h"
#endif
#ifndef  _CDI_LIMITS_H
#  include "cdi_limits.h"
#endif
#ifndef  _SERVICE_H
#  include "service.h"
#endif
#ifndef  _EXTRA_H
#  include "extra.h"
#endif
#ifndef  _IEG_H
#  include "ieg.h"
#endif
#ifndef RESOURCE_HANDLE_H
#  include "resource_handle.h"
#endif


#define check_parg(arg)  if ( arg == 0 ) Warning("Argument '" #arg "' not allocated!")

#if defined (__xlC__) /* performance problems on IBM */
#ifndef DBL_IS_NAN
#  define DBL_IS_NAN(x)     ((x) != (x))
#endif
#else
#ifndef DBL_IS_NAN
#if  defined  (HAVE_DECL_ISNAN)
#  define DBL_IS_NAN(x)     (isnan(x))
#elif  defined  (FP_NAN)
#  define DBL_IS_NAN(x)     (fpclassify(x) == FP_NAN)
#else
#  define DBL_IS_NAN(x)     ((x) != (x))
#endif
#endif
#endif

#ifndef DBL_IS_EQUAL
/*#define DBL_IS_EQUAL(x,y) (!(x < y || y < x)) */
#  define DBL_IS_EQUAL(x,y) (DBL_IS_NAN(x)||DBL_IS_NAN(y)?(DBL_IS_NAN(x)&&DBL_IS_NAN(y)?1:0):!(x < y || y < x))
#endif

#ifndef IS_EQUAL
#  define IS_NOT_EQUAL(x,y) (x < y || y < x)
#  define IS_EQUAL(x,y)     (!IS_NOT_EQUAL(x,y))
#endif


#define  FALSE  0
#define  TRUE   1

#define  TYPE_REC  0
#define  TYPE_VAR  1

#define  MEMTYPE_DOUBLE  1
#define  MEMTYPE_FLOAT   2

typedef struct
{
  void     *buffer;
  size_t    buffersize;
  off_t     position;
  int       recsize;
  int       size;
  int       dataread;
  int       param;
  int       level;
  int       date;
  int       time;
  int       gridID;
  int       zaxisID;
  int       used;
  int       nrec;
  int       varID;
  int       levelID;
  int       recid;
  int       prec;
  int       sec0[2];
  int       sec1[1024];
  int       sec2[4096];
  int       sec3[2];
  int       sec4[512];
  void     *exsep;
}
Record;


typedef struct
{
  off_t     position;
  size_t    size;
  int       zip;
  int       param;
  int       ilevel;
  int       ilevel2;
  int       ltype;
  short     used;
  short     varID;
  short     levelID;
  char      varname[32]; /* needed for grib decoding with GRIB_API */
}
record_t;


typedef struct {
  record_t *records;
  int       recordSize;  /* number of allocated records           */
  int      *recIDs;      /* IDs of non constant records           */
  int       nrecs;       /* number of used records                */
                         /* tsID=0 nallrecs                       */
                         /* tsID>0 number of non constant records */
  int       nallrecs;    /* number of all records                 */
  int       curRecID;    /* current record ID                     */
  long      next;
  off_t     position;    /* timestep file position                */
  taxis_t   taxis;
}
tsteps_t;


typedef struct {
  int       ncvarid;
  int       nlevs;
  int      *level;       /* record IDs */
  int      *lindex;      /* level index */
  int       defmiss;     /* TRUE if missval is defined in file */

  int       isUsed;
  int       gridID;
  int       zaxisID;
  int       tsteptype;   /* TSTEP_* */
}
svarinfo_t;


typedef struct {
  int       ilev;
  int       mlev;
  int       ilevID;
  int       mlevID;
}
VCT;


typedef struct {
  int         self;
  int         accesstype;   /* TYPE_REC or TYPE_VAR */
  int         accessmode;
  int         filetype;
  int         byteorder;
  int         fileID;
  int         dimgroupID;
  int         filemode;
  off_t       numvals;
  char       *filename;
  Record     *record;
  int        nrecs;        /* number of records                  */
  int         nvars;        /* number of variables                */
  int         varlocked;    /* variables locked                   */
  svarinfo_t *vars;
  int         varsAllocated;
  int         varinit;
  int         curTsID;      /* current timestep ID */
  int         rtsteps;      /* number of tsteps accessed       */
  long        ntsteps;      /* number of tsteps : only set if all records accessed */
  int         numTimestep;  /* number of tsteps : only set if all records accessed */
  tsteps_t   *tsteps;
  int         tstepsTableSize;
  int         tstepsNextID;
  basetime_t  basetime;
  int         ncmode;
  int         vlistID;
  int         xdimID[MAX_GRIDS_PS];
  int         ydimID[MAX_GRIDS_PS];
  int         zaxisID[MAX_ZAXES_PS];
  int         ncxvarID[MAX_GRIDS_PS];
  int         ncyvarID[MAX_GRIDS_PS];
  int         ncavarID[MAX_GRIDS_PS];
  int         historyID;
  int         globalatts;
  int         localatts;
  VCT         vct;
  int         unreduced;
  int         sortname;
  int         have_missval;
  int         comptype;      // compression type
  int         complevel;     // compression level
  int         curfile;
  int         nfiles;
  char      **fnames;
#if defined (GRIBCONTAINER2D)
  void      **gribContainers;
#else
  void       *gribContainers;
#endif
  int         vlistIDorig;
  /* only used by MPI-parallelized version of library */
  int       ownerRank;    // MPI rank of owner process
  /* ---------------------------------- */
  /* Local change: 2013-02-18, FP (DWD) */
  /* ---------------------------------- */

  void *gh; // grib handle
}
stream_t;


extern int CDI_Debug;      /* If set to 1, debuggig (default 0)            */
extern int cdiGribApiDebug;
extern double cdiDefaultMissval;
extern int cdiDefaultInstID;
extern int cdiDefaultModelID;
extern int cdiDefaultTableID;
extern int cdiDefaultLeveltype;
//extern int cdiNcMissingValue;
extern int cdiNcChunksizehint;
extern int cdiChunkType;
extern int cdiSplitLtype105;

extern char *cdiPartabPath;
extern int   cdiPartabIntern;

stream_t *stream_to_pointer(int idx);
stream_t *stream_new_entry(void);
void stream_delete_entry(stream_t *streamptr);
void stream_check_ptr(const char *caller, stream_t *streamptr);

int     streamInqFileID(int streamID);

void    gridDefHasDims(int gridID, int hasdims);
int     gridInqHasDims(int gridID);
const char *gridNamePtr(int gridtype);
char   *zaxisNamePtr(int leveltype);
int     zaxisInqLevelID(int zaxisID, double level);

void    streamCheckID(const char *caller, int streamID);

void    streamDefineTaxis(int streamID);

int     streamsNewEntry(int filetype);
void    streamsInitEntry(int streamID);
void    cdiStreamSetupVlist(stream_t *streamptr, int vlistID, int vlistIDorig);
int     stream_new_var(stream_t *streamptr, int gridID, int zaxisID);

int     tstepsNewEntry(stream_t *streamptr);

const char *strfiletype(int filetype);

void    cdi_generate_vars(stream_t *streamptr);

void    vlist_check_contents(int vlistID);

void    cdi_create_records(stream_t *streamptr, int tsID);

int     recordNewEntry(stream_t *streamptr, int tsID);

void    cdiCreateTimesteps(stream_t *streamptr);

void    recordInitEntry(record_t *record);

void    cdiCheckZaxis(int zaxisID);

void    cdiPrintDatatypes(void);

void    cdiDefAccesstype(int streamID, int type);
int     cdiInqAccesstype(int streamID);

void    streamDefDimgroupID(int streamID, int dimgroupID);
int     streamInqDimgroupID(int streamID);

int     getByteswap(int byteorder);

int     streamSize ();
void    streamGetIndexList ( int, int * );


void  cdiInitialize(void);

void stream_write_record(int streamID, int memtype, const void *data, int nmiss);

void uuid2str(const char *uuid, char *uuidstr);
void str2uuid(const char *uuidstr, char *uuid);

#define CDI_UNIT_PA   1
#define CDI_UNIT_HPA  2
#define CDI_UNIT_MM   3
#define CDI_UNIT_CM   4
#define CDI_UNIT_DM   5
#define CDI_UNIT_M    6

struct streamAssoc
{
  int streamID, vlistID, vlistIDorig;
};

struct streamAssoc
streamUnpack(char * unpackBuffer, int unpackBufferSize,
             int * unpackBufferPos, int originNamespace, void *context);

int
cdiStreamOpenDefaultDelegate(const char *filename, const char *filemode,
                             int filetype, stream_t *streamptr,
                             int recordBufIsToBeCreated);

void
cdiStreamDefVlist_(int streamID, int vlistID);
void
cdiStreamWriteVar_(int streamID, int varID, int memtype, const void *data,
                   int nmiss);
void
cdiStreamwriteVarChunk_(int streamID, int varID, int memtype,
                        const int rect[][2], const void *data, int nmiss);
void
cdiStreamCloseDefaultDelegate(stream_t *streamptr,
                              int recordBufIsToBeDeleted);

int cdiStreamDefTimestep_(stream_t *streamptr, int tsID);

void cdiStreamSync_(stream_t *streamptr);

char *cdiUnitNamePtr(int cdi_unit);

extern const resOps streamOps;

#endif  /* _CDI_INT_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
