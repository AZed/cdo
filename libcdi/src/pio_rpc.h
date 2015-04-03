#ifndef PIO_RPC_H
#define PIO_RPC_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef USE_MPI

#include <mpi.h>
#include <yaxt.h>

typedef enum
  {
    FINALIZE,
    RESOURCES,
    WINCREATE,
    WRITETS,
    COLLBUFTX,
    COLLBUFNMISS,
  } command;

#define MAXWINBUFFERSIZE ((size_t)2048 * 1024 * 1024)

enum
{
  numRPCFuncs = 4,
  STREAMOPEN = -1,
  STREAMDEFVLIST = -2,
  STREAMCLOSE = -3,
  STREAMDEFTIMESTEP = -4,
  HEADERSIZEMARKER = -numRPCFuncs - 1,
  PARTDESCMARKER = -numRPCFuncs - 2,
};
enum { MAXDATAFILENAME = 256, MINFUNCID = -numRPCFuncs, MAXFUNCID = -1 };
extern const char * const funcMap[numRPCFuncs];

struct headerSize
{
  int sizeID, numDataEntries, numRPCEntries;
};

struct dataRecord
{
  int streamID, varID, offset, nmiss;
};

struct funcCallDesc
{
  int funcID;
  union
  {
    struct
    {
      int streamID, vlistID;
    } streamChange;
    struct
    {
      int streamID, tsID, offset;
    } streamNewTimestep;
    struct
    {
      int fnamelen, offset, filetype;
    } newFile;
  } funcArgs;
};

/* Describes offset and ID of serialized partition descriptor.
 * partDescMarker == PARTDESCMARKER, always. */
struct partDescRecord
{
  int partDescMarker, offset;
  Xt_uid uid;
};

union winHeaderEntry
{
  struct headerSize headerSize;
  struct dataRecord dataRecord;
  struct funcCallDesc funcCall;
  struct partDescRecord partDesc;
};

/* round size to next multiple of factor */
static inline size_t
roundUpToMultiple(size_t size, size_t factor)
{
  return (size + factor - 1)/factor * factor;
}

enum
{
  /* align window base addresses and sizes to this value */
  PIO_WIN_ALIGN = sizeof (double),
};

#endif

#endif
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
