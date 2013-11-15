#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include "cdi.h"
#include "namespace.h"
#include "resource_handle.h"
#include "pio_util.h"
#include "serialize.h"
#include "error.h"
#include "cdf_int.h"
#include "file.h"
#include "cdi_int.h"
#include "stream_cdf.h"

static int nNamespaces = 1;
static int activeNamespace = 0;

#ifdef HAVE_LIBNETCDF
#define CDI_NETCDF_SWITCHES                     \
  { .func = (void (*)()) nc__create },          \
  { .func = (void (*)()) cdf_def_var_serial },  \
  { .func = (void (*)()) cdfDefTimestep },      \
  { .func = (void (*)()) cdfDefVars }

#else
#define CDI_NETCDF_SWITCHES
#endif

#define defaultSwitches {                                   \
    { .func = (void (*)()) cdiAbortC_serial },              \
    { .func = (void (*)()) serializeGetSizeInCore },        \
    { .func = (void (*)()) serializePackInCore },           \
    { .func = (void (*)()) serializeUnpackInCore },         \
    { .func = (void (*)()) fileOpen_serial },               \
    { .func = (void (*)()) fileWrite },                     \
    { .func = (void (*)()) fileClose_serial },              \
    { .func = (void (*)()) cdiStreamOpenDefaultDelegate },  \
    { .func = (void (*)()) cdiStreamDefVlist_ },            \
    { .func = (void (*)()) cdiStreamWriteVar_ },            \
    { .func = (void (*)()) cdiStreamwriteVarChunk_ },       \
    { .data = NULL },                                       \
    { .func = (void (*)()) cdiStreamCloseDefaultDelegate }, \
    { .func = (void (*)()) cdiStreamDefTimestep_ }, \
    { .func = (void (*)()) cdiStreamSync_ },                \
    CDI_NETCDF_SWITCHES                        \
    }

struct namespace
{
  statusCode resStage;
  union namespaceSwitchValue switches[NUM_NAMESPACE_SWITCH];
} initialNamespace = {
  .resStage = STAGE_DEFINITION,
  .switches = defaultSwitches
};

struct namespace *namespaces = &initialNamespace;

static int namespacesSize = 1;

#if  defined  (HAVE_LIBPTHREAD)
#  include <pthread.h>

static pthread_mutex_t namespaceMutex = PTHREAD_MUTEX_INITIALIZER;

#  define NAMESPACE_LOCK()         pthread_mutex_lock(&namespaceMutex)
#  define NAMESPACE_UNLOCK()       pthread_mutex_unlock(&namespaceMutex)

#else

#  define NAMESPACE_LOCK()
#  define NAMESPACE_UNLOCK()

#endif


enum {
  intbits = sizeof(int) * CHAR_BIT,
  nspbits = 4,
  idxbits = intbits - nspbits,
  nspmask = (( 1 << nspbits ) - 1) << idxbits,
  idxmask = ( 1 << idxbits ) - 1,
};

enum {
  NUM_NAMESPACES = 1 << nspbits,
  NUM_IDX = 1 << idxbits,
};


#if 0
void namespaceShowbits ( int n, char *name )
{
  int i;
  unsigned mask;
  char bitvalues[intbits + 1];

  mask = 1;
  for ( i = 0; i < intbits; i++ )
    {
      bitvalues[i] = ((unsigned)n & mask) ? '1':'0';
      mask <<= 1;
    }
  bitvalues[intbits] = '\0';
  fprintf (stdout, "%s: %s\n", name, bitvalues );
}
#endif


int namespaceIdxEncode ( namespaceTuple_t tin )
{
  xassert ( tin.nsp < NUM_NAMESPACES && tin.idx < NUM_IDX);
  return ( tin.nsp << idxbits ) + tin.idx;
}

int namespaceIdxEncode2 ( int nsp, int idx )
{
  xassert(nsp < NUM_NAMESPACES && idx < NUM_IDX);
  return ( nsp << idxbits ) + idx;
}


namespaceTuple_t namespaceResHDecode ( int resH )
{
  namespaceTuple_t tin;

  tin.idx = resH & idxmask;
  tin.nsp = (int)(((unsigned)( resH & nspmask )) >> idxbits);

  return tin;
}

int
namespaceNew()
{
  int newNamespaceID = -1;
  NAMESPACE_LOCK();
  if (namespacesSize > nNamespaces)
    {
      /* namespace is already available and only needs reinitialization */
      for (int i = 0; i < namespacesSize; ++i)
        if (namespaces[i].resStage == STAGE_UNUSED)
          {
            newNamespaceID = i;
            break;
          }
    }
  else if (namespacesSize == 1)
    {
      /* make room for additional namespace */
      struct namespace *newNameSpaces
        = xmalloc((namespacesSize + 1) * sizeof (namespaces[0]));
      memcpy(newNameSpaces, namespaces, sizeof (namespaces[0]));
      namespaces = newNameSpaces;
      ++namespacesSize;
      newNamespaceID = 1;
    }
  else if (namespacesSize < NUM_NAMESPACES)
    {
      /* make room for additional namespace */
      newNamespaceID = namespacesSize;
      namespaces
        = xrealloc(namespaces, (namespacesSize + 1) * sizeof (namespaces[0]));
      ++namespacesSize;
    }
  else /* implicit: namespacesSize >= NUM_NAMESPACES */
    {
      NAMESPACE_UNLOCK();
      return -1;
    }
  xassert(newNamespaceID >= 0 && newNamespaceID < NUM_NAMESPACES);
  ++nNamespaces;
  namespaces[newNamespaceID].resStage = STAGE_DEFINITION;
  memcpy(namespaces[newNamespaceID].switches,
         (union namespaceSwitchValue[NUM_NAMESPACE_SWITCH])defaultSwitches,
         sizeof (namespaces[newNamespaceID].switches));
  reshListCreate(newNamespaceID);
  NAMESPACE_UNLOCK();
  return newNamespaceID;
}

void
namespaceDelete(int namespaceID)
{
  NAMESPACE_LOCK();
  xassert(namespaceID < namespacesSize && nNamespaces);
  reshListDestruct(namespaceID);
  namespaces[namespaceID].resStage = STAGE_UNUSED;
  --nNamespaces;
  NAMESPACE_UNLOCK();
}

void namespaceCleanup ( void )
{
  if ( nNamespaces > 1 )
    {
      initialNamespace = namespaces[0];
      free(namespaces);
      namespaces = &initialNamespace;
      nNamespaces = 1;
    }
}


int namespaceGetNumber ()
{
  return nNamespaces;
}


void pioNamespaceSetActive ( int nId )
{
  xassert(nId < namespacesSize && nId >= 0
          && namespaces[nId].resStage != STAGE_UNUSED);
  activeNamespace = nId;
}


int namespaceGetActive ()
{
  return activeNamespace;
}

int namespaceAdaptKey ( int key, int nspTarget )
{
  namespaceTuple_t tin;
  int nsp;

  if ( key == CDI_UNDEFID ) return CDI_UNDEFID;

  tin.idx = key & idxmask;
  tin.nsp = (int)(((unsigned)( key & nspmask )) >> idxbits);

  xassert ( tin.nsp == nspTarget );

  nsp = namespaceGetActive ();

  return namespaceIdxEncode2 ( nsp, tin.idx );
}


int namespaceAdaptKey2 ( int key )
{
  namespaceTuple_t tin;
  int nsp;

  if ( key == CDI_UNDEFID ) return CDI_UNDEFID;

  tin.idx = key & idxmask;
  tin.nsp = (int)(((unsigned)( key & nspmask )) >> idxbits);

  nsp = namespaceGetActive ();

  return namespaceIdxEncode2 ( nsp, tin.idx );
}


void namespaceDefResStatus ( statusCode argResStatus )
{
  int nsp = namespaceGetActive ();
  namespaces[nsp].resStage = argResStatus;
}


statusCode namespaceInqResStatus ( void )
{
  int nsp = namespaceGetActive ();
  return namespaces[nsp].resStage;
}

void namespaceSwitchSet(enum namespaceSwitch sw, union namespaceSwitchValue value)
{
  xassert(sw > NSSWITCH_NO_SUCH_SWITCH && sw < NUM_NAMESPACE_SWITCH);
  int nsp = namespaceGetActive();
  namespaces[nsp].switches[sw] = value;
}

union namespaceSwitchValue namespaceSwitchGet(enum namespaceSwitch sw)
{
  xassert(sw > NSSWITCH_NO_SUCH_SWITCH && sw < NUM_NAMESPACE_SWITCH);
  int nsp = namespaceGetActive();
  return namespaces[nsp].switches[sw];
}

void cdiReset(void)
{
  NAMESPACE_LOCK();
  for (int namespaceID = 0; namespaceID < namespacesSize; ++namespaceID)
    namespaceDelete(namespaceID);
  namespaces = &initialNamespace;
  namespacesSize = 1;
  nNamespaces = 1;
  activeNamespace = 0;
  NAMESPACE_UNLOCK();
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
