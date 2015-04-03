#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 /* PTHREAD_MUTEX_RECURSIVE */
#endif

#include <stdlib.h>
#include <stdio.h>

#include "dmemory.h"
#include "resource_handle.h"
#include "namespace.h"
#include "serialize.h"
#include "cdi.h"
#include "error.h"
#include "file.h"
#include "resource_unpack.h"
#include "institution.h"
#include "model.h"

enum { MIN_LIST_SIZE = 128 };

static void listInitialize(void);

typedef struct listElem {
  union
  {
    /* free-list management data */
    struct
    {
      int next, prev;
    } free;
    /* holding an actual value */
    struct
    {
      const resOps *ops;
      void         *val;//ptr
    } v;
  } res;
  int           status;
} listElem_t;

struct resHList_t
{
  int size, freeHead, hasDefaultRes;
  listElem_t *resources;
};

static struct resHList_t *resHList;

static int resHListSize = 0;

#if  defined  (HAVE_LIBPTHREAD)
#  include <pthread.h>

static pthread_once_t  listInitThread = PTHREAD_ONCE_INIT;
static pthread_mutex_t listMutex;

#  define LIST_LOCK()         pthread_mutex_lock(&listMutex)
#  define LIST_UNLOCK()       pthread_mutex_unlock(&listMutex)
#  define LIST_INIT(init0)         do {                         \
    pthread_once(&listInitThread, listInitialize);              \
    pthread_mutex_lock(&listMutex);                             \
    if ((init0) && (!resHList || !resHList[0].resources))       \
      reshListCreate(0);                                        \
    pthread_mutex_unlock(&listMutex);                           \
  } while (0)



#else

static int listInit = 0;

#  define LIST_LOCK()
#  define LIST_UNLOCK()
#  define LIST_INIT(init0)        do {                          \
  if ( !listInit )                                              \
    {                                                           \
      listInitialize();                                         \
      if ((init0) && (!resHList || !resHList[0].resources))     \
        reshListCreate(0);                                      \
      listInit = 1;                                             \
    }                                                           \
  } while(0)

#endif

/**************************************************************/

static void
listInitResources(int nsp)
{
  xassert(nsp < resHListSize && nsp >= 0);
  int size = resHList[nsp].size = MIN_LIST_SIZE;
  xassert(resHList[nsp].resources == NULL);
  resHList[nsp].resources = (listElem_t*) xcalloc(MIN_LIST_SIZE, sizeof(listElem_t));
  listElem_t *p = resHList[nsp].resources;

  for (int i = 0; i < size; i++ )
    {
      p[i].res.free.next = i + 1;
      p[i].res.free.prev = i - 1;
      p[i].status = RESH_UNUSED;
    }

  p[size-1].res.free.next = -1;
  resHList[nsp].freeHead = 0;
  int oldNsp = namespaceGetActive();
  namespaceSetActive(nsp);
  instituteDefaultEntries();
  modelDefaultEntries();
  namespaceSetActive(oldNsp);
}

static inline void
reshListClearEntry(int i)
{
  resHList[i].size = 0;
  resHList[i].resources = NULL;
  resHList[i].freeHead = -1;
}

void
reshListCreate(int namespaceID)
{
  LIST_INIT(namespaceID != 0);
  LIST_LOCK();
  if (resHListSize <= namespaceID)
    {
      resHList = (struct resHList_t*) xrealloc(resHList, (namespaceID + 1) * sizeof (resHList[0]));
      for (int i = resHListSize; i <= namespaceID; ++i)
        reshListClearEntry(i);
      resHListSize = namespaceID + 1;
    }
  listInitResources(namespaceID);
  LIST_UNLOCK();
}


/**************************************************************/

void
reshListDestruct(int namespaceID)
{
  LIST_LOCK();
  xassert(resHList && namespaceID >= 0 && namespaceID < resHListSize);
  int callerNamespaceID = namespaceGetActive();
  namespaceSetActive(namespaceID);
  if (resHList[namespaceID].resources)
    {
      for ( int j = 0; j < resHList[namespaceID].size; j++ )
        {
          listElem_t *listElem = resHList[namespaceID].resources + j;
          if (listElem->status != RESH_UNUSED)
            listElem->res.v.ops->valDestroy(listElem->res.v.val);
        }
      free(resHList[namespaceID].resources);
      reshListClearEntry(namespaceID);
    }
  if (resHList[callerNamespaceID].resources)
    namespaceSetActive(callerNamespaceID);
  LIST_UNLOCK();
}


static void listDestroy ( void )
{
  LIST_LOCK();
  for (int i = resHListSize; i > 0; --i)
    if (resHList[i-1].resources)
      namespaceDelete(i-1);
  free(resHList);
  resHList = NULL;
  cdiReset();
  LIST_UNLOCK();
}

/**************************************************************/

static
void listInitialize ( void )
{
#if  defined  (HAVE_LIBPTHREAD)
  pthread_mutexattr_t ma;
  pthread_mutexattr_init(&ma);
  pthread_mutexattr_settype(&ma, PTHREAD_MUTEX_RECURSIVE);
  /* initialize global API mutex lock */
  pthread_mutex_init ( &listMutex, &ma);
  pthread_mutexattr_destroy(&ma);
#endif
  /* file is special and has its own table, which needs to be
   * created, before we register the listDestroy exit handler */
  {
    int null_id;
    null_id = fileOpen_serial("/dev/null", "r");
    if (null_id != -1)
      fileClose_serial(null_id);
  }
  atexit ( listDestroy );
}

/**************************************************************/

static
void listSizeExtend()
{
  int nsp = namespaceGetActive ();
  int oldSize = resHList[nsp].size;
  int newListSize = oldSize + MIN_LIST_SIZE;

  resHList[nsp].resources = (listElem_t*) xrealloc(resHList[nsp].resources,
                                                   newListSize * sizeof(listElem_t));

  listElem_t *r = resHList[nsp].resources;
  for (int i = oldSize; i < newListSize; ++i)
    {
      r[i].res.free.next = i + 1;
      r[i].res.free.prev = i - 1;
      r[i].status = RESH_UNUSED;
    }

  if (resHList[nsp].freeHead != -1)
    r[resHList[nsp].freeHead].res.free.next
      = newListSize - 1;
  r[newListSize-1].res.free.next = resHList[nsp].freeHead;
  r[oldSize].res.free.prev = -1;
  resHList[nsp].freeHead = oldSize;
  resHList[nsp].size = newListSize;
}

/**************************************************************/

static void
reshPut_(int nsp, int entry, void *p, const resOps *ops)
{
  listElem_t *newListElem = resHList[nsp].resources + entry;
  int next = newListElem->res.free.next;
  if (next != -1)
    resHList[nsp].resources[next].res.free.prev = -1;
  resHList[nsp].freeHead = next;
  newListElem->res.v.val = p;
  newListElem->res.v.ops = ops;
  newListElem->status = RESH_ASSIGNED;
}

int reshPut ( void *p, const resOps *ops )
{
  xassert ( p && ops );

  LIST_INIT(1);

  LIST_LOCK();

  int nsp = namespaceGetActive ();

  if ( resHList[nsp].freeHead == -1) listSizeExtend();
  int entry = resHList[nsp].freeHead;
  cdiResH resH = namespaceIdxEncode2(nsp, entry);
  reshPut_(nsp, entry, p, ops);

  LIST_UNLOCK();

  return resH;
}

/**************************************************************/

static void
reshRemove_(int nsp, int idx)
{
  int curFree = resHList[nsp].freeHead;
  listElem_t *r = resHList[nsp].resources;
  r[idx].res.free.next = curFree;
  if (curFree != -1)
    r[curFree].res.free.prev = idx;
  r[idx].status = RESH_UNUSED;
  resHList[nsp].freeHead = idx;
}

void reshRemove ( cdiResH resH, const resOps * ops )
{
  int nsp;
  namespaceTuple_t nspT;

  LIST_INIT(1);

  LIST_LOCK();

  nsp = namespaceGetActive ();

  nspT = namespaceResHDecode ( resH );

  xassert ( nspT.nsp == nsp
            && nspT.idx >= 0
            && nspT.idx < resHList[nsp].size
            && resHList[nsp].resources[nspT.idx].status != RESH_UNUSED
            && resHList[nsp].resources[nspT.idx].res.v.ops
            && resHList[nsp].resources[nspT.idx].res.v.ops == ops );

  reshRemove_(nsp, nspT.idx);

  LIST_UNLOCK();
}

/**************************************************************/

void reshReplace(cdiResH resH, void *p, const resOps *ops)
{
  xassert(p && ops);
  LIST_INIT(1);
  LIST_LOCK();
  int nsp = namespaceGetActive();
  namespaceTuple_t nspT = namespaceResHDecode(resH);
  while (resHList[nsp].size <= nspT.idx)
    listSizeExtend();
  listElem_t *q = resHList[nsp].resources + nspT.idx;
  if (q->status != RESH_UNUSED)
    {
      q->res.v.ops->valDestroy(q->res.v.val);
      reshRemove_(nsp, nspT.idx);
    }
  reshPut_(nsp, nspT.idx, p, ops);
  LIST_UNLOCK();
}


static listElem_t *
reshGetElem(const char *caller, cdiResH resH, const resOps *ops)
{
  listElem_t *listElem;
  int nsp;
  namespaceTuple_t nspT;
  xassert ( ops );

  LIST_INIT(1);

  LIST_LOCK();

  nsp = namespaceGetActive ();

  nspT = namespaceResHDecode ( resH );

  if (nspT.nsp == nsp &&
      nspT.idx >= 0 &&
      nspT.idx < resHList[nsp].size)
    {
      listElem = resHList[nsp].resources + nspT.idx;
      LIST_UNLOCK();
    }
  else
    {
      LIST_UNLOCK();
      xabortC(caller, "Invalid namespace %d or index %d for resource handle %d when using namespace %d of size %d!",
              nspT.nsp, nspT.idx, (int)resH, nsp, resHList[nsp].size);
    }

  if ( !(listElem && listElem->res.v.ops == ops) )
    xabortC(caller, "Invalid resource handle %d, list element not found!",
            (int)resH);
  return listElem;
}


void *reshGetValue(const char * caller, cdiResH resH, const resOps * ops)
{
  return reshGetElem(caller, resH, ops)->res.v.val;
}

/**************************************************************/

void reshGetResHListOfType ( int c, int * resHs, const resOps * ops )
{
  int i, j = 0, nsp;

  xassert ( resHs && ops );

  LIST_INIT(1);

  LIST_LOCK();

  nsp = namespaceGetActive ();

  for ( i = 0; i < resHList[nsp].size && j < c; i++ )
    if (resHList[nsp].resources[i].status != RESH_UNUSED
        && resHList[nsp].resources[i].res.v.ops == ops)
      resHs[j++] = namespaceIdxEncode2(nsp, i);

  LIST_UNLOCK();
}

enum cdiApplyRet
cdiResHApply(enum cdiApplyRet (*func)(int id, void *res, const resOps *p,
                                      void *data), void *data)
{
  xassert(func);

  LIST_INIT(1);

  LIST_LOCK();

  int nsp = namespaceGetActive ();
  enum cdiApplyRet ret = CDI_APPLY_GO_ON;
  for (int i = 0; i < resHList[nsp].size && ret > 0; ++i)
    if (resHList[nsp].resources[i].status != RESH_UNUSED)
      ret = func(namespaceIdxEncode2(nsp, i),
                 resHList[nsp].resources[i].res.v.val,
                 resHList[nsp].resources[i].res.v.ops, data);
  LIST_UNLOCK();
  return ret;
}


enum cdiApplyRet
cdiResHFilterApply(const resOps *p,
                   enum cdiApplyRet (*func)(int id, void *res, void *data),
                   void *data)
{
  xassert(p && func);

  LIST_INIT(1);

  LIST_LOCK();

  int nsp = namespaceGetActive ();
  enum cdiApplyRet ret = CDI_APPLY_GO_ON;
  listElem_t *r = resHList[nsp].resources;
  for (int i = 0; i < resHList[nsp].size && ret > 0; ++i)
    if (r[i].status != RESH_UNUSED && r[i].res.v.ops == p)
      ret = func(namespaceIdxEncode2(nsp, i), r[i].res.v.val,
                 data);
  LIST_UNLOCK();
  return ret;
}




/**************************************************************/

int reshCountType ( const resOps * ops )
{
  int i, nsp, countType = 0;

  xassert ( ops );

  LIST_INIT(1);

  LIST_LOCK();

  nsp = namespaceGetActive ();

  listElem_t *r = resHList[nsp].resources;
  for ( i = 0; i < resHList[nsp].size; i++ )
    countType += (r[i].status != RESH_UNUSED && r[i].res.v.ops == ops);

  LIST_UNLOCK();

  return countType;
}

/**************************************************************/

int
reshResourceGetPackSize(int resH, const resOps *ops, void *context)
{
  listElem_t *curr = reshGetElem(__func__, resH, ops);
  return curr->res.v.ops->valGetPackSize(curr->res.v.val, context);
}

void
reshPackResource(int resH, const resOps *ops,
                 void *buf, int buf_size, int *position, void *context)
{
  listElem_t *curr = reshGetElem(__func__, resH, ops);
  curr->res.v.ops->valPack(curr->res.v.val, buf, buf_size, position, context);
}


static int getPackBufferSize(void *context)
{
  int nsp, i;
  int intpacksize, packBufferSize = 0;

  nsp = namespaceGetActive ();

  /* pack start marker, namespace and sererator marker */
  packBufferSize += 3 * (intpacksize = serializeGetSize(1, DATATYPE_INT, context));

  /* pack resources, type marker and seperator marker */
  listElem_t *r = resHList[nsp].resources;
  for ( i = 0; i < resHList[nsp].size; i++)
    if (r[i].status == RESH_ASSIGNED)
      {
        xassert ( r[i].res.v.ops );

        /* message plus frame of 2 ints */
        packBufferSize +=
          r[i].res.v.ops->valGetPackSize(r[i].res.v.val, context)
          + 2 * intpacksize;
      }
  /* end marker */
  packBufferSize += intpacksize;

  return packBufferSize;
}

/**************************************************************/

void reshPackBufferDestroy ( char ** buffer )
{
  if ( buffer ) free ( *buffer );
}

/**************************************************************/

void reshPackBufferCreate(char **packBuffer, int *packBufferSize, void *context)
{
  int i, packBufferPos = 0;
  int start = START, end = END, sep = SEPARATOR, type;

  xassert ( packBuffer );

  LIST_LOCK();

  int nsp = namespaceGetActive ();

  int pBSize = *packBufferSize = getPackBufferSize(context);
  char *pB = *packBuffer = (char*) xcalloc(1, *packBufferSize);

  {
    int header[3] = { start, nsp, sep };
    serializePack(header, 3,  DATATYPE_INT, pB, pBSize, &packBufferPos, context);
  }

  listElem_t *r = resHList[nsp].resources;
  for ( i = 0; i < resHList[nsp].size; i++ )
    if ( r[i].status == RESH_ASSIGNED)
      {
        listElem_t * curr = r + i;
        xassert ( curr->res.v.ops );

        type = curr->res.v.ops->valTxCode ();

        if ( ! type ) continue;

        serializePack( &type, 1, DATATYPE_INT, * packBuffer,
                       * packBufferSize, &packBufferPos, context);

        curr->res.v.ops->valPack(curr->res.v.val,
                                 pB, pBSize, &packBufferPos, context);

        serializePack(&sep, 1, DATATYPE_INT, pB, pBSize, &packBufferPos, context);

        curr->status = RESH_CLOSED;
      }

  LIST_UNLOCK();

  serializePack(&end, 1,  DATATYPE_INT, pB, pBSize, &packBufferPos, context);
}

/**************************************************************/

/* for thread safety this feature would have to be integrated in reshPut */

void reshSetStatus ( cdiResH resH, const resOps * ops, int status )
{
  int nsp;
  namespaceTuple_t nspT;
  listElem_t * listElem;

  xassert(ops && status != RESH_UNUSED);

  LIST_INIT(1);

  LIST_LOCK();

  nsp = namespaceGetActive ();

  nspT = namespaceResHDecode ( resH );

  xassert ( nspT.nsp == nsp &&
            nspT.idx >= 0 &&
            nspT.idx < resHList[nsp].size );

  listElem = resHList[nsp].resources + nspT.idx;

  xassert ( listElem &&
            listElem->res.v.ops == ops );

  listElem->status = status;

  LIST_UNLOCK();
}

/**************************************************************/

int reshGetStatus ( cdiResH resH, const resOps * ops )
{
  int nsp;
  namespaceTuple_t nspT;

  xassert ( ops );

  LIST_INIT(1);

  LIST_LOCK();

  nsp = namespaceGetActive ();

  nspT = namespaceResHDecode ( resH );

  xassert ( nspT.nsp == nsp &&
            nspT.idx >= 0 &&
            nspT.idx < resHList[nsp].size );

  listElem_t *listElem = resHList[nsp].resources + nspT.idx;

  const resOps *elemOps = listElem->res.v.ops;

  LIST_UNLOCK();

  xassert(listElem && elemOps == ops);

  return listElem->status;
}

/**************************************************************/

void reshLock ()
{
  LIST_LOCK();
}

/**************************************************************/

void reshUnlock ()
{
  LIST_UNLOCK();
}

/**************************************************************/

int reshListCompare ( int nsp0, int nsp1 )
{
  LIST_INIT(1);
  LIST_LOCK();

  xassert(resHListSize > nsp0 && resHListSize > nsp1 &&
          nsp0 >= 0 && nsp1 >= 0);

  int valCompare = 0;
  int i, listSizeMin = (resHList[nsp0].size <= resHList[nsp1].size)
    ? resHList[nsp0].size : resHList[nsp1].size;
  listElem_t *resources0 = resHList[nsp0].resources,
    *resources1 = resHList[nsp1].resources;
  for (i = 0; i < listSizeMin; i++)
    {
      int occupied0 = resources0[i].status != RESH_UNUSED,
        occupied1 = resources1[i].status != RESH_UNUSED;
      /* occupation mismatch ? */
      int diff = occupied0 ^ occupied1;
      valCompare |= (diff << cdiResHListOccupationMismatch);
      if (!diff && occupied0)
        {
          /* both occupied, do resource types match? */
          diff = (resources0[i].res.v.ops != resources1[i].res.v.ops
                  || resources0[i].res.v.ops == NULL);
          valCompare |= (diff << cdiResHListResourceTypeMismatch);
          if (!diff)
            {
              /* types match, does content match also? */
              diff
                = resources0[i].res.v.ops->valCompare(resources0[i].res.v.val,
                                                      resources1[i].res.v.val);
              valCompare |= (diff << cdiResHListResourceContentMismatch);
            }
        }
    }
  /* find resources in nsp 0 beyond end of nsp 1 */
  for (int j = listSizeMin; j < resHList[nsp0].size; ++j)
    valCompare |= ((resources0[j].status != RESH_UNUSED)
                   << cdiResHListOccupationMismatch);
  /* find resources in nsp 1 beyond end of nsp 0 */
  for (; i < resHList[nsp1].size; ++i)
    valCompare |= ((resources1[i].status != RESH_UNUSED)
                   << cdiResHListOccupationMismatch);

  LIST_UNLOCK();

  return valCompare;
}

/**************************************************************/

void reshListPrint(FILE *fp)
{
  int i, j, temp;
  listElem_t * curr;

  LIST_INIT(1);


  temp = namespaceGetActive ();

  fprintf ( fp, "\n\n##########################################\n#\n#  print " \
            "global resource list \n#\n" );

  for ( i = 0; i < namespaceGetNumber (); i++ )
    {
      namespaceSetActive ( i );

      fprintf ( fp, "\n" );
      fprintf ( fp, "##################################\n" );
      fprintf ( fp, "#\n" );
      fprintf ( fp, "# namespace=%d\n", i );
      fprintf ( fp, "#\n" );
      fprintf ( fp, "##################################\n\n" );

      fprintf ( fp, "resHList[%d].size=%d\n", i, resHList[i].size );

      for ( j = 0; j < resHList[i].size; j++ )
        {
          curr = resHList[i].resources + j;
          if (curr->status != RESH_UNUSED)
            {
              curr->res.v.ops->valPrint(curr->res.v.val, fp);
              fprintf ( fp, "\n" );
            }
        }
    }
  fprintf ( fp, "#\n#  end global resource list" \
            "\n#\n##########################################\n\n" );

  namespaceSetActive ( temp );
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
