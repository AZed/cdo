#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 /* PTHREAD_MUTEX_RECURSIVE */
#endif

#include <stdlib.h>
#include <stdio.h>

#include "resource_handle.h"
#include "pio_util.h"
#include "namespace.h"
#include "serialize.h"
#include "cdi.h"
#include "error.h"
#include "file.h"
#include "resource_unpack.h"

enum { MIN_LIST_SIZE = 128 };

static void listInitialize(void);

// ATTENTION: not thread safe yet, namespaces are set in model!

#if  defined  (HAVE_LIBPTHREAD)
#  include <pthread.h>

static pthread_once_t  listInitThread = PTHREAD_ONCE_INIT;
static pthread_mutex_t listMutex;

#  define LIST_LOCK()         pthread_mutex_lock(&listMutex)
#  define LIST_UNLOCK()       pthread_mutex_unlock(&listMutex)
#  define LIST_INIT()         pthread_once(&listInitThread, listInitialize)

#else

static int listInit = 0;

#  define LIST_LOCK()
#  define LIST_UNLOCK()
#  define LIST_INIT()        do {                               \
  if ( !listInit )                                              \
    {                                                           \
      listInitialize();                                         \
      listInit = 1;                                             \
    }                                                           \
  } while(0)

#endif


typedef struct listElem {
  cdiResH       resH;//idx
  int next;
  resOps      * ops;
  void        * val;//ptr
  int           status;
} listElem_t;

static struct
{
  int size, freeHead;
  listElem_t *resources;
} *resHList;

static int resHListSize = 0;

/**************************************************************/

static void
listInitResources(int nsp)
{
  xassert(nsp < resHListSize && nsp >= 0);
  int size = resHList[nsp].size = MIN_LIST_SIZE;
  xassert(resHList[nsp].resources == NULL);
  resHList[nsp].resources = xcalloc(MIN_LIST_SIZE, sizeof (listElem_t));
  listElem_t *p = resHList[nsp].resources;

  for (int i = 0; i < size; i++ )
    {
      p[i].resH   = namespaceIdxEncode2(nsp, i);
      p[i].next   = i + 1;
      p[i].ops    = NULL;
      p[i].val    = NULL;
      p[i].status = RESH_UNDEFID;
    }

  p[resHList[nsp].size-1].next = -1;
  resHList[nsp].freeHead = 0;
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
  LIST_LOCK();
  if (resHListSize <= namespaceID)
    {
      resHList = xrealloc(resHList, (namespaceID + 1) * sizeof (resHList[0]));
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
  pioNamespaceSetActive(namespaceID);
  if (resHList[namespaceID].resources)
    {
      for ( int j = 0; j < resHList[namespaceID].size; j++ )
        {
          listElem_t *listElem = resHList[namespaceID].resources + j;
          if (listElem->val)
            listElem->ops->valDestroy(listElem->val);
        }
      free(resHList[namespaceID].resources);
      reshListClearEntry(namespaceID);
    }
  if (resHList[callerNamespaceID].resources)
    pioNamespaceSetActive(callerNamespaceID);
  LIST_UNLOCK();
}


static void listDestroy ( void )
{
  LIST_LOCK();
  for (int i = 0; i < resHListSize; ++i)
    if (resHList[i].resources)
      namespaceDelete(i);
  free(resHList);
  resHList = NULL;
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
  // create default namespace
  reshListCreate(0);
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

  resHList[nsp].resources =
    xrealloc(resHList[nsp].resources,
             newListSize * sizeof (resHList[0].resources[0]));

  for (int i = oldSize; i < newListSize; ++i)
    {
      resHList[nsp].resources[i].resH   = namespaceIdxEncode2 ( nsp, i );
      resHList[nsp].resources[i].next   = i + 1;
      resHList[nsp].resources[i].ops    = NULL;
      resHList[nsp].resources[i].val    = NULL;
      resHList[nsp].resources[i].status = RESH_UNDEFID;
    }

  resHList[nsp].resources[newListSize-1].next = resHList[nsp].freeHead;
  resHList[nsp].freeHead = oldSize;
  resHList[nsp].size = newListSize;
}

/**************************************************************/

int reshPut ( void *p, resOps *ops )
{
  cdiResH resH = -1, nsp;
  listElem_t * newListElem;

  xassert ( p && ops );

  LIST_INIT();

  LIST_LOCK();

  nsp = namespaceGetActive ();

  if ( resHList[nsp].freeHead == -1) listSizeExtend();
  newListElem               = resHList[nsp].resources + resHList[nsp].freeHead;
  resHList[nsp].freeHead    = newListElem->next;
  newListElem->next         = -1;
  resH                      = newListElem->resH;
  newListElem->val          = p;
  newListElem->ops          = ops;
  newListElem->status       = ASSIGNED;

  LIST_UNLOCK();

  return resH;
}

/**************************************************************/

void reshRemove ( cdiResH resH, resOps * ops )
{
  int nsp;
  namespaceTuple_t nspT;

  LIST_INIT();

  LIST_LOCK();

  nsp = namespaceGetActive ();

  nspT = namespaceResHDecode ( resH );

  xassert ( nspT.nsp == nsp &&
            nspT.idx >= 0 &&
            nspT.idx < resHList[nsp].size &&
            resHList[nsp].resources[nspT.idx].ops &&
            resHList[nsp].resources[nspT.idx].ops == ops );

  resHList[nsp].resources[nspT.idx].next   = resHList[nsp].freeHead;
  resHList[nsp].resources[nspT.idx].ops    = NULL;
  resHList[nsp].resources[nspT.idx].val    = NULL;
  resHList[nsp].resources[nspT.idx].status = RESH_UNDEFID;
  resHList[nsp].freeHead                   = nspT.idx;

  LIST_UNLOCK();
}

/**************************************************************/

static listElem_t *
reshGetElem(const char *caller, cdiResH resH, resOps *ops)
{
  listElem_t *listElem;
  int nsp;
  namespaceTuple_t nspT;
  xassert ( ops );

  LIST_INIT();

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

  if ( !(listElem && listElem->ops == ops) )
    xabortC(caller, "Invalid resource handle %d, list element not found!",
            (int)resH);
  return listElem;
}


void *reshGetValue(const char * caller, cdiResH resH, resOps * ops)
{
  return reshGetElem(caller, resH, ops)->val;
}

/**************************************************************/

void reshGetResHListOfType ( int c, int * resHs, resOps * ops )
{
  int i, j = 0, nsp;

  xassert ( resHs && ops );

  LIST_INIT();

  LIST_LOCK();

  nsp = namespaceGetActive ();

  for ( i = 0; i < resHList[nsp].size; i++ )
    if ( resHList[nsp].resources[i].val && resHList[nsp].resources[i].ops )
      if ( resHList[nsp].resources[i].ops == ops )
        {
          resHs[j++] = namespaceIdxEncode2 ( nsp, i );
          if ( j == c ) break;
        }

  LIST_UNLOCK();
}

enum cdiApplyRet
cdiResHApply(enum cdiApplyRet (*func)(int id, void *res, const resOps *p,
                                      void *data), void *data)
{
  xassert(func);

  LIST_INIT();

  LIST_LOCK();

  int nsp = namespaceGetActive ();
  enum cdiApplyRet ret = CDI_APPLY_GO_ON;
  for (int i = 0; i < resHList[nsp].size && ret > 0; ++i)
    if (resHList[nsp].resources[i].val && resHList[nsp].resources[i].ops)
      ret = func(namespaceIdxEncode2(nsp, i), resHList[nsp].resources[i].val,
                 resHList[nsp].resources[i].ops, data);
  LIST_UNLOCK();
  return ret;
}


enum cdiApplyRet
cdiResHFilterApply(const resOps *p,
                   enum cdiApplyRet (*func)(int id, void *res, void *data),
                   void *data)
{
  xassert(p && func);

  LIST_INIT();

  LIST_LOCK();

  int nsp = namespaceGetActive ();
  enum cdiApplyRet ret = CDI_APPLY_GO_ON;
  for (int i = 0; i < resHList[nsp].size && ret > 0; ++i)
    if (resHList[nsp].resources[i].val && resHList[nsp].resources[i].ops
        && resHList[nsp].resources[i].ops == p)
      ret = func(namespaceIdxEncode2(nsp, i), resHList[nsp].resources[i].val,
                 data);
  LIST_UNLOCK();
  return ret;
}




/**************************************************************/

int reshCountType ( resOps * ops )
{
  int i, nsp, countType = 0;

  xassert ( ops );

  LIST_INIT();

  LIST_LOCK();

  nsp = namespaceGetActive ();

  for ( i = 0; i < resHList[nsp].size; i++ )
    if ( resHList[nsp].resources[i].val )
      if ( resHList[nsp].resources[i].ops == ops )
        countType++;

  LIST_UNLOCK();

  return countType;
}

/**************************************************************/

int
reshResourceGetPackSize(int resH, resOps *ops, void *context)
{
  listElem_t *curr = reshGetElem(__func__, resH, ops);
  return curr->ops->valGetPackSize(curr->val, context);
}

void
reshPackResource(int resH, resOps *ops,
                 void *buf, int buf_size, int *position, void *context)
{
  listElem_t *curr = reshGetElem(__func__, resH, ops);
  curr->ops->valPack(curr->val, buf, buf_size, position, context);
}


static int getPackBufferSize(void *context)
{
  int nsp, i;
  int intpacksize, packBufferSize = 0;
  listElem_t * curr;

  nsp = namespaceGetActive ();

  /* pack start marker, namespace and sererator marker */
  packBufferSize += (intpacksize = serializeGetSize(3, DATATYPE_INT, context));

  /* pack resources, type marker and seperator marker */
  for ( i = 0; i < resHList[nsp].size; i++ )
    if ( resHList[nsp].resources[i].val )
      if ( resHList[nsp].resources[i].status == ASSIGNED )
        {
          curr = resHList[nsp].resources + i;
          xassert ( curr->ops );

          /* message plus frame of 2 ints */
          packBufferSize += curr->ops->valGetPackSize(curr->val, context)
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
  listElem_t * curr;

  xassert ( packBuffer );

  LIST_LOCK();

  int nsp = namespaceGetActive ();

  int pBSize = *packBufferSize = getPackBufferSize(context);
  char *pB = *packBuffer = xcalloc(1, *packBufferSize);

  {
    int header[3] = { start, nsp, sep };
    serializePack(header, 3,  DATATYPE_INT, pB, pBSize, &packBufferPos, context);
  }

  for ( i = 0; i < resHList[nsp].size; i++ )
    if ( resHList[nsp].resources[i].val )
      if ( resHList[nsp].resources[i].status == ASSIGNED )
        {
          curr = resHList[nsp].resources + i;
          xassert ( curr->ops );

          type = curr->ops->valTxCode ();

          if ( ! type ) continue;

          serializePack( &type, 1, DATATYPE_INT, * packBuffer,
			    * packBufferSize, &packBufferPos, context);

          curr->ops->valPack(curr->val, pB, pBSize, &packBufferPos, context);

          serializePack(&sep, 1, DATATYPE_INT, pB, pBSize, &packBufferPos, context);

          curr->status = CLOSED;
        }

  LIST_UNLOCK();

  serializePack(&end, 1,  DATATYPE_INT, pB, pBSize, &packBufferPos, context);
}

/**************************************************************/

/* for thread safety this feature would have to be integrated in reshPut */

void reshSetStatus ( cdiResH resH, resOps * ops, int status )
{
  int nsp;
  namespaceTuple_t nspT;
  listElem_t * listElem;

  xassert ( ops );

  LIST_INIT();

  LIST_LOCK();

  nsp = namespaceGetActive ();

  nspT = namespaceResHDecode ( resH );

  xassert ( nspT.nsp == nsp &&
            nspT.idx >= 0 &&
            nspT.idx < resHList[nsp].size );

  listElem = resHList[nsp].resources + nspT.idx;

  xassert ( listElem &&
            listElem->ops == ops );

  listElem->status = status;

  LIST_UNLOCK();
}

/**************************************************************/

int reshGetStatus ( cdiResH resH, resOps * ops )
{
  int nsp;
  namespaceTuple_t nspT;
  listElem_t * listElem;

  xassert ( ops );

  LIST_INIT();

  LIST_LOCK();

  nsp = namespaceGetActive ();

  nspT = namespaceResHDecode ( resH );

  xassert ( nspT.nsp == nsp &&
            nspT.idx >= 0 &&
            nspT.idx < resHList[nsp].size );

  listElem = resHList[nsp].resources + nspT.idx;

  LIST_UNLOCK();

  xassert ( listElem &&
            listElem->ops == ops );

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
  int valCompare = 0;
  int i;


  LIST_INIT();
  LIST_LOCK();

  xassert(resHListSize > xmaxInt ( nsp0, nsp1 ) &&
          xminInt ( nsp0, nsp1 ) >= 0 );

  for ( i = 0; i < resHList[nsp0].size; i++ )
    {
      listElem_t *listElem0 = resHList[nsp0].resources + i,
        *listElem1 = resHList[nsp1].resources + i;
      if ( listElem0->val )
	{
	  if ( i >= resHList[nsp1].size )
	    {
              valCompare = 1;
	      xdebug("%s %d", "namespace active length mismatch at resource",
                     i);
	      break;
	    }

	  if ( !listElem1->val )
	    {
              valCompare = 1;
	      xdebug("%s %d", "namespace occupation mismatch at resource", i);
              break;
	    }

	  if ( listElem0->ops != listElem1->ops || listElem0->ops == NULL )
	    {
              valCompare = 1;
	      xdebug("%s %d", "resource type mismatch at resource", i);
              break;
	    }

	  valCompare = listElem0->ops->valCompare(listElem0->val,
                                                  listElem1->val);
          if (valCompare)
            break;
	}
      else if ( listElem1->val )
        {
          valCompare = 1;
          xdebug("namespace 1 has value at empty place %d of namespace 0",
                 i);
          break;
        }
    }

  if (!valCompare)
    {
      for ( ; i < resHList[nsp1].size; i++ )
        valCompare = valCompare || resHList[nsp1].resources[i].val != NULL;
      if (valCompare)
        xdebug("%s", "extra elements in second namespace");
    }

  LIST_UNLOCK();

  return valCompare;
}

/**************************************************************/

void reshListPrint(FILE *fp)
{
  int i, j, temp;
  listElem_t * curr;

  LIST_INIT();


  temp = namespaceGetActive ();

  fprintf ( fp, "\n\n##########################################\n#\n#  print " \
            "global resource list \n#\n" );

  for ( i = 0; i < namespaceGetNumber (); i++ )
    {
      pioNamespaceSetActive ( i );

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
          if ( curr->ops && curr->val )
            {
              curr->ops->valPrint (( void * ) curr->val, fp );
              fprintf ( fp, "\n" );
            }
        }
    }
  fprintf ( fp, "#\n#  end global resource list" \
            "\n#\n##########################################\n\n" );

  pioNamespaceSetActive ( temp );
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
