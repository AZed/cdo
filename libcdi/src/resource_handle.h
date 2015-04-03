#ifndef RESOURCE_HANDLE_H
#define RESOURCE_HANDLE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>

/*
 * CDI internal handling of resource handles given to user code
 */

/*
 * for reasons of compatibility with cfortran.h, the handle type is: int
 */
typedef int cdiResH;

/* return 0 on equality, not 0 otherwise */
typedef int    ( * valCompareFunc     )( void *, void * );
typedef void   ( * valDestroyFunc     )( void * );
typedef void   ( * valPrintFunc       )( void *, FILE * );
typedef int    ( * valGetPackSizeFunc )( void *, void *context );
typedef void   ( * valPackFunc        )( void *, void *buf, int size, int *pos, void *context );
typedef int    ( * valTxCodeFunc      )( void );

typedef struct {
  valCompareFunc     valCompare;
  valDestroyFunc     valDestroy;
  valPrintFunc       valPrint;
  valGetPackSizeFunc valGetPackSize;
  valPackFunc        valPack;
  valTxCodeFunc      valTxCode;
}resOps;

enum {
  RESH_UNUSED,             /* resource holds no value */
  RESH_ASSIGNED,           /* resource is user-assigned */
  RESH_PRE_ASSIGNED,       /* resource is pre-assigned by library */
  RESH_CLOSED,             /* resource is marked immutable */
};

void   reshListCreate(int namespaceID);
void   reshListDestruct(int namespaceID);
int    reshPut ( void *, const resOps * );
void reshReplace(cdiResH resH, void *p, const resOps *ops);
void   reshRemove ( cdiResH, const resOps * );

int    reshCountType ( const resOps * );

void * reshGetValue(const char *, cdiResH, const resOps * );
#define reshGetVal(resH, ops)  reshGetValue(__func__, resH, ops)

void   reshGetResHListOfType ( int, int *, const resOps * );

enum cdiApplyRet {
  CDI_APPLY_ERROR = -1,
  CDI_APPLY_STOP,
  CDI_APPLY_GO_ON,
};
enum cdiApplyRet
cdiResHApply(enum cdiApplyRet (*func)(int id, void *res, const resOps *p,
                                      void *data), void *data);
enum cdiApplyRet
cdiResHFilterApply(const resOps *p,
                   enum cdiApplyRet (*func)(int id, void *res,
                                            void *data),
                   void *data);

void   reshPackBufferCreate ( char **, int *, void *context );
void   reshPackBufferDestroy ( char ** );
int    reshResourceGetPackSize(int resh, const resOps *ops, void *context);
void   reshPackResource(int resh, const resOps *ops,
                        void *buf, int buf_size, int *position, void *context);
void   reshSetStatus ( cdiResH, const resOps *, int );
int    reshGetStatus ( cdiResH, const resOps * );

void   reshLock   ( void );
void   reshUnlock ( void );

enum reshListMismatch {
  cdiResHListOccupationMismatch,
  cdiResHListResourceTypeMismatch,
  cdiResHListResourceContentMismatch,
};

int reshListCompare(int nsp0, int nsp1);
void reshListPrint(FILE *fp);

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