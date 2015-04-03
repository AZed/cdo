#ifndef _DMEMORY_H
#define _DMEMORY_H

#include <stdlib.h>


/*
 * if DEBUG_MEMORY is defined setenv MEMORY_DEBUG to debug memory
 */

#define  DEBUG_MEMORY

#define  WITH_CALLER_NAME

extern size_t  memTotal(void);
extern void    memDebug(int debug);
extern void    memExitOnError(void);

#if  defined  DEBUG_MEMORY

extern void   *Realloc(const char *caller, const char *file, int line, void *ptr, size_t size);
extern void   *Calloc (const char *caller, const char *file, int line, size_t nmemb, size_t size);
extern void   *Malloc (const char *caller, const char *file, int line, size_t size);
extern void    Free   (const char *caller, const char *file, int line, void *ptr);

#if  defined  calloc
#  undef  calloc
#endif

#if  defined  WITH_CALLER_NAME
#  define  realloc(x, y)  Realloc(func, __FILE__, __LINE__, x, y)
#  define   calloc(x, y)   Calloc(func, __FILE__, __LINE__, x, y)
#  define   malloc(x)      Malloc(func, __FILE__, __LINE__, x)
#  define     free(x)        Free(func, __FILE__, __LINE__, x)
#else
#  define  realloc(x, y)  Realloc((void *) NULL, __FILE__, __LINE__, x, y)
#  define   calloc(x, y)   Calloc((void *) NULL, __FILE__, __LINE__, x, y)
#  define   malloc(x)      Malloc((void *) NULL, __FILE__, __LINE__, x)
#  define     free(x)        Free((void *) NULL, __FILE__, __LINE__, x)
#endif

#endif /* DEBUG_MEMORY */

#endif /* _DMEMORY_H */
