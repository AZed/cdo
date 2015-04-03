#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include <errno.h>


#if ! defined (HAVE_CONFIG_H)
#if ! defined (HAVE_MALLOC_H)
#  if defined (SX)
#    define  HAVE_MALLOC_H
#  endif
#endif
#endif

#if  defined  (HAVE_MALLOC_H)
#    include <malloc.h>
#endif


#define  MALLOC_FUNC   0
#define  CALLOC_FUNC   1
#define  REALLOC_FUNC  2
#define  FREE_FUNC     3

#undef   UNDEFID
#define  UNDEFID  -1

#define  MAXNAME  32   /* Min = 8, for  "unknown" ! */

int dmemory_ExitOnError = 0;

typedef struct
{
  void     *ptr;
  int       item;
  size_t    size;
  size_t    nobj;
  int       mtype;
  int       line;
  char      file[MAXNAME];
  char      caller[MAXNAME];
}
MemTable;

static MemTable *memTable;
static int     memTableSize  = 0;
static long    memAccess     = 0;

static size_t  MemObjs       = 0;
static size_t  MaxMemObjs    = 0;
static size_t  MemUsed       = 0;
static size_t  MaxMemUsed    = 0;

static int     MEM_Debug     = 0;   /* If set to 1, debugging */

void memDebug(int debug)
{
  MEM_Debug = debug;
}

static
void memInternalProblem(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  printf("\n");
   fprintf(stderr, "Internal problem (%s) : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  exit(EXIT_FAILURE);
}

static
void memError(const char *caller, const char *file, int line, size_t size)
{
  printf("\n");
  fprintf(stderr, "Error (%s) : Allocation of %lu bytes failed. [ line %d file %s ]\n",
	  caller, (unsigned long) size, line, file);

  if ( errno )
    perror("System error message ");
	
  exit(EXIT_FAILURE);
}

static
void memListPrintEntry(int mtype, int item, size_t size, void *ptr,
		       const char *caller, const char *file, int line)
{
  switch (mtype)
    {
    case MALLOC_FUNC:
      fprintf(stderr, "[%-7s ", "Malloc");
      break;
    case CALLOC_FUNC:
      fprintf(stderr, "[%-7s ", "Calloc");
      break;
    case REALLOC_FUNC:
      fprintf(stderr, "[%-7s ", "Realloc");
      break;
    case FREE_FUNC:
      fprintf(stderr, "[%-7s ", "Free");
      break;
    }

   fprintf(stderr, "memory item %3d ", item);
   fprintf(stderr, "(%6lu byte) ", (unsigned long) size);
   fprintf(stderr, "at %p", ptr);
   if ( file != NULL )
     {
       fprintf(stderr, " line %4d", line);
       fprintf(stderr, " file %s", file);
     }    
   if ( caller != NULL )
     fprintf(stderr, " (%s)", caller);     
   fprintf(stderr, "]\n");
}

static
void memListPrintTable(void)
{
  int memID, item, item1, item2 = 0;

  if ( MemObjs ) fprintf(stderr, "\nMemory table:\n");

  /* find maximum item */
  for ( memID = 0; memID < memTableSize; memID++ )
    if ( memTable[memID].item != UNDEFID )
      if ( memTable[memID].item > item2 ) item2 = memTable[memID].item;

  /* find minimum item */
  item1 = item2;
  for ( memID = 0; memID < memTableSize; memID++ )
    if ( memTable[memID].item != UNDEFID )
      if ( memTable[memID].item < item1 ) item1 = memTable[memID].item;

  for ( item = item1; item <= item2; item++ )
    for ( memID = 0; memID < memTableSize; memID++ )
      {
	if ( memTable[memID].item == item )
	  memListPrintEntry(memTable[memID].mtype, memTable[memID].item,
			    memTable[memID].size*memTable[memID].nobj,
			    memTable[memID].ptr, memTable[memID].caller,
			    memTable[memID].file, memTable[memID].line);
      }

  if ( MemObjs )
    {
      fprintf(stderr, "  Memory access             : %6u\n", (unsigned) memAccess);
      fprintf(stderr, "  Maximum objects           : %6u\n", (unsigned) memTableSize);
      fprintf(stderr, "  Objects used              : %6u\n", (unsigned) MaxMemObjs);
      fprintf(stderr, "  Objects in use            : %6u\n", (unsigned) MemObjs);
      fprintf(stderr, "  Memory allocated          : ");
      if (MemUsed > 1024*1024*1024)
	fprintf(stderr, " %5d GB\n",   (int) (MemUsed/(1024*1024*1024)));
      else if (MemUsed > 1024*1024)
	fprintf(stderr, " %5d MB\n",   (int) (MemUsed/(1024*1024)));
      else if (MemUsed > 1024)
	fprintf(stderr, " %5d KB\n",   (int) (MemUsed/(1024)));
      else
	fprintf(stderr, " %5d Byte\n", (int)  MemUsed);
    }

  if ( MaxMemUsed )
    {
      fprintf(stderr, "  Maximum memory allocated  : ");
      if (MaxMemUsed > 1024*1024*1024)
	fprintf(stderr, " %5d GB\n",   (int) (MaxMemUsed/(1024*1024*1024)));
      else if (MaxMemUsed > 1024*1024)
	fprintf(stderr, " %5d MB\n",   (int) (MaxMemUsed/(1024*1024)));
      else if (MaxMemUsed > 1024)
	fprintf(stderr, " %5d KB\n",   (int) (MaxMemUsed/(1024)));
      else
	fprintf(stderr, " %5d Byte\n", (int)  MaxMemUsed);
    }
}

static
void memGetDebugLevel(void)
{
  char *debugLevel;

  debugLevel = getenv("MEMORY_DEBUG");

  if ( debugLevel )
    {
      if ( isdigit((int) debugLevel[0]) )
	MEM_Debug = atoi(debugLevel);

      if ( MEM_Debug )
	atexit(memListPrintTable);
    }
}

static
void memInit(void)
{
  static int initDebugLevel = 0;

  if ( ! initDebugLevel )
    {
      memGetDebugLevel();
      initDebugLevel = 1;
    }  
}

static
int memListDeleteEntry(void *ptr, size_t *size)
{
  int memID = 0;
  int item = UNDEFID;

  for ( memID = 0; memID < memTableSize; memID++ )
    {
      if ( memTable[memID].item == UNDEFID ) continue;
      if ( memTable[memID].ptr == ptr ) break;
    }

  if ( memID != memTableSize )
    {
      MemObjs--;
      MemUsed -= memTable[memID].size * memTable[memID].nobj;
      *size = memTable[memID].size * memTable[memID].nobj;
       item = memTable[memID].item;
       memTable[memID].item   = UNDEFID;
    }

  return (item);
}

static
void memTableInitEntry(int memID)
{
  if ( memID < 0 || memID >= memTableSize )
    memInternalProblem(__func__, "memID %d undefined!", memID);

  memTable[memID].ptr    = NULL;
  memTable[memID].item   = UNDEFID;
  memTable[memID].size   = 0;
  memTable[memID].nobj   = 0;
  memTable[memID].mtype  = UNDEFID;
  memTable[memID].line   = UNDEFID;
}

static
int memListNewEntry(int mtype, void *ptr, size_t size, size_t nobj,
		    const char *caller, const char *file, int line)
{
  static int item = 0;
  size_t memSize = 0;
  int memID = 0;
  size_t len;
  int i;

  /*
    Look for a free slot in memTable.
    (Create the table the first time through).
  */
  if ( memTableSize == 0 )
    {
      memTableSize = 8;
      memSize  = memTableSize*sizeof(MemTable);
      memTable = (MemTable *) malloc(memSize);
      if( memTable == NULL ) memError(__func__, __FILE__, __LINE__, memSize);

      for( i = 0; i < memTableSize; i++ )
	memTableInitEntry(i);
    }
  else
    {
      while( memID < memTableSize )
	{
	  if ( memTable[memID].item == UNDEFID ) break;
	  memID++;
	}
    }
  /*
    If the table overflows, double its size.
  */
  if ( memID == memTableSize )
    {
      memTableSize = 2*memTableSize;
      memSize  = memTableSize*sizeof(MemTable);
      memTable = (MemTable *) realloc(memTable, memSize);
      if( memTable == NULL ) memError(__func__, __FILE__, __LINE__, memSize);

      for( i = memID; i < memTableSize; i++ )
	memTableInitEntry(i);
    }

  memTable[memID].item  = item;
  memTable[memID].ptr   = ptr;
  memTable[memID].size  = size;
  memTable[memID].nobj  = nobj;
  memTable[memID].mtype = mtype;
  memTable[memID].line  = line;

  if ( file )
    {
      len = strlen(file);
      if ( len > MAXNAME-1 ) len = MAXNAME-1;
    
      (void) memcpy(memTable[memID].file, file, len);
      memTable[memID].file[len] = '\0';
    }
  else
    {
      (void) strcpy(memTable[memID].file, "unknown");
    }

  if ( caller )
    {
      len = strlen(caller);
      if ( len > MAXNAME-1 ) len = MAXNAME-1;

      (void) memcpy(memTable[memID].caller, caller, len);
      memTable[memID].caller[len] = '\0';
    }
  else
    {
      (void) strcpy(memTable[memID].caller, "unknown");
    }

  MaxMemObjs++;
  MemObjs++;
  MemUsed += size*nobj;
  if ( MemUsed > MaxMemUsed ) MaxMemUsed = MemUsed;

  return (item++);
}

static
int memListChangeEntry(void *ptrold, void *ptr, size_t size,
		       const char *caller, const char *file, int line)
{
  int item = UNDEFID;
  int memID = 0;
  size_t len;
  size_t sizeold;

  while( memID < memTableSize )
    {
      if ( memTable[memID].item != UNDEFID )
	if ( memTable[memID].ptr == ptrold ) break;
      memID++;
    }

  if ( memID == memTableSize )
    {
      if ( ptrold != NULL )
	memInternalProblem(__func__, "Item at %p not found.", ptrold);
    }
  else
    {
      item = memTable[memID].item;

      sizeold = memTable[memID].size*memTable[memID].nobj;
      
      memTable[memID].ptr   = ptr;
      memTable[memID].size  = size;
      memTable[memID].nobj  = 1;
      memTable[memID].mtype = REALLOC_FUNC;
      memTable[memID].line  = line;

      if ( file )
	{
	  len = strlen(file);
	  if ( len > MAXNAME-1 ) len = MAXNAME-1;

	  (void) memcpy(memTable[memID].file, file, len);
	  memTable[memID].file[len] = '\0';
	}
      else
	{
	  (void) strcpy(memTable[memID].file, "unknown");
	}

      if ( caller )
	{
	  len = strlen(caller);
	  if ( len > MAXNAME-1 ) len = MAXNAME-1;

	  (void) memcpy(memTable[memID].caller, caller, len);
	  memTable[memID].caller[len] = '\0';
	}
      else
	{
	  (void) strcpy(memTable[memID].caller, "unknown");
	}

      MemUsed -= sizeold;
      MemUsed += size;
      if ( MemUsed > MaxMemUsed ) MaxMemUsed = MemUsed;
    }

  return (item);
}


void *Calloc(const char *caller, const char *file, int line, size_t nobjs, size_t size)
{
  void *ptr = NULL;
  int item = UNDEFID;

  memInit();

  if ( nobjs*size > 0 )
    {
      ptr = calloc(nobjs, size);

      if ( MEM_Debug )
	{
	  memAccess++;

	  if ( ptr )
	    item = memListNewEntry(CALLOC_FUNC, ptr, size, nobjs, caller, file, line);

	  memListPrintEntry(CALLOC_FUNC, item, size*nobjs, ptr, caller, file, line);
	}

      if ( ptr == NULL && dmemory_ExitOnError )
	memError(caller, file, line, size*nobjs);
    }
  else
    fprintf(stderr, "Warning (%s) : Allocation of 0 bytes! [ line %d file %s ]\n", caller, line, file);

  return(ptr);
}


void *Malloc(const char *caller, const char *file, int line, size_t size)
{
  void *ptr = NULL;
  int item = UNDEFID;

  memInit();

  if ( size > 0 )
    {
      ptr = malloc(size);

      if ( MEM_Debug )
	{
	  memAccess++;

	  if ( ptr )
	    item = memListNewEntry(MALLOC_FUNC, ptr, size, 1, caller, file, line);

	  memListPrintEntry(MALLOC_FUNC, item, size, ptr, caller, file, line);
	}

      if ( ptr == NULL && dmemory_ExitOnError )
	memError(caller, file, line, size);
    }
  else
    fprintf(stderr, "Warning (%s) : Allocation of 0 bytes! [ line %d file %s ]\n", caller, line, file);

  return (ptr);
}


void *Realloc(const char *caller, const char *file, int line, void *ptrold, size_t size)
{
  void *ptr = NULL;
  int item = UNDEFID;

  memInit();

  if ( size > 0 )
    {
      ptr = realloc(ptrold, size);

      if ( MEM_Debug )
	{
	  memAccess++;

	  if ( ptr )
	    {
	      item = memListChangeEntry(ptrold, ptr, size, caller, file, line);

	      if ( item == UNDEFID )
		item = memListNewEntry(REALLOC_FUNC, ptr, size, 1, caller, file, line);
	    }

	  memListPrintEntry(REALLOC_FUNC, item, size, ptr, caller, file, line);
	}

      if ( ptr == NULL && dmemory_ExitOnError )
	memError(caller, file, line, size);
    }
  else
    fprintf(stderr, "Warning (%s) : Allocation of 0 bytes! [ line %d file %s ]\n", caller, line, file);

  return (ptr);
}


void Free(const char *caller, const char *file, int line, void *ptr)
{
  int item;
  size_t size;

  memInit();

  if ( MEM_Debug )
    {
      if ( (item = memListDeleteEntry(ptr, &size)) >= 0 )
	{
	  memListPrintEntry(FREE_FUNC, item, size, ptr, caller, file, line);
	}
      else
	{
	  if ( ptr )
	    fprintf(stderr, "%s info: memory entry at %p not found. [line %4d file %s (%s)]\n",
		    __func__, ptr, line, file, caller);
	}
    }

  free(ptr);
}


size_t memTotal(void)
{
  size_t memtotal = 0;
#if  defined  (HAVE_MALLINFO)
  struct mallinfo meminfo = mallinfo();
  if ( MEM_Debug )
    {
      fprintf(stderr, "arena      %8ld (non-mmapped space allocated from system)\n", (unsigned long) meminfo.arena);
      fprintf(stderr, "ordblks    %8ld (number of free chunks)\n", (unsigned long) meminfo.ordblks);
      fprintf(stderr, "smblks     %8ld (number of fastbin blocks)\n", (unsigned long) meminfo.smblks);
      fprintf(stderr, "hblks      %8ld (number of mmapped regions)\n", (unsigned long) meminfo.hblks);
      fprintf(stderr, "hblkhd     %8ld (space in mmapped regions)\n", (unsigned long) meminfo.hblkhd);
      fprintf(stderr, "usmblks    %8ld (maximum total allocated space)\n", (unsigned long) meminfo.usmblks);
      fprintf(stderr, "fsmblks    %8ld (maximum total allocated space)\n", (unsigned long) meminfo.fsmblks);
      fprintf(stderr, "uordblks   %8ld (total allocated space)\n", (unsigned long) meminfo.uordblks);
      fprintf(stderr, "fordblks   %8ld (total free space)\n", (unsigned long) meminfo.fordblks);
      fprintf(stderr, "Memory in use:   %8ld bytes\n", (unsigned long) meminfo.usmblks + meminfo.uordblks);
      fprintf(stderr, "Total heap size: %8ld bytes\n", (unsigned long) meminfo.arena);

      /* malloc_stats(); */
    }
  memtotal = meminfo.arena;
#endif

  return (memtotal);
}


void memExitOnError(void)
{
  dmemory_ExitOnError = 1;
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
