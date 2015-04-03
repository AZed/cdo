#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "dmemory.h"
#include "cdi.h"
#include "stream_int.h"

#undef  UNDEFID
#define UNDEFID -1

int ECHAM4 = UNDEFID;
int ECHAM5 = UNDEFID;
int COSMO  = UNDEFID;

typedef struct
{
  int      self;
  int      used;  
  int      instID;  
  int      modelgribID;
  char    *name;
}
model_t;


static int  MODEL_Debug = 0;   /* If set to 1, debugging */

static int _model_max = MAX_MODELS;

static void model_initialize(void);

static int _model_init = FALSE;

#if  defined  (HAVE_LIBPTHREAD)
#  include <pthread.h>

static pthread_once_t _model_init_thread = PTHREAD_ONCE_INIT;
static pthread_mutex_t _model_mutex;

#  define MODEL_LOCK           pthread_mutex_lock(&_model_mutex);
#  define MODEL_UNLOCK         pthread_mutex_unlock(&_model_mutex);
#  define MODEL_INIT                               \
   if ( _model_init == FALSE ) pthread_once(&_model_init_thread, model_initialize);

#else

#  define MODEL_LOCK
#  define MODEL_UNLOCK
#  define MODEL_INIT                               \
   if ( _model_init == FALSE ) model_initialize();

#endif


typedef struct _modelPtrToIdx {
  int idx;
  model_t *ptr;
  struct _modelPtrToIdx *next;
} modelPtrToIdx;


static modelPtrToIdx *_modelList  = NULL;
static modelPtrToIdx *_modelAvail = NULL;


static
void model_list_new(void)
{
  assert(_modelList == NULL);

  _modelList = (modelPtrToIdx *) malloc(_model_max*sizeof(modelPtrToIdx));
}

static
void model_list_delete(void)
{
  if ( _modelList ) free(_modelList);
}

static
void model_init_pointer(void)
{
  int  i;
  
  for ( i = 0; i < _model_max; i++ )
    {
      _modelList[i].next = _modelList + i + 1;
      _modelList[i].idx  = i;
      _modelList[i].ptr  = 0;
    }

  _modelList[_model_max-1].next = 0;

  _modelAvail = _modelList;
}

static
model_t *model_to_pointer(int idx)
{
  model_t *modelptr = NULL;

  MODEL_INIT

  if ( idx >= 0 && idx < _model_max )
    {
      MODEL_LOCK

      modelptr = _modelList[idx].ptr;

      MODEL_UNLOCK
    }
  else
    Error("model index %d undefined!", idx);

  return (modelptr);
}

/* Create an index from a pointer */
static
int model_from_pointer(model_t *ptr)
{
  int      idx = -1;
  modelPtrToIdx *newptr;

  if ( ptr )
    {
      MODEL_LOCK

      if ( _modelAvail )
	{
	  newptr       = _modelAvail;
	  _modelAvail  = _modelAvail->next;
	  newptr->next = 0;
	  idx	       = newptr->idx;
	  newptr->ptr  = ptr;
      
	  if ( MODEL_Debug )
	    Message("Pointer %p has idx %d from model list", ptr, idx);
	}
      else
	Warning("Too many open models (limit is %d)!", _model_max);

      MODEL_UNLOCK
    }
  else
    Error("Internal problem (pointer %p undefined)", ptr);

  return (idx);
}

static
void model_init_entry(model_t *modelptr)
{
  modelptr->self        = model_from_pointer(modelptr);

  modelptr->used        = 1;

  modelptr->instID      = UNDEFID;
  modelptr->modelgribID = UNDEFID;
  modelptr->name        = NULL;
}

static
model_t *model_new_entry(void)
{
  model_t *modelptr;

  modelptr = (model_t *) malloc(sizeof(model_t));

  if ( modelptr ) model_init_entry(modelptr);

  return (modelptr);
}

static
void model_delete_entry(model_t *modelptr)
{
  int idx;

  idx = modelptr->self;

  MODEL_LOCK

  free(modelptr);

  _modelList[idx].next = _modelAvail;
  _modelList[idx].ptr  = 0;
  _modelAvail          = &_modelList[idx];

  MODEL_UNLOCK

  if ( MODEL_Debug )
    Message("Removed idx %d from model list", idx);
}

int modelDef(int instID, int modelgribID, const char *name);

static
void model_defaults(void)
{
  int instID;

  instID  = institutInq(  0,   0, "ECMWF", NULL);
  /* (void)    modelDef(instID, 131, "ERA15"); */
  /* (void)    modelDef(instID, 199, "ERA40"); */

  instID  = institutInq(  0,   0, "MPIMET", NULL);
  ECHAM5  = modelDef(instID,  64, "ECHAM5.4");
  (void)    modelDef(instID,  63, "ECHAM5.3");
  (void)    modelDef(instID,  62, "ECHAM5.2");
  (void)    modelDef(instID,  61, "ECHAM5.1");

  instID  = institutInq( 98, 255, "MPIMET", NULL);
  (void)    modelDef(instID,  60, "ECHAM5.0");
  ECHAM4  = modelDef(instID,  50, "ECHAM4");
  (void)    modelDef(instID, 110, "MPIOM1");

  instID  = institutInq(  0,   0, "DWD", NULL);
  (void)    modelDef(instID, 149, "GME");

  instID  = institutInq(  0,   0, "MCH", NULL);
  //(void)  = modelDef(instID, 137, "COSMO");
  COSMO   = modelDef(instID, 255, "COSMO");

  instID  = institutInq(  0,   1, "NCEP", NULL);
  (void)    modelDef(instID,  80, "T62L28MRF");
}

static
void model_initialize(void)
{
  char *env;

#if  defined  (HAVE_LIBPTHREAD)
  /* initialize global API mutex lock */
  pthread_mutex_init(&_model_mutex, NULL);
#endif

  env = getenv("MODEL_DEBUG");
  if ( env ) MODEL_Debug = atoi(env);

  model_list_new();
  atexit(model_list_delete);

  MODEL_LOCK

  model_init_pointer();

  MODEL_UNLOCK

  _model_init = TRUE;

  model_defaults();
}

static
void model_check_ptr(const char *caller, model_t *modelptr)
{
  if ( modelptr == NULL )
    Errorc("model undefined!");
}

int modelSize(void)
{
  int modelsize = 0;
  int i;
  
  MODEL_INIT

  MODEL_LOCK

  for ( i = 0; i < _model_max; i++ )
    if ( _modelList[i].ptr ) modelsize++;

  MODEL_UNLOCK

  return (modelsize);
}


int modelInq(int instID, int modelgribID, char *name)
{
  int modelID = UNDEFID;
  size_t len;
  int found;
  int model_size;
  model_t *modelptr;

  MODEL_INIT

  model_size = modelSize();

  for( modelID = 0; modelID < model_size; modelID++ )
    {
      modelptr = model_to_pointer(modelID);

      if ( modelptr->used )
	{
	  if ( name )
	    {
	      found = 1;
	      if ( instID      != -1 && modelptr->instID      != instID )      found = 0;
	      if ( modelgribID !=  0 && modelptr->modelgribID != modelgribID ) found = 0;

	      if ( found )
		{
		  if ( modelptr->name )
		    {
		      len = strlen(modelptr->name);
		      if ( memcmp(modelptr->name, name, len) == 0 ) break;
		      len = strlen(name);
		      if ( memcmp(modelptr->name, name, len) == 0 ) break;
		    }
		}
	    }
	  else
	    {
	      if ( modelptr->instID      == instID &&
		   modelptr->modelgribID == modelgribID ) break;
	    }
	}
    }

  if ( modelID == model_size ) modelID = UNDEFID;

  return (modelID);
}


int modelDef(int instID, int modelgribID, const char *name)
{
  int modelID = UNDEFID;
  model_t *modelptr;

  MODEL_INIT

  /*
  modelID = modelInq(instID, modelgribID, name);
  */
  if ( modelID == UNDEFID )
    {
      modelptr = model_new_entry();
      if ( ! modelptr ) Error("No memory");

      modelID = modelptr->self;

      modelptr->instID      = instID;
      modelptr->modelgribID = modelgribID;

      if ( name ) modelptr->name = strdupx(name);
    }

  return (modelID);
}


int modelInqInstitut(int modelID)
{
  int instID = UNDEFID;
  model_t *modelptr;

  MODEL_INIT

  if ( modelID != UNDEFID )
    {
      modelptr = model_to_pointer(modelID);

      model_check_ptr(__func__, modelptr);
  
      instID = modelptr->instID;
    }

  return (instID);
}


int modelInqGribID(int modelID)
{
  int modelgribID = UNDEFID;
  model_t *modelptr;

  MODEL_INIT

  if ( modelID != UNDEFID )
    {
      modelptr = model_to_pointer(modelID);

      model_check_ptr(__func__, modelptr);
  
      modelgribID = modelptr->modelgribID;
    }

  return (modelgribID);
}


char *modelInqNamePtr(int modelID)
{
  char *name = NULL;
  model_t *modelptr;

  MODEL_INIT

  if ( modelID != UNDEFID )
    {
      modelptr = model_to_pointer(modelID);

      model_check_ptr(__func__, modelptr);
  
      if ( modelptr->name )
	name = modelptr->name;
    }

  return (name);
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
