/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/


#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include "pipe.h"
#include "pstream_int.h"
#include "cdi.h"
#include "cdo.h"
#include "error.h"
#include "dmemory.h"

#if  defined  (HAVE_LIBPTHREAD)

static int PipeDebug = 0;

static void pipe_init(PIPE *pipe)
{
  static char func[] = "pipe_init";
  pthread_mutexattr_t m_attr;
  pthread_condattr_t c_attr;

  pthread_mutexattr_init(&m_attr);
  pthread_condattr_init(&c_attr);
  /*
#if defined (_POSIX_THREAD_PROCESS_SHARED)
  if ( PipeDebug )
    {
      Message(func, "setpshared mutexattr to PTHREAD_PROCESS_SHARED");
      Message(func, "setpshared condattr to PTHREAD_PROCESS_SHARED");
    }

  pthread_mutexattr_setpshared(&m_attr, PTHREAD_PROCESS_SHARED);
  pthread_condattr_setpshared(&c_attr, PTHREAD_PROCESS_SHARED);

  if ( PipeDebug )
    {
      int pshared;
      pthread_mutexattr_getpshared(&m_attr, &pshared);
      if ( pshared == PTHREAD_PROCESS_SHARED )
	Message(func, "getpshared mutexattr is PTHREAD_PROCESS_SHARED");
      else if ( pshared == PTHREAD_PROCESS_PRIVATE )
	Message(func, "getpshared mutexattr is PTHREAD_PROCESS_PRIVATE");

      pthread_condattr_getpshared(&c_attr, &pshared);
      if ( pshared == PTHREAD_PROCESS_SHARED )
	Message(func, "getpshared condattr is PTHREAD_PROCESS_SHARED");
      else if ( pshared == PTHREAD_PROCESS_PRIVATE )
	Message(func, "getpshared condattr is PTHREAD_PROCESS_PRIVATE");
    }
#else
  if ( PipeDebug )
    Message(func, "_POSIX_THREAD_PROCESS_SHARED undefined");
#endif
  */
  pipe->EOP     = 0;

  pipe->recIDr  = -1;
  pipe->recIDw  = -1;
  pipe->tsIDr   = -1;
  pipe->tsIDw   = -1;

  pipe->nvals   = 0;
  pipe->nmiss   = 0;
  pipe->data    = NULL;
  pipe->hasdata = 0;
  pipe->usedata = TRUE;
  pipe->pstreamptr_in = 0;

  pipe->mutex = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(pipe->mutex, &m_attr);

  pipe->tsDef = (pthread_cond_t *) malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->tsDef, &c_attr);
  pipe->tsInq = (pthread_cond_t *) malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->tsInq, &c_attr);

  pipe->recDef = (pthread_cond_t *) malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->recDef, &c_attr);
  pipe->recInq = (pthread_cond_t *) malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->recInq, &c_attr);
  
  pipe->vlistDef = (pthread_cond_t *) malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->vlistDef, &c_attr);
  pipe->isclosed = (pthread_cond_t *) malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->isclosed, &c_attr);

  pipe->writeCond = (pthread_cond_t *) malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->writeCond, &c_attr);

  pipe->readCond = (pthread_cond_t *) malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->readCond, &c_attr);

  pthread_mutexattr_destroy(&m_attr);
  pthread_condattr_destroy(&c_attr);
}


PIPE *pipeNew()
{
  static char func[] = "pipeNew";  
  PIPE *pipe;

  pipe = (PIPE *) malloc(sizeof(PIPE));

  pipe_init(pipe);

  return (pipe);
}


void pipeDelete(PIPE *pipe)
{
  static char func[] = "pipeDelete";  

  if ( pipe )
    {
      if ( pipe->mutex )     free(pipe->mutex);
      if ( pipe->tsDef )     free(pipe->tsDef);
      if ( pipe->tsInq )     free(pipe->tsInq);
      if ( pipe->recDef )    free(pipe->recDef);
      if ( pipe->recInq )    free(pipe->recInq);
      if ( pipe->vlistDef )  free(pipe->vlistDef);
      if ( pipe->isclosed )  free(pipe->isclosed);
      if ( pipe->writeCond ) free(pipe->writeCond);
      if ( pipe->readCond )  free(pipe->readCond);
      free(pipe);
    }
}


void pipeDefVlist(PSTREAM *pstreamptr, int vlistID)
{
  static char func[] = "pipeDefVlist";
  char *pname = pstreamptr->name;
  PIPE *pipe;

  if ( PipeDebug )
    Message(func, "%s pstreamID %d", pname, pstreamptr->self);

  pipe = pstreamptr->pipe;

  pthread_mutex_lock(pipe->mutex);
  pstreamptr->vlistID = vlistID;
  pthread_mutex_unlock(pipe->mutex);  

  pthread_cond_signal(pipe->vlistDef);
}


int pipeInqVlist(PSTREAM *pstreamptr)
{
  static char func[] = "pipeInqVlist";
  char *pname = pstreamptr->name;
  PIPE *pipe;
  int vlistID;

  if ( PipeDebug )
    Message(func, "%s pstreamID %d", pname, pstreamptr->self);

  pipe = pstreamptr->pipe;

  pthread_mutex_lock(pipe->mutex);
  while ( pstreamptr->vlistID == -1 )
    {
      if ( PipeDebug ) Message(func, "%s wait of vlist", pname);
      pthread_cond_wait(pipe->vlistDef, pipe->mutex);
    }
  vlistID = pstreamptr->vlistID;
  pthread_mutex_unlock(pipe->mutex);

  return (vlistID);
}


int pipeInqTimestep(PSTREAM *pstreamptr, int tsID)
{
  static char func[] = "pipeInqTimestep";
  char *pname = pstreamptr->name;
  PIPE *pipe;
  int nrecs;

  if ( PipeDebug )
    Message(func, "%s pstreamID %d", pname, pstreamptr->self);

  pipe = pstreamptr->pipe;

  pthread_mutex_lock(pipe->mutex);
  pipe->usedata = FALSE;
  pipe->recIDr = -1;
  if ( tsID != pipe->tsIDr+1 )
    {
      if ( ! (tsID == pipe->tsIDr && pipe->tsIDr == pipe->tsIDw && pipe->recIDr == -1) ) 
	Error(func, "%s unexpected tsID %d %d %d", pname, tsID, pipe->tsIDr+1, pipe->tsIDw);
    }
  
  pipe->tsIDr = tsID;
  while ( pipe->tsIDw != tsID )
    {
      if ( pipe->EOP == TRUE )
	{
	  if ( PipeDebug ) Message(func, "%s EOP", pname);
	  break;
	}
      if ( pipe->hasdata )
	{
	  if ( PipeDebug ) Message(func, "%s has data", pname);
	  pipe->hasdata = 0;
	  pipe->data = NULL;
	  pthread_cond_signal(pipe->readCond);
	}
      else
	if ( PipeDebug ) Message(func, "%s has no data", pname);

      pthread_cond_signal(pipe->recInq); /* o.k. ??? */

      if ( PipeDebug ) Message(func, "%s wait of tsDef", pname);
      pthread_cond_wait(pipe->tsDef, pipe->mutex);
    }

  if ( pipe->EOP == TRUE )
    nrecs = 0;
  else
    nrecs = pipe->nrecs;

  pthread_mutex_unlock(pipe->mutex);

  pthread_cond_signal(pipe->tsInq);

  return (nrecs);
}


void pipeDefTimestep(PSTREAM *pstreamptr, int tsID)
{
  static char func[] = "pipeDefTimestep";
  char *pname = pstreamptr->name;
  int nrecs;
  PIPE *pipe;

  if ( PipeDebug )
    Message(func, "%s pstreamID %d", pname, pstreamptr->self);

  pipe = pstreamptr->pipe;

  pthread_mutex_lock(pipe->mutex);
  pipe->recIDw = -1;
  pipe->tsIDw++;
  if ( tsID != pipe->tsIDw )
    Error(func, "unexpected tsID %d(%d) for %s", tsID, pipe->tsIDw, pname);

  if ( tsID == 0 )
    nrecs = vlistNrecs(pstreamptr->vlistID);
  else
    {
      int vlistID, varID;
      vlistID = pstreamptr->vlistID;
      nrecs = 0;
      for ( varID = 0; varID < vlistNvars(vlistID); varID++ )
	if ( vlistInqVarTime(vlistID, varID) == TIME_VARIABLE )
	  nrecs += zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
    }

  pipe->nrecs = nrecs;
  pthread_mutex_unlock(pipe->mutex);
  pthread_cond_signal(pipe->tsDef);
  /*
sleep(1);
*/
  pthread_mutex_lock(pipe->mutex);
  while ( pipe->tsIDr < tsID )
    {
      if ( pipe->EOP == TRUE )
	{
	  if ( PipeDebug ) Message(func, "EOP");
	  break;
	}
      if ( PipeDebug ) Message(func, "%s wait of tsInq", pname);
      pthread_cond_wait(pipe->tsInq, pipe->mutex);
    }
  pthread_mutex_unlock(pipe->mutex);
}


int pipeInqRecord(PSTREAM *pstreamptr, int *varID, int *levelID)
{
  static char func[] = "pipeInqRecord";
  char *pname = pstreamptr->name;
  PIPE *pipe;
  int condSignal = FALSE;

  if ( PipeDebug )
    Message(func, "%s pstreamID %d", pname, pstreamptr->self);
 
  pipe = pstreamptr->pipe;

  pthread_mutex_lock(pipe->mutex);

  if ( PipeDebug )
    Message(func, "%s has no data %d %d", pname, pipe->recIDr, pipe->recIDw);

  if ( pipe->hasdata || pipe->usedata )
    {
      pipe->hasdata = 0;
      pipe->data = NULL;
      pipe->usedata = FALSE;
      condSignal = TRUE;
    }
 	  
  pthread_mutex_unlock(pipe->mutex);
  if ( condSignal ) pthread_cond_signal(pipe->readCond);

  pthread_mutex_lock(pipe->mutex);

  pipe->usedata = TRUE;

  pipe->recIDr++;
  
  if ( PipeDebug )
    Message(func, "%s recID %d %d", pname, pipe->recIDr, pipe->recIDw);  

  while ( pipe->recIDw != pipe->recIDr )
    {
      if ( pipe->EOP == TRUE )
	{
	  if ( PipeDebug ) Message(func, "EOP");
	  break;
	}
      if ( PipeDebug ) Message(func, "%s wait of recDef", pname);
      pthread_cond_wait(pipe->recDef, pipe->mutex);
    }

  if ( pipe->EOP == TRUE )
    {
      *varID   = -1;
      *levelID = -1;
    }
  else
    {
      *varID   = pipe->varID;
      *levelID = pipe->levelID;
    }

  pthread_mutex_unlock(pipe->mutex);

  pthread_cond_signal(pipe->recInq);

  return (0);
}


void pipeDefRecord(PSTREAM *pstreamptr, int varID, int levelID)
{
  static char func[] = "pipeDefRecord";
  char *pname = pstreamptr->name;
  PIPE *pipe;
  int condSignal = FALSE;

  if ( PipeDebug )
    Message(func, "%s pstreamID %d", pname, pstreamptr->self);

  pipe = pstreamptr->pipe;

  pthread_mutex_lock(pipe->mutex);

  if ( PipeDebug )
    Message(func, "%s has data %d %d", pname, pipe->recIDr, pipe->recIDw);
	  
  if ( pipe->hasdata )
    {
      pipe->hasdata = 0;
      pipe->data = NULL;
      condSignal = TRUE;
    }

  pthread_mutex_unlock(pipe->mutex);
  if ( condSignal ) pthread_cond_signal(pipe->readCond);

  pthread_mutex_lock(pipe->mutex);
  pipe->usedata = TRUE;
  pipe->recIDw++;

  if ( PipeDebug )
    Message(func, "%s recID %d %d", pname, pipe->recIDr, pipe->recIDw);  

  pipe->varID = varID;
  pipe->levelID = levelID;
  pthread_mutex_unlock(pipe->mutex);
  pthread_cond_signal(pipe->recDef);

  pthread_mutex_lock(pipe->mutex);
  while ( pipe->recIDr < pipe->recIDw )
    {
      if ( pipe->tsIDw != pipe->tsIDr ) break;
      if ( pipe->EOP == TRUE ) break;
      if ( PipeDebug ) Message(func, "%s wait of recInq %d", pname, pipe->recIDr);
      pthread_cond_wait(pipe->recInq, pipe->mutex);
    }
  pthread_mutex_unlock(pipe->mutex);
}


void pipeCopyRecord(PSTREAM *pstreamptr_out, PSTREAM *pstreamptr_in)
{
  static char func[] = "pipeCopyRecord";
  char *ipname = pstreamptr_in->name;
  char *opname = pstreamptr_out->name;
  PIPE *pipe;

  if ( PipeDebug )
    Message(func, "%s pstreamIDin %d", ipname, pstreamptr_in->self);
  if ( PipeDebug )
    Message(func, "%s pstreamIDout %d", opname, pstreamptr_out->self);

  pipe = pstreamptr_out->pipe;

  pthread_mutex_lock(pipe->mutex);
  pipe->hasdata = 2; /* pipe */
  pipe->pstreamptr_in = pstreamptr_in;
  pthread_mutex_unlock(pipe->mutex);

  pthread_cond_signal(pipe->writeCond);

  pthread_mutex_lock(pipe->mutex);
  while ( pipe->hasdata )
    {
      if ( pipe->usedata == FALSE ) break;

      if ( pipe->recIDw != pipe->recIDr ) break;

      if ( pipe->EOP == TRUE )
	{
	  if ( PipeDebug ) Message(func, "EOP");
	  break;
	}
      if ( PipeDebug ) Message(func, "%s wait of read record", opname);
      pthread_cond_wait(pipe->readCond, pipe->mutex);
    }
  pthread_mutex_unlock(pipe->mutex);
}


void pipeReadRecord(PSTREAM *pstreamptr, double *data, int *nmiss)
{
  static char func[] = "pipeReadRecord";
  char *pname = pstreamptr->name;
  PIPE *pipe;

  *nmiss = 0;
  if ( PipeDebug )
    Message(func, "%s pstreamID %d", pname, pstreamptr->self);

  pipe = pstreamptr->pipe;

  pthread_mutex_lock(pipe->mutex);
  while ( pipe->hasdata == 0 )
    {
      if ( PipeDebug ) Message(func, "%s wait of writeCond", pname);
      pthread_cond_wait(pipe->writeCond, pipe->mutex);
    }

  if ( pipe->hasdata == 2 )
    {
      PSTREAM *pstreamptr_in;

      pstreamptr_in = pipe->pstreamptr_in;
      /*
      if ( pstreamptr_in->ispipe )
	{
	  pstreamptr_in = pstreamptr_in->pipe->pstreamptr_in;
	  if ( pstreamptr_in->ispipe )
	    Error(func, "istream is pipe");
	}
      */
      pstreamptr = pstreamptr_in;
      while ( pstreamptr_in->ispipe )
	{
	  if ( PipeDebug ) fprintf(stderr, "%s: istream %d is pipe\n", func, pstreamptr_in->self);
	  pstreamptr    = pstreamptr_in;
	  pstreamptr_in = pstreamptr_in->pipe->pstreamptr_in;
	  if ( pstreamptr_in == 0 ) break;
	}

      if ( pstreamptr_in == 0 )
	{
	  if ( PipeDebug ) fprintf(stderr, "pstreamID = %d\n", pstreamptr->self);
	  if ( pstreamptr->pipe->hasdata == 1 )
	    {
	      int vlistID, gridsize;

	      if ( ! pstreamptr->pipe->data )
		Error(func, "No data pointer for %s", pname);

	      vlistID = pstreamptr->vlistID;
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID, pstreamptr->pipe->varID));
	      memcpy(data, pstreamptr->pipe->data, gridsize*sizeof(double));
	      *nmiss = pstreamptr->pipe->nmiss;
	    }
	  else
	    Error(func, "Internal problem! istream undefined");
	}
      else
	{
	  if ( PipeDebug ) fprintf(stderr, "%s: istream %d is file\n", func, pstreamptr_in->self);
	  streamReadRecord(pstreamptr_in->fileID, data, nmiss);
	}
    }
  else if ( pipe->hasdata == 1 )
    {
      int vlistID, gridsize;

      if ( ! pipe->data )
	Error(func, "No data pointer for %s", pname);

      vlistID = pstreamptr->vlistID;
      gridsize = gridInqSize(vlistInqVarGrid(vlistID, pipe->varID));
      pipe->nvals += gridsize;
      memcpy(data, pipe->data, gridsize*sizeof(double));
      *nmiss = pipe->nmiss;
    }
  else
    {
      Error(func, "data type %d not implemented", pipe->hasdata);
    }

  if ( PipeDebug ) Message(func, "%s read record %d", pname, pipe->recIDr);
 

  pipe->hasdata = 0;
  pipe->data = NULL;
  pthread_mutex_unlock(pipe->mutex);

  pthread_cond_signal(pipe->readCond);
}


void pipeWriteRecord(PSTREAM *pstreamptr, double *data, int nmiss)
{
  static char func[] = "pipeWriteRecord";
  char *pname = pstreamptr->name;
  PIPE *pipe;

  if ( PipeDebug )
    Message(func, "%s pstreamID %d", pname, pstreamptr->self);

  pipe = pstreamptr->pipe;
  /*
  if ( ! pipe->usedata ) return;
  */
  pthread_mutex_lock(pipe->mutex);
  pipe->hasdata = 1; /* data pointer */
  pipe->data    = data;
  pipe->nmiss   = nmiss;
  pthread_mutex_unlock(pipe->mutex);
  pthread_cond_signal(pipe->writeCond);

  if ( PipeDebug ) Message(func, "%s write record %d", pname, pipe->recIDw);

  pthread_mutex_lock(pipe->mutex);
  while ( pipe->hasdata )
    {
      if ( pipe->usedata == FALSE ) break;
      /*     
      printf("ts ids %d %d\n", pipe->tsIDw, pipe->tsIDr);
      printf("rec ids %d %d\n", pipe->recIDw, pipe->recIDr);
      */
      if ( pipe->recIDw != pipe->recIDr ) break;

      if ( pipe->EOP == TRUE )
	{
	  if ( PipeDebug ) Message(func, "EOP");
	  break;
	}
      if ( PipeDebug ) Message(func, "%s wait of readCond", pname);
      pthread_cond_wait(pipe->readCond, pipe->mutex);
    }
  pthread_mutex_unlock(pipe->mutex);
}


void pipeDebug(int debug)
{
  PipeDebug = debug;
}

#endif
