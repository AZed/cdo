/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2010 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef _PIPE_H
#define _PIPE_H

#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <sys/types.h>

#if  defined  (HAVE_LIBPTHREAD)

#include <pthread.h>
#include <pthread_debug.h>

#endif

typedef struct {
  short       check_datarange;
  int         gridsize;
  int         datatype;
  double      missval;
  double      addoffset;
  double      scalefactor;
} varlist_t;


struct _PSTREAM {
  int            self;
  int            mode;
  int            fileID;
  int            vlistID;
  int            tsID;
  int            filetype;
  int            ispipe;
  int            isopen;
  int            tsID0;
  int            mfiles;
  int            nfiles;
  int            varID;           /* next varID defined with streamDefVar */
  char          *name;
  char         **mfnames;
  varlist_t     *varlist;
#if  defined  (HAVE_LIBPTHREAD)
  struct _PIPE  *pipe;
  pthread_t     rthreadID; /* read  thread ID */
  pthread_t     wthreadID; /* write thread ID */
#endif
};

typedef struct _PSTREAM PSTREAM;


#if  defined  (HAVE_LIBPTHREAD)

struct _PIPE {
  int     nrecs, EOP;
  int     varID, levelID;
  int     recIDr, recIDw, tsIDr, tsIDw;
  int     hasdata, usedata;
  int     nmiss;
  double *data;
  struct _PSTREAM *pstreamptr_in;
  /* unsigned long */ off_t nvals;
  pthread_mutex_t *mutex;
  pthread_cond_t *tsDef, *tsInq, *vlistDef, *isclosed;
  pthread_cond_t *recDef, *recInq;
  pthread_cond_t *writeCond, *readCond;
};

typedef struct _PIPE PIPE;

PIPE *pipeNew(void);
void  pipeDelete(PIPE *pipe);

void  pipeDebug(int debug);

void  pipeDefVlist(PSTREAM *pstreamptr, int vlistID);
int   pipeInqVlist(PSTREAM *pstreamptr);

void  pipeDefTimestep(PSTREAM *pstreamptr, int tsID);
int   pipeInqTimestep(PSTREAM *pstreamptr, int tsID);

void  pipeDefRecord(PSTREAM *pstreamptr, int  varID, int  levelID);
int   pipeInqRecord(PSTREAM *pstreamptr, int *varID, int *levelID);

void  pipeReadRecord(PSTREAM *pstreamptr, double *data, int *nmiss);
void  pipeWriteRecord(PSTREAM *pstreamptr, double *data, int nmiss);
void  pipeCopyRecord(PSTREAM *pstreamptr_dest, PSTREAM *pstreamptr_src);

#endif

#endif  /* _PIPE_H */
