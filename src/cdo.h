/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef _CDO_H
#define _CDO_H

#include <stdio.h>
#include <stdbool.h>
#include "dmemory.h"
#include "util.h"
#include "text.h"

/* dummy use of unused parameters to silence compiler warnings */
#define  UNUSED(x) (void)x

#undef   TRUE
#define  TRUE   1
#undef   FALSE
#define  FALSE  0

#undef   MIN
#define  MIN(a,b)  ((a) < (b) ? (a) : (b))
#undef   MAX
#define  MAX(a,b)  ((a) > (b) ? (a) : (b))

#define  ADD_PLURAL(n)  ((n)>1 ? "s" : "")

#define  UNCHANGED_RECORD  (processSelf() == 0 && cdoStreamName(0)->argv[0][0] != '-' && cdoRegulargrid == FALSE && cdoDefaultFileType == -1 && cdoDefaultDataType == -1 && cdoDefaultByteorder == -1 )


extern int ompNumThreads;

extern int stdin_is_tty;
extern int stdout_is_tty;
extern int stderr_is_tty;

extern int cdoDefaultFileType;
extern int cdoDefaultDataType;
extern int cdoDefaultByteorder;
extern int cdoDefaultTableID;
extern int cdoDefaultInstID;

extern int cdoLockIO;
extern int cdoCheckDatarange;

extern int cdoSilentMode;
extern int cdoOverwriteMode;
extern int cdoRegulargrid;
extern int cdoBenchmark;
extern int cdoTimer;
extern int cdoVerbose;
extern int cdoDebug;
extern int cdoCompress;
extern int cdoInteractive;
extern int cdoParIO;

extern int cdoCompType;
extern int cdoCompLevel;

extern int cdoChunkType;

extern int cdoExpMode;

extern int CDO_Color;
extern int CDO_Use_FFTW;
extern int cdoDiag;

extern int cdoNumVarnames;
extern char **cdoVarnames;

int cdo_omp_get_thread_num(void);



#endif  /* _CDO_H */
