/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2009 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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
#include "dmemory.h"

#undef   TRUE
#define  TRUE   1
#undef   FALSE
#define  FALSE  0

#undef   MIN
#define  MIN(a,b)  ((a) < (b) ? (a) : (b))
#undef   MAX
#define  MAX(a,b)  ((a) > (b) ? (a) : (b))
#undef   NINT
#define  NINT(x)   ((x) < 0 ? (int)((x)-0.5) : (int)((x)+0.5))

#define  UNCHANGED_RECORD  (processSelf() == 0 && *cdoStreamName(0) != '-' && cdoRegulargrid == FALSE && cdoDefaultDataType == -1 && cdoDefaultByteorder == -1 )


typedef struct {
  int      argc;
  char   **argv;
}
ARGUMENT;


extern int ompNumThreads;

extern int cdoDefaultFileType;
extern int cdoDefaultDataType;
extern int cdoDefaultByteorder;
extern int cdoDefaultTableID;
extern int cdoDefaultInstID;

extern int cdoSilentMode;
extern int cdoRegulargrid;
extern int cdoBenchmark;
extern int cdoTimer;
extern int cdoVerbose;
extern int cdoDebug;
extern int cdoCompress;
extern int cdoInteractive;
extern int cdoParIO;

extern int cdoZtype;
extern int cdoZlevel;

extern int cdoExpMode;

extern int cdoDisableFilesuffix;
extern int cdoDiag;

void    cdiError(int cdiErrno, const char *fmt, ...);
void    cdoAbort(const char *fmt, ...);
void    cdoWarning(const char *fmt, ...);
void    cdoPrint(const char *fmt, ...);

int  timer_new(char *text);
void timer_report(void);
void timer_start(int it);
void timer_stop(int it);

void    timerStart(int timer);
void    timerStop(int timer);
void    timerClear(int timer);
void    timerPrint(int timer);
void    timersPrint(void);
void    timersInit(void);

void    operatorInputArg(const char *enter);
int     operatorArgc(void);
char  **operatorArgv(void);
void    operatorCheckArgc(int numargs);

const char *cdoStreamName(int cnt);

void    cdoInitialize(void *argument);
void    cdoFinish(void);

int     cdoStreamCnt(void);
int     cdoOperatorAdd(const char *name, int func, int intval, const char *enter);
int     cdoOperatorID(void);
int     cdoOperatorFunc(int operID);
int     cdoOperatorIntval(int operID);
const char *cdoOperatorName(int operID);
const char *cdoOperatorEnter(int operID);

int     cdoFiletype(void);

void    cdoInqHistory(int fileID);
void    cdoDefHistory(int fileID, char *histstring);

int     cdoDefineGrid(const char *gridfile);
int     cdoDefineZaxis(const char *zaxisfile);

int     vlistIsSzipped(int vlistID);
void    vlistCompare(int vlistID1, int vlistID2, int function);

int  gridWeights(int gridID, double *weights);
int  gridGenArea(int gridID, double *area);
void gaussaw(double pa[], double pw[], int nlat);
void genXbounds(long xsize, long ysize, const double *grid_center_lon, double *grid_corner_lon, double dlon);
void genYbounds(long xsize, long ysize, const double *grid_center_lat, double *grid_corner_lat);
void writeNCgrid(const char *gridfile, int gridID, int *imask);
void defineZaxis(const char *zaxisarg);
void cdiDefTableID(int tableID);
void gridGenXvals(int xsize, double xfirst, double xlast, double xinc, double *xvals);
void gridGenYvals(int gridtype, int ysize, double yfirst, double ylast, double yinc, double *yvals);

int gridFromName(const char *gridname);

#endif  /* _CDO_H */
