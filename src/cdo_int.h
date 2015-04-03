/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef _CDO_INT_H
#define _CDO_INT_H

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "timebase.h"
#include "field.h"
#include "functs.h"
#include "dmemory.h"
#include "process.h"
#include "const.h"

#define  OPENMP4  201307

#ifndef strdupx
#ifndef strdup
char *strdup(const char *s);
#endif
#define strdupx  strdup
/*
#define strdupx(s)			          \
({					      	  \
   const char *__old = (s);			  \
   size_t __len = strlen(__old) + 1;		  \
   char *__new = malloc(__len);	  \
   (char *) memcpy(__new, __old, __len);	  \
})
*/
#endif

#define strcompare(s1, s2)  (strncmp(s1, s2, strlen(s2)))


/* sxxxYYYYMMDDhhmm0 */
#define  DATE_LEN  31        /* YYYYMMDDhhmmss allocate DTLEN+1 !!!! */
#define  SET_DATE(dtstr, date, time)      (sprintf(dtstr, "%*d%*d", DATE_LEN-6, date, 6, time))
#define  DATE_IS_NEQ(dtstr1, dtstr2, len) (memcmp(dtstr1, dtstr2, len) != 0)


#if defined(__xlC__) /* performance problems on IBM */
#ifndef DBL_IS_NAN
#  define DBL_IS_NAN(x)     ((x) != (x))
#endif
#else
#ifndef DBL_IS_NAN
#if defined(HAVE_DECL_ISNAN)
#  define DBL_IS_NAN(x)     (isnan(x))
#elif defined(FP_NAN)
#  define DBL_IS_NAN(x)     (fpclassify(x) == FP_NAN)
#else
#  define DBL_IS_NAN(x)     ((x) != (x))
#endif
#endif
#endif

#ifndef DBL_IS_EQUAL
/*#define DBL_IS_EQUAL(x,y) (!(x < y || y < x)) */
#  define DBL_IS_EQUAL(x,y) (DBL_IS_NAN(x)||DBL_IS_NAN(y)?(DBL_IS_NAN(x)&&DBL_IS_NAN(y)?1:0):!(x < y || y < x))
#endif

#ifndef IS_EQUAL
#  define IS_NOT_EQUAL(x,y) (x < y || y < x)
#  define IS_EQUAL(x,y)     (!IS_NOT_EQUAL(x,y))
#endif


#ifndef  M_LN10
#define  M_LN10      2.30258509299404568402  /* log_e 10 */
#endif

#ifndef  M_PI
#define  M_PI        3.14159265358979323846  /* pi */
#endif


#define  IX2D(y,x,nx)  ((y)*(nx)+(x))

#define  MEMTYPE_DOUBLE  1
#define  MEMTYPE_FLOAT   2

#define  CDO_EXP_LOCAL   1
#define  CDO_EXP_REMOTE  2


enum {DATE_FIRST, DATE_LAST, DATE_MIDDLE};

void strtolower(char *str);

void print_pthread_info(void);

void cdoProcessTime(double *utime, double *stime);

void    setCommandLine(int argc, char **argv);
char   *commandLine(void);
int     readline(FILE *fp, char *line, int len);

int zaxis2ltype(int zaxisID);


int nfc2nlat(int nfc, int ntr);
int nlat2ntr(int nlat);
int nlat2ntr_linear(int nlat);
int ntr2nlat(int ntr);
int ntr2nlat_linear(int ntr);
int compNlon(int nlat);

void param2str(int param, char *paramstr, int maxlen);
void date2str(int date, char *datestr, int maxlen);
void time2str(int time, char *timestr, int maxlen);

const char * tunit2str(int tunits);
const char * calendar2str(int calendar);


typedef struct {
  int   date;
  int   time;
} datetime_t;

typedef struct
{
  datetime_t v;
  datetime_t b[2];
} dtinfo_t;

typedef struct {
  int   julday;
  int   secofday;
} juldate_t;


juldate_t juldate_encode(int calendar, int date, int time);
void      juldate_decode(int calendar, juldate_t juldate, int *date, int *time);
juldate_t juldate_sub(juldate_t juldate2, juldate_t juldate1);
juldate_t juldate_add_seconds(int seconds, juldate_t juldate);
double    juldate_to_seconds(juldate_t juldate);

void    get_timestat_date(int *tstat_date);
void    datetime_avg(int dpy, int ndates, datetime_t *datetime);
void    datetime_avg_dtinfo(int dpy, int ndates, dtinfo_t *dtinfo);
void    taxisInqDTinfo(int taxisID, dtinfo_t *dtinfo);
void    taxisDefDTinfo(int taxisID, dtinfo_t dtinfo);

int     days_per_month(int calendar, int year, int month);
int     days_per_year(int calendar, int year);
int     calendar_dpy(int calendar);

void    defineGrid(const char *gridarg);
void    defineInstitution(char *instarg);
int     defineTable(char *tablearg);

void    cdolog(const char *prompt, double cputime);
void    cdologs(int noper);
void    cdologo(int noper);
void    nospec(int vlistID);
void    gridWrite(FILE *fp, int gridID);

void openLock(void);
void openUnlock(void);

int  cdf_openread(const char *filename);

void printFiletype(int streamID, int vlistID);

void job_submit(const char *expname, const char *jobfilename, const char *jobname, const char *tmppath, const char *ftppath);

void minmaxval(long nvals, double *array, int *imiss, double *minval, double *maxval);

off_t filesize(const char *filename);


#endif  /* _CDO_INT_H */
