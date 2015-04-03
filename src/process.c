/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

#if  defined  (HAVE_PTHREAD_H)
#  include <pthread.h>
#endif

#include <stdio.h>
#include <string.h>

#include "cdo.h"
#include "cdo_int.h"
#include "error.h"
#include "modules.h"
#include "util.h"
#include "pstream_int.h"
#include "dmemory.h"


#define  MAX_PROCESS   128
#define  MAX_STREAM     64
#define  MAX_OPERATOR  128
#define  MAX_ARGC     4096


typedef struct {
  int         f1;
  int         f2;
  const char *name;
  const char *enter;
}
operator_t;

typedef struct {
#if  defined  (HAVE_LIBPTHREAD)
  pthread_t threadID;
#endif
  short      nchild;
  short      nstream;
  short      streams[MAX_STREAM];
  double     s_utime;
  double     s_stime;
  double     a_utime;
  double     a_stime;
  double     cputime;

  off_t      nvals;
  short      nvars;
  int        ntimesteps;
  short      streamCnt;
  char     **streamNames;
  char      *xoperator;
  char      *operatorName;
  char      *operatorArg;
  int        oargc;
  char      *oargv[MAX_ARGC];
  char       prompt[64];
  short      noper;
  operator_t operator[MAX_OPERATOR];
}
process_t;


static process_t Process[MAX_PROCESS];

static int NumProcess = 0;
static int NumProcessActive = 0;

#if  defined  (HAVE_LIBPTHREAD)
pthread_mutex_t processMutex = PTHREAD_MUTEX_INITIALIZER;
#endif


int processCreate(void)
{
  int processID;

#if  defined  (HAVE_LIBPTHREAD)
  pthread_mutex_lock(&processMutex);
#endif
  processID = NumProcess++;
  NumProcessActive++;
#if  defined  (HAVE_LIBPTHREAD)
  pthread_mutex_unlock(&processMutex);  
#endif

  if ( processID >= MAX_PROCESS )
    Error("Limit of %d processes reached!", MAX_PROCESS);

#if  defined  (HAVE_LIBPTHREAD)
  Process[processID].threadID     = pthread_self();
#endif
  Process[processID].nstream      = 0;
  Process[processID].nchild       = 0;

  cdoProcessTime(&Process[processID].s_utime, &Process[processID].s_stime);
  Process[processID].a_utime      = 0;
  Process[processID].a_stime      = 0;
  Process[processID].cputime      = 0;

  Process[processID].oargc        = 0;
  Process[processID].xoperator    = NULL;
  Process[processID].operatorName = NULL;
  Process[processID].operatorArg  = NULL;

  Process[processID].noper        = 0;

  return (processID);
}


int processSelf(void)
{
  int processID = 0;
#if  defined  (HAVE_LIBPTHREAD)
  pthread_t thID = pthread_self();

  pthread_mutex_lock(&processMutex);

  for ( processID = 0; processID < NumProcess; processID++ )
    if ( pthread_equal(Process[processID].threadID, thID) ) break;

  if ( processID == NumProcess )
    {
      if ( NumProcess > 0 )
	Error("Internal problem, process not found!");
      else
	processID = 0;
    }

  pthread_mutex_unlock(&processMutex);  

#endif

  return (processID);
}


int processNums(void)
{
  int pnums = 0;

#if  defined  (HAVE_LIBPTHREAD)
  pthread_mutex_lock(&processMutex);
#endif

  pnums = NumProcess;

#if  defined  (HAVE_LIBPTHREAD)
  pthread_mutex_unlock(&processMutex);  
#endif

  return (pnums);
}


int processNumsActive(void)
{
  int pnums = 0;

#if  defined  (HAVE_LIBPTHREAD)
  pthread_mutex_lock(&processMutex);
#endif

  pnums = NumProcessActive;

#if  defined  (HAVE_LIBPTHREAD)
  pthread_mutex_unlock(&processMutex);  
#endif

  return (pnums);
}


void processAddNvals(off_t nvals)
{
  int processID = processSelf();

  Process[processID].nvals += nvals;
}


off_t processInqNvals(int processID)
{
  off_t nvals = 0;

  nvals = Process[processID].nvals;

  return (nvals);
}


void processAddStream(int streamID)
{
  int processID = processSelf();
  int sindex;

  if ( pstreamIsPipe(streamID) ) Process[processID].nchild++;

  sindex = Process[processID].nstream++;

  if ( sindex >= MAX_STREAM )
    Error("limit of %d streams per process reached!", MAX_STREAM);

  Process[processID].streams[sindex] = streamID;
}


void processDefCputime(int processID, double cputime)
{
  Process[processID].cputime = cputime;
}


double processInqCputime(int processID)
{
  return (Process[processID].cputime);
}


void processStartTime(double *utime, double *stime)
{
  int processID = processSelf();

  *utime = Process[processID].s_utime;
  *stime = Process[processID].s_stime;
}


void processEndTime(double *utime, double *stime)
{
  *utime = Process[0].a_utime;
  *stime = Process[0].a_stime;
}


void processAccuTime(double utime, double stime)
{
  Process[0].a_utime += utime;
  Process[0].a_stime += stime;
}


int processInqStreamNum(void)
{
  int processID = processSelf();

  return (Process[processID].nstream);
}


int processInqChildNum(void)
{
  int processID = processSelf();

  return (Process[processID].nchild);
}


int processInqStreamID(int streamindex)
{
  int processID = processSelf();

  return (Process[processID].streams[streamindex]);
}


const char *processInqOpername2(int processID)
{
  return (Process[processID].operatorName);
}


const char *processInqOpername(void)
{
  int processID = processSelf();

  return (Process[processID].operatorName);
}


void processDefPrompt(char *opername)
{
  int processID = processSelf();
  extern char *Progname;

  if ( processID == 0 )
    sprintf(Process[processID].prompt, "%s %s", Progname, opername);
  else
    sprintf(Process[processID].prompt, "%s(%d) %s", Progname, processID+1, opername);
}


const char *processInqPrompt(void)
{
  int processID = processSelf();

  return (Process[processID].prompt);
}


int cdoStreamCnt(void)
{
  int processID = processSelf();
  int cnt;

  cnt = Process[processID].streamCnt;

  return (cnt);
}


const char *cdoStreamName(int cnt)
{
  int processID = processSelf();

  if ( cnt > Process[processID].streamCnt || cnt < 0 )
    Error("count %d out of range!", cnt);

  return (Process[processID].streamNames[cnt]);
}


const char *processOperator(void)
{
  int processID = processSelf();

  return (Process[processID].xoperator);
}

static
char *getOperatorArg(const char *xoperator)
{
  char *commapos;
  char *operatorArg = NULL;
  size_t len;

  if ( xoperator )
    {
      commapos = strchr(xoperator, ',');

      if ( commapos )
	{
	  len = strlen(commapos+1);
	  if ( len )
	    {
	      operatorArg = (char *) malloc(len+1);
	      strcpy(operatorArg, commapos+1);
	    }
	}
    }

  return (operatorArg);
}

static int skipInputStreams(int argc, char *argv[], int globArgc, int nstreams);

static
int getGlobArgc(int argc, char *argv[], int globArgc)
{
  int streamInCnt;
  int streamOutCnt;
  char *opername;
  char *comma_position;
  const char *caller = processInqPrompt();
  /*
  { int i;
  for ( i = 0; i < argc; i++ )
    printf("%d %d %s\n", globArgc, i, argv[i]);
  }
  */
  opername = &argv[globArgc][1];
  comma_position = strchr(opername, ',');
  if ( comma_position ) *comma_position = 0;

  streamInCnt  = operatorStreamInCnt(opername);
  streamOutCnt = operatorStreamOutCnt(opername);

  if ( streamInCnt == -1 )
    {
      /*
      int i;

      for ( i = globArgc+1; i < argc; i++ )
	if ( argv[i][0] == '-' ) break;

      printf("%d %d %d\n", i, argc, globArgc);
      if ( i < argc )
      */
      streamInCnt = 1;
      /*
      Errorc("Unlimited input streams not allowed in CDO pipes (Operator %s)!", opername);
      */
    }

  if ( streamOutCnt > 1 )
    Errorc("More than one output stream not allowed in CDO pipes (Operator %s)!", opername);

  globArgc++;

  if ( streamInCnt > 0 )
    globArgc = skipInputStreams(argc, argv, globArgc, streamInCnt);
  if ( comma_position ) *comma_position = ',';

  return (globArgc);
}

static
int skipInputStreams(int argc, char *argv[], int globArgc, int nstreams)
{
  const char *caller = processInqPrompt();

  while ( nstreams > 0 )
    {
      if ( globArgc >= argc )
	{
	  Errorc("Too few arguments. Check command line!");
	  break;
	}
      if ( argv[globArgc][0] == '-' )
	{
	  globArgc = getGlobArgc(argc, argv, globArgc);
	}
      else
	globArgc++;

      nstreams--;
    }

  return (globArgc);
}

static
int getStreamCnt(int argc, char *argv[])
{
  int streamCnt = 0;
  int globArgc = 1;

  while ( globArgc < argc )
    {
      if ( argv[globArgc][0] == '-' )
	{
	  globArgc = getGlobArgc(argc, argv, globArgc);
	}
      else
	globArgc++;

      streamCnt++;
    }

  return (streamCnt);
}

static
void setStreamNames(int argc, char *argv[])
{
  int processID = processSelf();
  int i;
  int globArgc = 1;
  int globArgcStart;
  char *streamname;
  size_t len;

  while ( globArgc < argc )
    {
      //     printf("arg %d %d %s\n", argc, globArgc, argv[globArgc]);
      if ( argv[globArgc][0] == '-' )
	{
	  globArgcStart = globArgc;

	  globArgc = getGlobArgc(argc, argv, globArgc);
	  //	  printf("globArgc %d\n", globArgc);
	  len = 0;
	  for ( i = globArgcStart; i < globArgc; i++ ) len += strlen(argv[i]) + 1;
	  streamname = (char *) malloc(len);
	  memcpy(streamname, argv[globArgcStart], len);
	  for ( i = 1; i < (int) len-1; i++ ) if ( streamname[i] == '\0' ) streamname[i] = ' ';
	  Process[processID].streamNames[Process[processID].streamCnt++] = streamname;
	  //	  printf("streamname1: %s\n", streamname);
	}
      else
	{
	  len = strlen(argv[globArgc]) + 1;
	  streamname = (char *) malloc(len);
	  strcpy(streamname, argv[globArgc]);
	  Process[processID].streamNames[Process[processID].streamCnt++] = streamname;
	  //	  printf("streamname2: %s\n", streamname);
	  globArgc++;
	}
    }
}

static
void checkStreamCnt(void)
{
  int processID = processSelf();
  int streamInCnt, streamOutCnt;
  int streamCnt = 0;
  int i, j;
  int obase = FALSE;
  const char *caller = processInqPrompt();

  streamInCnt  = operatorStreamInCnt(Process[processID].operatorName);
  streamOutCnt = operatorStreamOutCnt(Process[processID].operatorName);

  if ( streamOutCnt == -1 )
    {
      streamOutCnt = 1;
      obase = TRUE;
    }

  if ( streamInCnt == -1 && streamOutCnt == -1 )
    Errorc("I/O stream counts unlimited no allowed!");
    
  if ( streamInCnt == -1 )
    {
      streamInCnt = Process[processID].streamCnt - streamOutCnt;
      if ( streamInCnt < 1 ) Errorc("Input streams missing!");
    }

  if ( streamOutCnt == -1 )
    {
      streamOutCnt = Process[processID].streamCnt - streamInCnt;
      if ( streamInCnt < 1 ) Errorc("Output streams missing!");
    }

  streamCnt = streamInCnt + streamOutCnt;

  if ( Process[processID].streamCnt > streamCnt )
    Errorc("Too many streams!"
	   " Operator needs %d input and %d output streams.", streamInCnt, streamOutCnt);

  if ( Process[processID].streamCnt < streamCnt )
    Errorc("Too few streams specified!"
	   " Operator needs %d input and %d output streams.", streamInCnt, streamOutCnt);


  for ( i = streamInCnt; i < streamCnt; i++ )
    {
      if ( Process[processID].streamNames[i][0] == '-' )
	{
	  Errorc("Output file name %s must not begin with \"-\"!\n",
		 Process[processID].streamNames[i]);
	}
      else if ( !obase )
	{
	  for ( j = 0; j < streamInCnt; j++ ) /* does not work with files in pipes */
	    if ( strcmp(Process[processID].streamNames[i], Process[processID].streamNames[j]) == 0 )
	      Errorc("Output file name %s is equal to input file name"
		     " on position %d!\n", Process[processID].streamNames[i], j+1);
	}
    }  
}

#define  MAX_ARGV  8192

static
void setStreams(const char *argument)
{
  int processID = processSelf();
  int streamCnt;
  int i;
  int argc = 0;
  char *argv[MAX_ARGV];
  char *string;
  size_t arglen;

  arglen = 1 + strlen(argument);
  string = (char *) malloc(arglen);
  strcpy(string, argument);

  argv[argc++] = string;
  for ( i = 1; i < (int) arglen-1; i++ )
    {
      if ( string[i] == ' ' )
	{
	  string[i] = '\0';
	  argv[argc++] = &string[i+1];
	  if ( argc >= MAX_ARGV )
	    Error("Internal problem! More than %d arguments.", argc);
	}
    }

  streamCnt = getStreamCnt(argc, argv);

  Process[processID].nvals = 0;
  Process[processID].nvars = 0;
  Process[processID].ntimesteps = 0;

  Process[processID].streamCnt  = 0; /* filled in setStreamNames */
  if ( streamCnt )
    Process[processID].streamNames = (char **) malloc(streamCnt*sizeof(char *));
  for ( i = 0; i < streamCnt; i++ ) Process[processID].streamNames[i] = NULL;

  setStreamNames(argc, argv);

  checkStreamCnt();

  if ( Process[processID].streamCnt != streamCnt )
    Error("Internal problem with stream count %d %d", Process[processID].streamCnt, streamCnt);
  /*
  for ( i = 0; i < streamCnt; i++ )
    fprintf(stderr, "stream %d %s\n", i+1, Process[processID].streamNames[i]);
  */

  free(argv[0]);
}


void processDefArgument(const char *argument)
{
  int processID = processSelf();
  char *operatorArg;
  char *commapos;
  int oargc = 0;
  char **oargv = Process[processID].oargv;

  /*printf("argument: %s\n", argument);*/
  Process[processID].xoperator    = getOperator(argument);
  Process[processID].operatorName = getOperatorName(Process[processID].xoperator);
  Process[processID].operatorArg  = getOperatorArg(Process[processID].xoperator);
  operatorArg = Process[processID].operatorArg;

  if ( operatorArg )
    {
      oargv[oargc++] = operatorArg;
      /*printf("%d %s\n", oargc, operatorArg);*/

      commapos = operatorArg;
      while ( (commapos = strchr(commapos, ',')) != NULL )
	{
	  *commapos++ = '\0';
	  if ( strlen(commapos) )
	    {
	      if ( oargc >= MAX_ARGC )
		cdoAbort("Too many parameter (limit=%d)!", MAX_ARGC);

	      oargv[oargc++] = commapos;
	    }
	}
      Process[processID].oargc = oargc;
    }

  processDefPrompt(Process[processID].operatorName);

  setStreams(argument);
}

void processDefVarNum(int nvars, int streamID)
{
  int processID = processSelf();

  /*  if ( streamID == Process[processID].streams[0] ) */
    Process[processID].nvars += nvars;
}


int processInqVarNum(void)
{
  int processID = processSelf();

  return (Process[processID].nvars);
}


void processDefTimesteps(int streamID)
{
  int processID = processSelf();
  /*
  int i;
  printf("streamID %d %d %d %d\n", streamID, Process[processID].streams[0], Process[processID].streams[1], processID);

  for ( i = 0; i < Process[processID].nstream; i++)
    printf("streamID %d %d %d %d << \n", processID, Process[processID].nstream, i, Process[processID].streams[i]);
  */
  /*  if ( streamID == Process[processID].streams[0] )*/
    Process[processID].ntimesteps++;
}


int processInqTimesteps(void)
{
  int processID = processSelf();

  return (Process[processID].ntimesteps);
}


void processDelete(void)
{
  int processID = processSelf();

  //fprintf(stderr, "delete processID %d\n", processID);
#if  defined  (HAVE_LIBPTHREAD)
  pthread_mutex_lock(&processMutex);
#endif
  NumProcessActive--;
#if  defined  (HAVE_LIBPTHREAD)
  pthread_mutex_unlock(&processMutex);  
#endif
}


int operatorArgc(void)
{
  int processID = processSelf();

  return (Process[processID].oargc);
}


char **operatorArgv(void)
{
  int processID = processSelf();

  return (Process[processID].oargv);
}


void operatorCheckArgc(int numargs)
{
  int processID = processSelf();
  int argc = Process[processID].oargc;

  if ( argc < numargs )
    cdoAbort("Too few arguments! Need %d found %d.", numargs, argc);
  else if ( argc > numargs )
    cdoAbort("Too many arguments! Need %d found %d.", numargs, argc);
}


void operatorInputArg(const char *enter)
{
  char line[1024];
  char *pline = line;
  int processID = processSelf();
  size_t pos, len, linelen;
  int oargc;
  int lreadline;

  oargc = Process[processID].oargc;

  if ( oargc ) return;

  while ( oargc == 0 )
    {
      lreadline = 1;

      if ( enter ) fprintf(stderr, "%-16s : Enter %s > ", processInqPrompt(), enter);

      while ( lreadline )
	{
	  readline(stdin, pline, 1024);

	  lreadline = 0;
	  while ( 1 )
	    {
	      pos = 0;
	      while ( pline[pos] == ' ' || pline[pos] == ',' ) pos++;
	      pline += pos;
	      linelen = strlen(pline);
	      if ( linelen > 0 )
		{
		  if ( pline[0] == '\\' )
		    {
		      lreadline = 1;
		      break;
		    }
		  len = 0;
		  while ( pline[len] != ' '  && pline[len] != ',' &&
			  pline[len] != '\\' && len < linelen ) len++;

		  Process[processID].oargv[oargc] = (char *) malloc(len+1);
		  memcpy(Process[processID].oargv[oargc], pline, len);
		  Process[processID].oargv[oargc][len] = '\0';
		  oargc++;

		  pline += len;
		}
	      else
		break;
	    }
	}
    }

  Process[processID].oargc = oargc;
}


int cdoOperatorAdd(const char *name, int f1, int f2, const char *enter)
{
  int processID = processSelf();
  int operID = Process[processID].noper;

  if ( operID < MAX_OPERATOR )
    {
      Process[processID].operator[operID].f1     = f1;
      Process[processID].operator[operID].f2     = f2;
      Process[processID].operator[operID].name   = name;
      Process[processID].operator[operID].enter  = enter;

      Process[processID].noper++;
    }
  else
    {
      cdoAbort("Maximum of %d operators reached!", MAX_OPERATOR);
    }

  return (operID);
}


int cdoOperatorID(void)
{
  int processID = processSelf();
  int operID = -1;

  if ( Process[processID].noper > 0 )
    {
      for ( operID = 0; operID < Process[processID].noper; operID++ )
	if ( Process[processID].operator[operID].name )
	  if ( strcmp(Process[processID].operatorName, Process[processID].operator[operID].name) == 0 ) break;

      if ( operID == Process[processID].noper )
	cdoAbort("Operator not callable by this name!");
    }
  else
    {
      cdoAbort("Operator not initialized!");
    }

  return (operID);
}


int cdoOperatorF1(int operID)
{
  int processID = processSelf();

  return (Process[processID].operator[operID].f1);
}


int cdoOperatorF2(int operID)
{
  int processID = processSelf();

  return (Process[processID].operator[operID].f2);
}


const char *cdoOperatorName(int operID)
{
  int processID = processSelf();

  return (Process[processID].operator[operID].name);
}


const char *cdoOperatorEnter(int operID)
{
  int processID = processSelf();

  return (Process[processID].operator[operID].enter);
}


int cdoStreamNumber()
{
  int processID = processSelf();

  return (operatorStreamNumber(Process[processID].operatorName));
}
