#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <cdi.h>
#include "cdo.h"
#include "process.h"

static int _ExitOnError   = 1;	/* If set to 1, exit on error       */

void cdiError(int cdiErrno, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  printf("\n");
   fprintf(stderr, "%s: ", processInqPrompt());
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  fprintf(stderr, "%s\n", cdiStringError(cdiErrno));

  if ( _ExitOnError ) exit(EXIT_FAILURE);
}


void cdoAbort(const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  printf("\n");
   fprintf(stderr, "%s (Abort): ", processInqPrompt());
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  if ( _ExitOnError ) exit(EXIT_FAILURE);
}


void cdoWarning(const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);

   fprintf(stderr, "%s (Warning): ", processInqPrompt());
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);
}


void cdoPrint(const char *fmt, ...)
{
  va_list args;

  if ( ! cdoSilentMode )
    {
      va_start(args, fmt);

      fprintf(stderr, "%s: ", processInqPrompt());
      vfprintf(stderr, fmt, args);
      fprintf(stderr, "\n");

      va_end(args);
    }
}
