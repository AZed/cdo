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
  set_text_color(stderr, RESET, RED);
   fprintf(stderr, "%s: ", processInqPrompt());
  reset_text_color(stderr);
  set_text_color(stderr, RESET, BLACK);
  vfprintf(stderr, fmt, args);
  reset_text_color(stderr);
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
  set_text_color(stderr, RESET, RED);
   fprintf(stderr, "%s (Abort): ", processInqPrompt());
  reset_text_color(stderr);
  set_text_color(stderr, RESET, BLACK);
  vfprintf(stderr, fmt, args);
  reset_text_color(stderr);
   fprintf(stderr, "\n");

  va_end(args);

  if ( _ExitOnError ) exit(EXIT_FAILURE);
}


void cdoWarning(const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);

  set_text_color(stderr, BRIGHT, YELLOW);
   fprintf(stderr, "%s (Warning): ", processInqPrompt());
  reset_text_color(stderr);
  set_text_color(stderr, RESET, BLACK);
  vfprintf(stderr, fmt, args);
  reset_text_color(stderr);
   fprintf(stderr, "\n");

  va_end(args);
}


void cdoPrint(const char *fmt, ...)
{
  va_list args;

  if ( ! cdoSilentMode )
    {
      va_start(args, fmt);

      set_text_color(stderr, RESET, GREEN);
      fprintf(stderr, "%s: ", processInqPrompt());
      reset_text_color(stderr);
      set_text_color(stderr, RESET, BLACK);
      vfprintf(stderr, fmt, args);
      reset_text_color(stderr);
      fprintf(stderr, "\n");

      va_end(args);
    }
}
