#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>


int _ExitOnError   = 1;	/* If set to 1, exit on error       */
int _Verbose = 1;	/* If set to 1, errors are reported */
int _Debug   = 0;       /* If set to 1, debugging           */


void SysError_(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  printf("\n");
   fprintf(stderr, "Error (%s) : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  if ( errno )
    perror("System error message ");
	
  exit(EXIT_FAILURE);
}


void Error_(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  printf("\n");
   fprintf(stderr, "Error (%s) : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  if ( _ExitOnError ) exit(EXIT_FAILURE);
}


void Warning_(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  if ( _Verbose )
    {
       fprintf(stderr, "Warning (%s) : ", caller);
      vfprintf(stderr, fmt, args);
       fprintf(stderr, "\n");
    }

  va_end(args);
}


void Message_(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

   fprintf(stdout, "%-18s : ", caller);
  vfprintf(stdout, fmt, args);
   fprintf(stdout, "\n");

  va_end(args);
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
