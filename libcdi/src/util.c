#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <float.h>
#include <sys/types.h>

#include "cdi.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "binary.h"


#undef  IsBigendian
#define IsBigendian()  ( u_byteorder.c[sizeof(long) - 1] )

void cdiPrintDatatypes(void)
{
  /* IsBigendian returns 1 for big endian byte order */
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};

  fprintf (stderr, "+-------------+-------+\n");
  fprintf (stderr, "| types       | bytes |\n");
  fprintf (stderr, "+-------------+-------+\n");
  fprintf (stderr, "| void *      |   %3d |\n", (int) sizeof(void *));
  fprintf (stderr, "+-------------+-------+\n");
  fprintf (stderr, "| char        |   %3d |\n", (int) sizeof(char));
  fprintf (stderr, "+-------------+-------+\n");
  fprintf (stderr, "| short       |   %3d |\n", (int) sizeof(short));
  fprintf (stderr, "| int         |   %3d |\n", (int) sizeof(int));
  fprintf (stderr, "| long        |   %3d |\n", (int) sizeof(long));
  fprintf (stderr, "| long long   |   %3d |\n", (int) sizeof(long long));
  fprintf (stderr, "| size_t      |   %3d |\n", (int) sizeof(size_t));
  fprintf (stderr, "| off_t       |   %3d |\n", (int) sizeof(off_t));
  fprintf (stderr, "+-------------+-------+\n");
  fprintf (stderr, "| float       |   %3d |\n", (int) sizeof(float));
  fprintf (stderr, "| double      |   %3d |\n", (int) sizeof(double));
  fprintf (stderr, "| long double |   %3d |\n", (int) sizeof(long double));
  fprintf (stderr, "+-------------+-------+\n\n");
#define XSTRING(x)	#x
#define STRING(x)	XSTRING(x)
  fprintf (stderr, "+-------------+-----------+\n");
  fprintf (stderr, "| INT32       | %-9s |\n", STRING(INT32));
  fprintf (stderr, "| INT64       | %-9s |\n", STRING(INT64));
  fprintf (stderr, "| FLT32       | %-9s |\n", STRING(FLT32));
  fprintf (stderr, "| FLT64       | %-9s |\n", STRING(FLT64));
  fprintf (stderr, "+-------------+-----------+\n");

  if ( IsBigendian() )
    fprintf (stderr, "\n  byte ordering is BIGENDIAN\n\n");
  else
    fprintf (stderr, "\n  byte ordering is LITTLEENDIAN\n\n");
}

static char uuidFmt[] = "%02hhx%02hhx%02hhx%02hhx-"
  "%02hhx%02hhx-%02hhx%02hhx-%02hhx%02hhx-"
  "%02hhx%02hhx%02hhx%02hhx%02hhx%02hhx";

enum {
  uuidNumHexChars = 36,
};

void uuid2str(const unsigned char *uuid, char *uuidstr)
{

  if ( uuid == NULL || uuidstr == NULL ) return;

  int iret = sprintf(uuidstr, uuidFmt,
                     uuid[0], uuid[1], uuid[2], uuid[3],
                     uuid[4], uuid[5], uuid[6], uuid[7],
                     uuid[8], uuid[9], uuid[10], uuid[11],
                     uuid[12], uuid[13], uuid[14], uuid[15]);

  if ( iret != uuidNumHexChars ) uuidstr[0] = 0;
}


int str2uuid(const char *uuidstr, unsigned char *uuid)
{
  if ( uuid == NULL || uuidstr == NULL || strlen(uuidstr) != uuidNumHexChars)
    return -1;

  int iret = sscanf(uuidstr, uuidFmt,
                    &uuid[0], &uuid[1], &uuid[2], &uuid[3],
                    &uuid[4], &uuid[5], &uuid[6], &uuid[7],
                    &uuid[8], &uuid[9], &uuid[10], &uuid[11],
                    &uuid[12], &uuid[13], &uuid[14], &uuid[15]);
  if ( iret != CDI_UUID_SIZE ) return -1;
  return iret;
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
