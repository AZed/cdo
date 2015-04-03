#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <float.h>
#include <sys/types.h>

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


void uuid2str(const char *uuid, char *uuidstr)
{
  int iret;
  unsigned int ui[16];

  if ( uuid == NULL ) return;
  if ( uuidstr == NULL ) return;

  uuidstr[0] = 0;

  for ( int i = 0; i < 16; ++i ) ui[i] = (unsigned char) uuid[i];

  iret = sprintf(uuidstr, "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x",
		 ui[0], ui[1], ui[2], ui[3], ui[4], ui[5], ui[6], ui[7],
		 ui[8], ui[9], ui[10], ui[11], ui[12], ui[13], ui[14], ui[15]);

  if ( iret != 36 ) uuidstr[0] = 0;
}


void str2uuid(const char *uuidstr, char *uuid)
{
  int iret;
  unsigned int ui[16];

  if ( uuid == NULL ) return;
  if ( uuidstr == NULL ) return;

  uuid[0] = 0;

  if ( strlen(uuidstr) != 36 ) return;

  iret = sscanf(uuidstr, "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x",
		&ui[0], &ui[1], &ui[2], &ui[3], &ui[4], &ui[5], &ui[6], &ui[7],
		&ui[8], &ui[9], &ui[10], &ui[11], &ui[12], &ui[13], &ui[14], &ui[15]);

  if ( iret != 16 ) return;

  for ( int i = 0; i < 16; ++i ) uuid[i] = ui[i];

  uuid[16] = 0;
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
