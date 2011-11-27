#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>

#include "error.h"
#include "binary.h"

void swap4byte(void *ptr, size_t size)
{
  INT32 *ptrtmp;
  int nval;

  nval = size;
  if ( nval < 0 ) nval = 0;
  ptrtmp = (INT32 *) ptr;

  if ( sizeof(INT32) == 4 )
    {
      while ( nval-- )
	{
	  *ptrtmp = (((*ptrtmp >> 24) & 0x00ff) | ((*ptrtmp & 0x00ff) << 24) |
		     ((*ptrtmp >>  8) & 0xff00) | ((*ptrtmp & 0xff00) <<  8));
	  ptrtmp++;
	}
    }
  else
    {
      Error("not implemented for %d byte data", sizeof(INT32));
    }
}

void swap8byte(void *ptr, size_t size)
{
  INT64 *ptrtmp;
  int nval;

  nval = size;
  if ( nval < 0 ) nval = 0;
  ptrtmp = (INT64 *) ptr;

  if ( sizeof(INT64) == 8 )
    {
      while ( nval-- )
	{
	  *ptrtmp = (((*ptrtmp >> 56) & 0x000000ff) | ((*ptrtmp & 0x000000ff) << 56) |
		     ((*ptrtmp >> 40) & 0x0000ff00) | ((*ptrtmp & 0x0000ff00) << 40) |
		     ((*ptrtmp >> 24) & 0x00ff0000) | ((*ptrtmp & 0x00ff0000) << 24) |
		     ((*ptrtmp >>  8) & 0xff000000) | ((*ptrtmp & 0xff000000) <<  8));
	  ptrtmp++;
	}
    }
  else
    {
      Error("not implemented for %d byte data", sizeof(INT64));
    }
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
