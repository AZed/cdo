#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#if  defined  (HAVE_LIBGRIB_API)
#  include <grib_api.h>
#endif

#include <stdio.h>

#include "cdi.h"
#include "stream_int.h"
#include "gribapi.h"
#include "dmemory.h"

#define XSTRING(x)	#x
#define STRING(x)	XSTRING(x)

static char gribapi_libvers[64] = "";

const char *gribapiLibraryVersion(void)
{
#if  defined  (HAVE_LIBGRIB_API)
  long version = grib_get_api_version();
  int major_version, minor_version, revision_version;

  major_version    = version/10000;
  minor_version    = (version-major_version*10000)/100;
  revision_version = (version-major_version*10000-minor_version*100);

  sprintf(gribapi_libvers, "%d.%d.%d",
	  major_version, minor_version, revision_version);
#endif

  return (gribapi_libvers);
}


void gribContainersNew(int streamID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

#if  defined  (HAVE_LIBCGRIBEX)
  if ( streamptr->filetype == FILETYPE_GRB )
    {
    }
  else
#endif
    {
      int i, editionNumber = 2;
      gribContainer_t *gribContainers;

      if ( streamptr->filetype == FILETYPE_GRB ) editionNumber = 1;

      gribContainers = (gribContainer_t *) malloc(streamptr->nvars*sizeof(gribContainer_t));
      streamptr->gribContainers = (void *) gribContainers;

      for ( i = 0; i < streamptr->nvars; ++i )
	{
	  gribContainers[i].gribHandle = gribHandleNew(editionNumber);
	  gribContainers[i].init = FALSE;
	}
    }
}


void gribContainersDelete(int streamID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( streamptr->gribContainers )
    {
      int i;
      gribContainer_t *gribContainers = (gribContainer_t *) streamptr->gribContainers;

      for ( i = 0; i < streamptr->nvars; ++i )
	{
	  gribHandleDelete(gribContainers[i].gribHandle);
	}
      free(streamptr->gribContainers);
      streamptr->gribContainers = NULL;
    }
}
