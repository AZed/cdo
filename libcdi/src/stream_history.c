#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "stream_cdf.h"


void streamDefHistory(int streamID, int length, const char *history)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( streamptr->filetype == FILETYPE_NC  ||
       streamptr->filetype == FILETYPE_NC2 ||
       streamptr->filetype == FILETYPE_NC4 ||
       streamptr->filetype == FILETYPE_NC4C )
    {
      char *histstring;
      size_t len;
      if ( history )
	{
	  len = strlen(history);
	  if ( len )
	    {
	      histstring = strdupx(history);
	      cdfDefHistory(streamptr, length, histstring);
	      free(histstring);
	    }
	}
    }
}


int streamInqHistorySize(int streamID)
{
  int size = 0;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( streamptr->filetype == FILETYPE_NC  ||
       streamptr->filetype == FILETYPE_NC2 ||
       streamptr->filetype == FILETYPE_NC4 ||
       streamptr->filetype == FILETYPE_NC4C )
    {
      size = cdfInqHistorySize(streamptr);
    }

  return (size);
}


void streamInqHistoryString(int streamID, char *history)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( streamptr->filetype == FILETYPE_NC  ||
       streamptr->filetype == FILETYPE_NC2 ||
       streamptr->filetype == FILETYPE_NC4 ||
       streamptr->filetype == FILETYPE_NC4C )
    {
      cdfInqHistoryString(streamptr, history);
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
