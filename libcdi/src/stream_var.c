#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <string.h>

#include "dmemory.h"
#include "error.h"

#include "cdi.h"
#include "cdi_int.h"


static
void streamvar_init_entry(stream_t *streamptr, int varID)
{
  streamptr->vars[varID].ncvarid      = CDI_UNDEFID;
  streamptr->vars[varID].defmiss      = 0;
  streamptr->vars[varID].nlevs        = 0;
  streamptr->vars[varID].level        = NULL;
  streamptr->vars[varID].lindex       = NULL;

  streamptr->vars[varID].gridID       = CDI_UNDEFID;
  streamptr->vars[varID].zaxisID      = CDI_UNDEFID;
  streamptr->vars[varID].tsteptype    = CDI_UNDEFID;
}

static
int streamvar_new_entry(stream_t *streamptr)
{
  int varID = 0;
  int streamvarSize;
  svarinfo_t *streamvar;

  streamvarSize = streamptr->varsAllocated;
  streamvar     = streamptr->vars;
  /*
    Look for a free slot in streamvar.
    (Create the table the first time through).
  */
  if ( ! streamvarSize )
    {
      int i;

      streamvarSize = 2;
      streamvar
        = (svarinfo_t *)xmalloc((size_t)streamvarSize * sizeof(svarinfo_t));
      if ( streamvar == NULL )
	{
          Message("streamvarSize = %d", streamvarSize);
	  SysError("Allocation of svarinfo_t failed");
	}

      for ( i = 0; i < streamvarSize; i++ )
	streamvar[i].isUsed = FALSE;
    }
  else
    {
      while ( varID < streamvarSize )
	{
	  if ( ! streamvar[varID].isUsed ) break;
	  varID++;
	}
    }
  /*
    If the table overflows, double its size.
  */
  if ( varID == streamvarSize )
    {
      int i;

      streamvarSize = 2*streamvarSize;
      streamvar
        = (svarinfo_t *)xrealloc(streamvar,
                                 (size_t)streamvarSize * sizeof (svarinfo_t));
      if ( streamvar == NULL )
	{
          Message("streamvarSize = %d", streamvarSize);
	  SysError("Reallocation of svarinfo_t failed");
	}
      varID = streamvarSize/2;

      for ( i = varID; i < streamvarSize; i++ )
	streamvar[i].isUsed = FALSE;
    }

  streamptr->varsAllocated = streamvarSize;
  streamptr->vars          = streamvar;

  streamvar_init_entry(streamptr, varID);

  streamptr->vars[varID].isUsed = TRUE;

  return (varID);
}


int stream_new_var(stream_t *streamptr, int gridID, int zaxisID)
{
  int varID;
  int *level;
  int *lindex;
  int nlevs;
  int levID;

  if ( CDI_Debug )
    Message("gridID = %d  zaxisID = %d", gridID, zaxisID);

  varID = streamvar_new_entry(streamptr);

  streamptr->nvars++;

  streamptr->vars[varID].gridID  = gridID;
  streamptr->vars[varID].zaxisID = zaxisID;

  nlevs = zaxisInqSize(zaxisID);

  level  = (int *)xmalloc((size_t)nlevs * sizeof (int));
  lindex = (int *)xmalloc((size_t)nlevs * sizeof (int));

  for ( levID = 0; levID < nlevs; levID++ )
    level[levID] = CDI_UNDEFID;

  for ( levID = 0; levID < nlevs; levID++ )
    lindex[levID] = levID;

  streamptr->vars[varID].nlevs  = nlevs;
  streamptr->vars[varID].level  = level;
  streamptr->vars[varID].lindex = lindex;

  return (varID);
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
