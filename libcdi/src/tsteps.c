#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "dmemory.h"

#include "cdi.h"
#include "stream_int.h"


static void tstepsInitEntry(int streamID, int tsID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  streamptr->tsteps[tsID].curRecID     = CDI_UNDEFID;
  streamptr->tsteps[tsID].position     = 0;
  streamptr->tsteps[tsID].records      = NULL;
  streamptr->tsteps[tsID].recordSize   = 0;
  streamptr->tsteps[tsID].nallrecs     = 0;
  streamptr->tsteps[tsID].recIDs       = NULL;
  streamptr->tsteps[tsID].nrecs        = 0;
  streamptr->tsteps[tsID].next         = 0;

  ptaxisInit(&streamptr->tsteps[tsID].taxis);
}

int tstepsNewEntry(int streamID)
{
  int tsID = 0;
  int tstepsTableSize;
  TSTEPS *tstepsTable;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  tsID            = streamptr->tstepsNextID++;
  tstepsTableSize = streamptr->tstepsTableSize;
  tstepsTable     = streamptr->tsteps;

  /*
    If the table overflows, double its size.
  */
  if ( tsID == tstepsTableSize )
    {
      if ( tstepsTableSize == 0 ) tstepsTableSize = 1;
      tstepsTableSize = 2*tstepsTableSize;
      tstepsTable = (TSTEPS *) realloc(tstepsTable, tstepsTableSize*sizeof(TSTEPS));
      if ( tstepsTable == NULL )
	{
          Message("tstepsTableSize = %d", tstepsTableSize);
	  SysError("Reallocation of TSTEPS failed");
	}
    }

  streamptr->tstepsTableSize = tstepsTableSize;
  streamptr->tsteps          = tstepsTable;

  tstepsInitEntry(streamID, tsID);

  streamptr->tsteps[tsID].taxis.used = TRUE;

  return (tsID);
}

void cdiCreateTimesteps(int streamID)
{
  int ntsteps;
  int tsID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( streamptr->ntsteps < 0 || streamptr->tstepsTableSize > 0 )
    return;

  if ( streamptr->ntsteps == 0 ) ntsteps = 1;    /* <<<<<-------- */
  else ntsteps = streamptr->ntsteps;

  streamptr->tsteps = (TSTEPS *) malloc(ntsteps*sizeof(TSTEPS));
  if ( streamptr->tsteps == NULL )
    SysError("Allocation of TSTEPS failed");

  streamptr->tstepsTableSize = ntsteps;
  streamptr->tstepsNextID    = ntsteps;

  for ( tsID = 0; tsID < ntsteps; tsID++ )
    {
      tstepsInitEntry(streamID, tsID);
      streamptr->tsteps[tsID].taxis.used = TRUE;
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
