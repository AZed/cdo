#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "dmemory.h"

#include "cdi.h"
#include "cdi_int.h"


static
void tstepsInitEntry(stream_t *streamptr, int tsID)
{
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


int tstepsNewEntry(stream_t *streamptr)
{
  int tsID = 0;
  int tstepsTableSize;
  tsteps_t *tstepsTable;

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
      tstepsTable = (tsteps_t *) realloc(tstepsTable, tstepsTableSize*sizeof(tsteps_t));
      if ( tstepsTable == NULL )
	{
          Message("tstepsTableSize = %d", tstepsTableSize);
	  SysError("Reallocation of tsteps_t failed");
	}
    }

  streamptr->tstepsTableSize = tstepsTableSize;
  streamptr->tsteps          = tstepsTable;

  tstepsInitEntry(streamptr, tsID);

  streamptr->tsteps[tsID].taxis.used = TRUE;

  return (tsID);
}


void cdiCreateTimesteps(stream_t *streamptr)
{
  int ntsteps;
  int tsID;

  if ( streamptr->ntsteps < 0 || streamptr->tstepsTableSize > 0 )
    return;

  if ( streamptr->ntsteps == 0 ) ntsteps = 1;    /* <<<<<-------- */
  else ntsteps = streamptr->ntsteps;

  streamptr->tsteps = (tsteps_t *) malloc(ntsteps*sizeof(tsteps_t));
  if ( streamptr->tsteps == NULL )
    SysError("Allocation of tsteps_t failed");

  streamptr->tstepsTableSize = ntsteps;
  streamptr->tstepsNextID    = ntsteps;

  for ( tsID = 0; tsID < ntsteps; tsID++ )
    {
      tstepsInitEntry(streamptr, tsID);
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
