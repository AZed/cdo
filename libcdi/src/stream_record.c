#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>

#include "dmemory.h"

#include "cdi.h"
#include "cdf_int.h"
#include "cdi_int.h"
#include "stream_grb.h"
#include "stream_cdf.h"
#include "stream_srv.h"
#include "stream_ext.h"
#include "stream_ieg.h"


void recordInitEntry(record_t *record)
{
  (*record).position = CDI_UNDEFID;
  (*record).size     = 0;
  (*record).param    = 0;
  (*record).ilevel   = CDI_UNDEFID;
  (*record).used     = FALSE;
  (*record).varID    = CDI_UNDEFID;
  (*record).levelID  = CDI_UNDEFID;
  memset((*record).varname, 0, sizeof((*record).varname));
}


int recordNewEntry(stream_t *streamptr, int tsID)
{
  int recordID = 0;
  int recordSize;
  record_t *records;

  recordSize = streamptr->tsteps[tsID].recordSize;
  records    = streamptr->tsteps[tsID].records;
  /*
    Look for a free slot in record.
    (Create the table the first time through).
  */
  if ( ! recordSize )
    {
      int i;
      recordSize = 1;   /*  <<<<----  */
      records = (record_t *) malloc(recordSize*sizeof(record_t));
      if ( records == NULL )
	{
          Message("recordSize = %d", recordSize);
	  SysError("Allocation of record_tTABLE failed");
	}

      for ( i = 0; i < recordSize; i++ )
	records[i].used = CDI_UNDEFID;
    }
  else
    {
      while ( recordID < recordSize )
	{
	  if ( records[recordID].used == CDI_UNDEFID ) break;
	  recordID++;
	}
    }
  /*
    If the table overflows, double its size.
  */
  if ( recordID == recordSize )
    {
      int i;

      recordSize = 2*recordSize;
      records    = (record_t *) realloc(records, recordSize*sizeof(record_t));
      if ( records == NULL )
	{
          Message("recordSize = %d", recordSize);
	  SysError("Reallocation of record_tTABLE failed");
	}
      recordID = recordSize/2;

      for ( i = recordID; i < recordSize; i++ )
	records[i].used = CDI_UNDEFID;
    }


  recordInitEntry(&records[recordID]);

  records[recordID].used = 1;

  streamptr->tsteps[tsID].recordSize = recordSize;
  streamptr->tsteps[tsID].records    = records;

  return (recordID);
}

static
void cdiInitRecord(stream_t *streamptr)
{
  streamptr->record = (Record *) malloc(sizeof(Record));

  streamptr->record->used       = 0;
  streamptr->record->nrec       = 0;
  streamptr->record->dataread   = 1;
  streamptr->record->param      = 0;
  streamptr->record->level      = 0;
  streamptr->record->date       = 0;
  streamptr->record->time       = 0;
  streamptr->record->gridID     = 0;
  streamptr->record->zaxisID    = 0;
  streamptr->record->buffer     = NULL;
  streamptr->record->buffersize = 0;
  streamptr->record->position   = 0;
  streamptr->record->varID      = 0;
  streamptr->record->levelID    = CDI_UNDEFID;
  streamptr->record->recid      = 0;
}


void streamInqRecord(int streamID, int *varID, int *levelID)
{
  /* int rec = 0; */
  int recID, tsID, rindex;
  int lindex;
  stream_t *streamptr;

  check_parg(varID);
  check_parg(levelID);

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  cdiDefAccesstype(streamID, TYPE_REC);

  if ( ! streamptr->record ) cdiInitRecord(streamptr);

  tsID   = streamptr->curTsID;
  rindex = streamptr->tsteps[tsID].curRecID + 1;

  if ( rindex >= streamptr->tsteps[tsID].nrecs )
    Error("record %d not available at timestep %d", rindex+1, tsID+1);

  recID  = streamptr->tsteps[tsID].recIDs[rindex];

  if ( recID == -1 || recID >= streamptr->tsteps[tsID].nallrecs )
    Error("Internal problem! tsID = %d recID = %d", tsID, recID);

  *varID   = streamptr->tsteps[tsID].records[recID].varID;
  lindex   = streamptr->tsteps[tsID].records[recID].levelID;

  *levelID = streamptr->vars[*varID].lindex[lindex];

  if ( CDI_Debug )
    Message("tsID = %d, recID = %d, varID = %d, levelID = %d\n",
	    tsID, recID, *varID, *levelID);

  streamptr->curTsID = tsID;
  streamptr->tsteps[tsID].curRecID = rindex;

  /*
  rec = recID + 1;
  filetype = streamptr->filetype;

  switch ( filetype )
    {
    case FILETYPE_GRB:
      {
        rec = grbInqRecord(streamID, varID, levelID);
	break;
      }
    case FILETYPE_SRV:
      {
        rec = srvInqRecord(streamptr, varID, levelID);
	break;
      }
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
	rec = cdfInqRecord(streamID, varID, levelID);
	break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }
  */
}


void streamDefRecord(int streamID, int varID, int levelID)
{
  int status = 0;
  int filetype;
  int param, gridID, zaxisID, level;
  int tsID;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  tsID = streamptr->curTsID;

  if ( tsID == CDI_UNDEFID )
    {
      tsID++;
      streamDefTimestep(streamID, tsID);
    }

  if ( ! streamptr->record ) cdiInitRecord(streamptr);

  vlistID = streamptr->vlistID;
  gridID  = vlistInqVarGrid(vlistID, varID);
  zaxisID = vlistInqVarZaxis(vlistID, varID);
  param   = vlistInqVarParam(vlistID, varID);
  level   = (int) zaxisInqLevel(zaxisID, levelID);

  streamptr->record->varID    = varID;
  streamptr->record->levelID  = levelID;
  streamptr->record->param    = param;
  streamptr->record->level    = level;
  streamptr->record->date     = streamptr->tsteps[tsID].taxis.vdate;
  streamptr->record->time     = streamptr->tsteps[tsID].taxis.vtime;
  streamptr->record->gridID   = gridID;
  streamptr->record->zaxisID  = zaxisID;
  streamptr->record->prec     = vlistInqVarDatatype(vlistID, varID);

  filetype = streamptr->filetype;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        /* FIXME: return value not inspected */
        status = grbDefRecord(streamptr);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        /* FIXME: return value not inspected */
        status = srvDefRecord(streamptr);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        /* FIXME: return value not inspected */
        status = extDefRecord(streamptr);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        /* FIXME: return value not inspected */
        status = iegDefRecord(streamptr);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
	if ( streamptr->accessmode == 0 ) cdfEndDef(streamptr);
	status = cdfDefRecord(streamptr);
	break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }
}


void streamReadRecord(int streamID, double *data, int *nmiss)
{
  int status = 0;
  int filetype;
  stream_t *streamptr;

  check_parg(data);
  check_parg(nmiss);

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  filetype = streamptr->filetype;

  *nmiss = 0;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        /* FIXME: return value not inspected */
        status = grbReadRecord(streamptr, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        /* FIXME: return value not inspected */
        status = srvReadRecord(streamptr, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        /* FIXME: return value not inspected */
        status = extReadRecord(streamptr, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        /* FIXME: return value not inspected */
        status = iegReadRecord(streamptr, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
        /* FIXME: return value not inspected */
	status = cdfReadRecord(streamptr, data, nmiss);
	break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }
}


void stream_write_record(int streamID, int memtype, const void *data, int nmiss)
{
  int status = 0;
  int filetype;
  stream_t *streamptr;

  check_parg(data);

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  filetype = streamptr->filetype;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        status = grb_write_record(streamptr, memtype, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("srvWriteRecord not implemented for memtype float!");
        status = srvWriteRecord(streamptr, (const double*) data);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("extWriteRecord not implemented for memtype float!");
        status = extWriteRecord(streamptr, (const double*) data);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("iegWriteRecord not implemented for memtype float!");
        status = iegWriteRecord(streamptr, (const double*) data);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
	cdf_write_record(streamptr, memtype, data, nmiss);
	break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }
}


void streamWriteRecord(int streamID, const double *data, int nmiss)
{
  stream_write_record(streamID, MEMTYPE_DOUBLE, (const void *) data, nmiss);
}

void streamWriteRecordF(int streamID, const float *data, int nmiss)
{
  stream_write_record(streamID, MEMTYPE_FLOAT, (const void *) data, nmiss);
}


void streamCopyRecord(int streamID2, int streamID1)
{
  int status = 0;
  int filetype = CDI_UNDEFID, filetype1, filetype2;
  stream_t *streamptr1;
  stream_t *streamptr2;

  streamptr1 = stream_to_pointer(streamID1);
  streamptr2 = stream_to_pointer(streamID2);

  stream_check_ptr(__func__, streamptr1);
  stream_check_ptr(__func__, streamptr2);

  filetype1 = streamptr1->filetype;
  filetype2 = streamptr2->filetype;

  if ( filetype1 == filetype2 ) filetype = filetype2;
  else
    {
      switch (filetype1)
        {
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
          switch (filetype2)
            {
            case FILETYPE_NC:
            case FILETYPE_NC2:
            case FILETYPE_NC4:
            case FILETYPE_NC4C:
              Warning("Streams have different file types (%s -> %s)!", strfiletype(filetype1), strfiletype(filetype2));
              filetype = filetype2;
              break;
            }
          break;
        }
    }

  if ( filetype == CDI_UNDEFID )
    Error("Streams have different file types (%s -> %s)!", strfiletype(filetype1), strfiletype(filetype2));

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        /* FIXME: return value not inspected */
	status = grbCopyRecord(streamptr2, streamptr1);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        /* FIXME: return value not inspected */
	status = srvCopyRecord(streamptr2, streamptr1);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        /* FIXME: return value not inspected */
	status = extCopyRecord(streamptr2, streamptr1);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        /* FIXME: return value not inspected */
	status = iegCopyRecord(streamptr2, streamptr1);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
        /* FIXME: return value not inspected */
	status = cdfCopyRecord(streamptr2, streamptr1);
	break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }
}


void cdi_create_records(stream_t *streamptr, int tsID)
{
  int nrecords, maxrecords;
  int nvars, varID, recID;
  record_t *records;
  int vlistID;

  if ( streamptr->tsteps[tsID].records ) return;

  vlistID  = streamptr->vlistID;

  if ( tsID == 0 )
    {
      maxrecords = 0;
      nvars = streamptr->nvars;
      for ( varID = 0; varID < nvars; varID++)
	maxrecords += streamptr->vars[varID].nlevs;
    }
  else
    maxrecords = streamptr->tsteps[0].recordSize;

  if ( tsID == 0 )
    {
      nrecords = maxrecords;
    }
  else if ( tsID == 1 )
    {
      nrecords = 0;
      maxrecords = streamptr->tsteps[0].recordSize;
      for ( recID = 0; recID < maxrecords; recID++ )
	{
	  varID = streamptr->tsteps[0].records[recID].varID;
	  if ( varID != -1 ) /* varID = -1 for write mode !!! */
	    if ( vlistInqVarTsteptype(vlistID, varID) == TSTEP_CONSTANT )
	      continue;
	  nrecords++;
	}
    }
  else
    nrecords = streamptr->tsteps[1].nallrecs;

  if ( maxrecords > 0 )
    records = (record_t *) malloc(maxrecords*sizeof(record_t));
  else
    records = NULL;

  streamptr->tsteps[tsID].records    = records;
  streamptr->tsteps[tsID].recordSize = maxrecords;
  streamptr->tsteps[tsID].nallrecs   = nrecords;

  if ( tsID == 0 )
    {
      for ( recID = 0; recID < maxrecords; recID++ )
	recordInitEntry(&streamptr->tsteps[tsID].records[recID]);
    }
  else
    {
      memcpy(streamptr->tsteps[tsID].records,
	     streamptr->tsteps[0].records,
	     maxrecords*sizeof(record_t));

      for ( recID = 0; recID < maxrecords; recID++ )
	{
	  varID = streamptr->tsteps[0].records[recID].varID;
	  if ( varID != -1 ) /* varID = -1 for write mode !!! */
	    if ( vlistInqVarTsteptype(vlistID, varID) != TSTEP_CONSTANT )
	      {
		streamptr->tsteps[tsID].records[recID].position = CDI_UNDEFID;
		streamptr->tsteps[tsID].records[recID].size     = 0;
		streamptr->tsteps[tsID].records[recID].used     = FALSE;
	      }
	}
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
