#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>

#include "dmemory.h"

#include "error.h"
#include "file.h"
#include "cdi.h"
#include "stream_int.h"
#include "varscan.h"
#include "datetime.h"
#include "service.h"
#include "vlist.h"


#undef  UNDEFID
#define UNDEFID  CDI_UNDEFID

#define SINGLE_PRECISION  4
#define DOUBLE_PRECISION  8

#if defined (HAVE_LIBSERVICE)


typedef struct {
  int param;
  int level;
} SRVCOMPVAR; 


int srvInqDatatype(int prec)
{
  int datatype;

  if ( prec == DOUBLE_PRECISION ) datatype = DATATYPE_FLT64;
  else                            datatype = DATATYPE_FLT32;

  return (datatype);
}


int srvDefDatatype(int datatype)
{
  int prec;

  if ( datatype == DATATYPE_CPX32 || datatype == DATATYPE_CPX64 )
    Error("CDI/SERVICE library does not support complex numbers!");

  if ( datatype != DATATYPE_FLT32 && datatype != DATATYPE_FLT64 )
    datatype = DATATYPE_FLT32;

  if ( datatype == DATATYPE_FLT64 ) prec = DOUBLE_PRECISION;
  else                              prec = SINGLE_PRECISION;

  return (prec);
}

/* not used
int srvInqRecord(int streamID, int *varID, int *levelID)
{
  int status;
  int fileID;
  int icode, ilevel;
  int zaxisID = -1;
  int header[8];
  int vlistID;
  srvrec_t *srvp;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);
  srvp    = streamptr->record->srvp;

  *varID   = -1;
  *levelID = -1;

  status = srvRead(fileID, srvp);
  if ( status != 0 ) return (0);

  srvInqHeader(srvp, header);

  icode  = header[0];
  ilevel = header[1];

  *varID = vlistInqVarID(vlistID, icode);

  if ( *varID == UNDEFID ) Error("Code %d undefined", icode);

  zaxisID = vlistInqVarZaxis(vlistID, *varID);

  *levelID = zaxisInqLevelID(zaxisID, (double) ilevel);
  
  return (1);
}
*/

int srvReadRecord(int streamID, double *data, int *nmiss)
{
  int vlistID, fileID;
  int status;
  int recID, vrecID, tsID;
  off_t recpos;
  int header[8];
  int varID, gridID;
  int i, size;
  double missval;
  srvrec_t *srvp;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);
  tsID    = streamptr->curTsID;
  vrecID  = streamptr->tsteps[tsID].curRecID;
  recID   = streamptr->tsteps[tsID].recIDs[vrecID];
  recpos  = streamptr->tsteps[tsID].records[recID].position;
  varID   = streamptr->tsteps[tsID].records[recID].varID;
  srvp    = streamptr->record->srvp;

  fileSetPos(fileID, recpos, SEEK_SET);

  status = srvRead(fileID, srvp);
  if ( status != 0 ) return (0);

  srvInqHeader(srvp, header);
  srvInqDataDP(srvp, data);

  missval = vlistInqVarMissval(vlistID, varID);
  gridID  = vlistInqVarGrid(vlistID, varID);
  size    = gridInqSize(gridID);

  streamptr->numvals += size;

  *nmiss = 0;
  for ( i = 0; i < size; i++ )
    if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
      {
	data[i] = missval;
	(*nmiss)++;
      }

  return (1);
}


int srvCopyRecord(int streamID2, int streamID1)
{
  int fileID1, fileID2;
  int tsID, recID, vrecID;
  long recsize;
  off_t recpos;
  int status = 0;
  char *buffer;
  stream_t *streamptr1;
  stream_t *streamptr2;

  streamptr1 = stream_to_pointer(streamID1);
  streamptr2 = stream_to_pointer(streamID2);

  stream_check_ptr(__func__, streamptr1);
  stream_check_ptr(__func__, streamptr2);

  fileID1 = streamInqFileID(streamID1);
  fileID2 = streamInqFileID(streamID2);

  tsID    = streamptr1->curTsID;
  vrecID  = streamptr1->tsteps[tsID].curRecID;
  recID   = streamptr1->tsteps[tsID].recIDs[vrecID];
  recpos  = streamptr1->tsteps[tsID].records[recID].position;
  recsize = streamptr1->tsteps[tsID].records[recID].size;

  fileSetPos(fileID1, recpos, SEEK_SET);

  buffer = (char *) malloc(recsize);

  fileRead(fileID1, buffer, recsize);

  fileWrite(fileID2, buffer, recsize);

  free(buffer);

  return (status);
}


int srvDefRecord(int streamID)
{
  int fileID;
  int gridID;
  int header[8];
  int status = 0;
  int xsize, ysize;
  int datatype;
  int pdis, pcat, pnum;
  srvrec_t *srvp;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  fileID = streamInqFileID(streamID);
  gridID = streamptr->record->gridID;
  srvp   = streamptr->record->srvp;

  cdiDecodeParam(streamptr->record->param, &pnum, &pcat, &pdis);
  header[0] = pnum;
  header[1] = streamptr->record->level;
  header[2] = streamptr->record->date;
  header[3] = streamptr->record->time;

  xsize = gridInqXsize(gridID);
  ysize = gridInqYsize(gridID);
  if ( xsize == 0 || ysize == 0 )
    {
      xsize = gridInqSize(gridID);
      ysize = 1;
    }
  if ( gridInqType(gridID) == GRID_UNSTRUCTURED ) ysize = 1;
  if ( gridInqSize(gridID) != xsize*ysize )
    Error("Internal problem with gridsize!");

  header[4] = xsize;
  header[5] = ysize;
  header[6] = 0;
  header[7] = 0;

  datatype = streamptr->record->prec;

  srvp->dprec = srvDefDatatype(datatype);

  srvDefHeader(srvp, header);

  return (status);
}


int srvWriteRecord(int streamID, const double *data)
{
  int fileID;
  int status = 0;
  srvrec_t *srvp;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  fileID = streamInqFileID(streamID);
  srvp   = streamptr->record->srvp;

  srvDefDataDP(srvp, data);

  srvWrite(fileID, srvp);

  return (status);
}


void srvAddRecord(int streamID, int param, int level, int xsize, int ysize,
		  long recsize, off_t position, int prec)
{
  int leveltype;
  int gridID = UNDEFID;
  int levelID = 0;
  int tsID, recID, varID;
  int datatype;
  record_t *record;
  grid_t grid;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID = streamInqVlist(streamID);
  tsID    = streamptr->curTsID;
  recID   = recordNewEntry(streamID, tsID);
  record  = &streamptr->tsteps[tsID].records[recID];

  (*record).size     = recsize;
  (*record).position = position;
  (*record).param    = param;
  (*record).ilevel   = level;

  memset(&grid, 0, sizeof(grid_t));
  grid.type  = GRID_GENERIC; 
  grid.size  = xsize*ysize;
  grid.xsize = xsize;
  grid.ysize = ysize;
  grid.xvals = NULL;
  grid.yvals = NULL;
  gridID = varDefGrid(vlistID, grid, 0);
  /*
  if ( level == 0 ) leveltype = ZAXIS_SURFACE;
  else              leveltype = ZAXIS_GENERIC;
  */
  leveltype = ZAXIS_GENERIC;

  datatype = srvInqDatatype(prec);

  varAddRecord(recID, param, gridID, leveltype, 0, level, 0,
	       datatype, &varID, &levelID, 0, 0, 0, NULL, NULL, NULL);

  (*record).varID   = varID;
  (*record).levelID = levelID;

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;

  if ( CDI_Debug )
    Message("varID = %d gridID = %d levelID = %d",
	    varID, gridID, levelID);
}


void srvCmpRecord(int streamID, int tsID, int recID, off_t position, int param,
		  int level, int xsize, int ysize)
{
  int varID = 0;
  int levelID = 0;
  record_t *record;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  record  = &streamptr->tsteps[tsID].records[recID];

  if ( param != (*record).param || level != (*record).ilevel )
    Error("inconsistent timestep");

  (*record).position = position;
  /*
  varID   = (*record).varID;
  levelID = (*record).levelID;

  streamptr->vars[varID].level[levelID] = recID;

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;
  */
  if ( CDI_Debug )
    Message("varID = %d levelID = %d", varID, levelID);
}

static
void srvScanTimestep1(int streamID)
{  
  int header[8];
  int prec = 0;
  int status;
  int fileID;
  int rxsize = 0, rysize = 0;
  int param = 0;
  int rcode = 0, rlevel = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0;
  int tsID;
  int varID;
  long recsize;
  off_t recpos;
  int nrecords, nrecs, recID;
  int taxisID = -1;
  taxis_t *taxis;
  int vlistID;
  SRVCOMPVAR compVar, compVar0;
  srvrec_t *srvp;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  streamptr->curTsID = 0;

  srvp  = streamptr->record->srvp;
  tsID  = tstepsNewEntry(streamID);
  taxis = &streamptr->tsteps[tsID].taxis;

  if ( tsID != 0 )
    Error("Internal problem! tstepsNewEntry returns %d", tsID);

  fileID = streamInqFileID(streamID);

  nrecs = 0;
  while ( TRUE )
    {
      recpos = fileGetPos(fileID);
      status = srvRead(fileID, srvp);
      if ( status != 0 )
	{
	  streamptr->ntsteps = 1;
	  break;
	}
      recsize = fileGetPos(fileID) - recpos;

      srvInqHeader(srvp, header);

      prec   = srvp->dprec;
      rcode  = header[0];
      rlevel = header[1];
      vdate  = header[2];
      vtime  = header[3];
      rxsize = header[4];
      rysize = header[5];

      param = cdiEncodeParam(rcode, 255, 255);

      if ( nrecs == 0 )
	{
	  datetime0.date = vdate;
	  datetime0.time = vtime;
	}
      else
	{
	  datetime.date = vdate;
	  datetime.time = vtime;
	  compVar.param = param;
          compVar.level = rlevel;
	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      compVar0.param = streamptr->tsteps[0].records[recID].param;
	      compVar0.level = streamptr->tsteps[0].records[recID].ilevel;

	      if ( memcmp(&compVar0, &compVar, sizeof(SRVCOMPVAR)) == 0 ) break;
	    }
	  if ( recID < nrecs ) break;
	  if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) )
	    Warning("Inconsistent verification time for code %d level %d", rcode, rlevel);
	}

      nrecs++;

      if ( CDI_Debug )
	Message("%4d%8d%4d%8d%8d%6d", nrecs, (int)recpos, rcode, rlevel, vdate, vtime);

      srvAddRecord(streamID, param, rlevel, rxsize, rysize, recsize, recpos, prec);
    }

  streamptr->rtsteps = 1;

  cdiGenVars(streamID);

  taxisID = taxisCreate(TAXIS_ABSOLUTE);
  taxis->type  = TAXIS_ABSOLUTE;
  taxis->vdate = datetime0.date;
  taxis->vtime = datetime0.time;

  vlistID = streamInqVlist(streamID);
  vlistDefTaxis(vlistID, taxisID);

  cdiCheckContents(streamID);

  nrecords = streamptr->tsteps[0].nallrecs;
  if ( nrecords < streamptr->tsteps[0].recordSize )
    {
      streamptr->tsteps[0].recordSize = nrecords;
      streamptr->tsteps[0].records =
	(record_t *) realloc(streamptr->tsteps[0].records, nrecords*sizeof(record_t));
    }

  streamptr->tsteps[0].recIDs = (int *) malloc(nrecords*sizeof(int));
  streamptr->tsteps[0].nrecs = nrecords;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[0].recIDs[recID] = recID;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamID);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = TRUE;
      streamptr->tsteps[tsID].position = recpos;
    }

  if ( streamptr->ntsteps == 1 )
    {
      if ( taxis->vdate == 0 && taxis->vtime == 0 )
	{
	  streamptr->ntsteps = 0;
	  for ( varID = 0; varID < streamptr->nvars; varID++ )
	    {
	      vlistDefVarTime(vlistID, varID, TIME_CONSTANT);
	    }
	}
    }
}

static
int srvScanTimestep2(int streamID)
{  
  int header[8];
  int status;
  int fileID;
  int param = 0;
  int rcode = 0, rlevel = 0, vdate = 0, vtime = 0;
  int tsID;
  int varID;
  long recsize;
  off_t recpos = 0;
  int nrecords, nrecs, recID, rindex;
  int nextstep;
  taxis_t *taxis;
  int vlistID;
  SRVCOMPVAR compVar, compVar0;
  srvrec_t *srvp;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  streamptr->curTsID = 1;

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);
  srvp    = streamptr->record->srvp;

  tsID = streamptr->rtsteps;
  if ( tsID != 1 )
    Error("Internal problem! unexpeceted timestep %d", tsID+1);

  taxis = &streamptr->tsteps[tsID].taxis;

  fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

  cdiCreateRecords(streamID, tsID);

  nrecords = streamptr->tsteps[0].nallrecs;
  streamptr->tsteps[1].recIDs = (int *) malloc(nrecords*sizeof(int));
  streamptr->tsteps[1].nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[1].recIDs[recID] = -1;

  for ( recID = 0; recID < nrecords; recID++ )
    {
      varID = streamptr->tsteps[0].records[recID].varID;
      streamptr->tsteps[tsID].records[recID].position = 
	streamptr->tsteps[0].records[recID].position;
      streamptr->tsteps[tsID].records[recID].size     = 
	streamptr->tsteps[0].records[recID].size;
    }

  for ( rindex = 0; rindex <= nrecords; rindex++ )
    {
      recpos = fileGetPos(fileID);
      status = srvRead(fileID, srvp);
      if ( status != 0 )
	{
	  streamptr->ntsteps = 2;
	  break;
	}
      recsize = fileGetPos(fileID) - recpos;

      srvInqHeader(srvp, header);

      rcode  = header[0];
      rlevel = header[1];
      vdate  = header[2];
      vtime  = header[3];

      param = cdiEncodeParam(rcode, 255, 255);

      if ( rindex == 0 )
	{
	  taxis->type  = TAXIS_ABSOLUTE;
	  taxis->vdate = vdate;
	  taxis->vtime = vtime;
	}

      compVar.param = param;
      compVar.level = rlevel;
      nextstep = FALSE;
      for ( recID = 0; recID < nrecords; recID++ )
	{
	  compVar0.param  = streamptr->tsteps[tsID].records[recID].param;
	  compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

	  if ( memcmp(&compVar0, &compVar, sizeof(SRVCOMPVAR)) == 0 )
	    {
	      if ( streamptr->tsteps[tsID].records[recID].used )
		{
		  nextstep = TRUE;
		}
	      else
		{
		  streamptr->tsteps[tsID].records[recID].used = TRUE;
		  streamptr->tsteps[tsID].recIDs[rindex] = recID;
		}
	      break;
	    }
	}
      if ( recID == nrecords )
	{
	  Warning("Code %d level %d not found at timestep %d", rcode, rlevel, tsID+1);
	  return (CDI_EUFSTRUCT);
	}

      if ( nextstep ) break;

      if ( CDI_Debug )
	Message("%4d%8d%4d%8d%8d%6d", rindex+1, (int)recpos, rcode, rlevel, vdate, vtime);

      streamptr->tsteps[tsID].records[recID].size = recsize;

      compVar0.param  = streamptr->tsteps[tsID].records[recID].param;
      compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

      if ( memcmp(&compVar0, &compVar, sizeof(SRVCOMPVAR)) != 0 )
	{
	  Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		  tsID, recID,
		  streamptr->tsteps[tsID].records[recID].param, param,
		  streamptr->tsteps[tsID].records[recID].ilevel, rlevel);
	  return (CDI_EUFSTRUCT);
	}

      streamptr->tsteps[1].records[recID].position = recpos;
    }

  nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    {
      if ( ! streamptr->tsteps[tsID].records[recID].used )
	{
	  varID = streamptr->tsteps[tsID].records[recID].varID;
	  vlistDefVarTime(vlistID, varID, TIME_CONSTANT);
	}
      else
	{
	  nrecs++;
	}
    }
  streamptr->tsteps[tsID].nrecs = nrecs;

  streamptr->rtsteps = 2;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamID);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = TRUE;
      streamptr->tsteps[tsID].position = recpos;
    }

  return (0);
}


int srvInqContents(int streamID)
{
  int fileID;
  int status = 0;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  fileID = streamInqFileID(streamID);

  streamptr->curTsID = 0;

  srvScanTimestep1(streamID);
 
  if ( streamptr->ntsteps == -1 ) status = srvScanTimestep2(streamID);

  fileSetPos(fileID, 0, SEEK_SET);

  return (status);
}

static
int srvScanTimestep(int streamID)
{
  int header[8];
  int status;
  int fileID;
  int tsID;
  int rxsize = 0, rysize = 0;
  int param = 0;
  int rcode = 0, rlevel = 0, vdate = 0, vtime = 0;
  long recsize = 0;
  off_t recpos = 0;
  int recID;
  taxis_t *taxis;
  int rindex, nrecs = 0;
  SRVCOMPVAR compVar, compVar0;
  srvrec_t *srvp;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  if ( CDI_Debug )
    {
      Message("streamID = %d", streamID);
      Message("cts = %d", streamptr->curTsID);
      Message("rts = %d", streamptr->rtsteps);
      Message("nts = %d", streamptr->ntsteps);
    }

  if ( streamptr->rtsteps == 0 )
    Error("Internal problem! Missing contents.");

  srvp  = streamptr->record->srvp;
  tsID  = streamptr->rtsteps;
  taxis = &streamptr->tsteps[tsID].taxis;

  if ( streamptr->tsteps[tsID].recordSize == 0 )
    {
      cdiCreateRecords(streamID, tsID);

      nrecs = streamptr->tsteps[1].nrecs;

      streamptr->tsteps[tsID].nrecs = nrecs;
      streamptr->tsteps[tsID].recIDs = (int *) malloc(nrecs*sizeof(int));
      for ( recID = 0; recID < nrecs; recID++ )
	streamptr->tsteps[tsID].recIDs[recID] = streamptr->tsteps[1].recIDs[recID];

      fileID = streamInqFileID(streamID);

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

      for ( rindex = 0; rindex <= nrecs; rindex++ )
	{
	  recpos = fileGetPos(fileID);
	  status = srvRead(fileID, srvp);
	  if ( status != 0 )
	    {
	      streamptr->ntsteps = streamptr->rtsteps + 1;
	      break;
	    }
	  recsize = fileGetPos(fileID) - recpos;

	  srvInqHeader(srvp, header);

	  rcode  = header[0];
	  rlevel = header[1];
	  vdate  = header[2];
	  vtime  = header[3];
	  rxsize = header[4];
	  rysize = header[5];

	  param = cdiEncodeParam(rcode, 255, 255);

	  // if ( rindex == nrecs ) break; gcc-4.5 internal compiler error
	  if ( rindex == nrecs ) continue;
	  recID = streamptr->tsteps[tsID].recIDs[rindex];

	  if ( rindex == 0 )
	    {
	      taxis->type  = TAXIS_ABSOLUTE;
	      taxis->vdate = vdate;
	      taxis->vtime = vtime;
	    }
	  /*
	  srvCmpRecord(streamID, tsID, nrecs, recpos, param, rlevel, rxsize, rysize);
	  */
	  compVar.param  = param;
          compVar.level  = rlevel;
	  compVar0.param = streamptr->tsteps[tsID].records[recID].param;
	  compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

	  if ( memcmp(&compVar0, &compVar, sizeof(SRVCOMPVAR)) != 0 )
	    {
	      Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		      tsID, recID,
		      streamptr->tsteps[tsID].records[recID].param, param,
		      streamptr->tsteps[tsID].records[recID].ilevel, rlevel);
	      Error("Invalid, unsupported or inconsistent record structure");
	    }

	  streamptr->tsteps[tsID].records[recID].position = recpos;
	  streamptr->tsteps[tsID].records[recID].size = recsize;

	  if ( CDI_Debug )
	    Message("%4d%8d%4d%8d%8d%6d", rindex, (int)recpos, rcode, rlevel, vdate, vtime);
	}

      streamptr->rtsteps++;

      if ( streamptr->ntsteps != streamptr->rtsteps )
	{
	  tsID = tstepsNewEntry(streamID);
	  if ( tsID != streamptr->rtsteps )
	    Error("Internal error. tsID = %d", tsID);

	  streamptr->tsteps[tsID-1].next   = 1;
	  streamptr->tsteps[tsID].position = recpos;
	}

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);
      streamptr->tsteps[tsID].position = recpos;
    }

  if ( nrecs > 0 && nrecs < streamptr->tsteps[tsID].nrecs )
    {
      Warning("Incomplete timestep. Stop scanning at timestep %d.", tsID);
      streamptr->ntsteps = tsID;
    }

  return (streamptr->ntsteps);
}


int srvInqTimestep(int streamID, int tsID)
{
  int ntsteps, nrecs;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  if ( tsID == 0 && streamptr->rtsteps == 0 )
    Error("Call to cdiInqContents missing!");

  if ( CDI_Debug )
    Message("tsID = %d rtsteps = %d", tsID, streamptr->rtsteps);
  
  ntsteps = UNDEFID;
  while ( ( tsID + 1 ) > streamptr->rtsteps && ntsteps == UNDEFID )
    ntsteps = srvScanTimestep(streamID);

  if ( tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID )
    {
      nrecs = 0;
    }
  else
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
    }

  return (nrecs);
}


void srvReadVarDP(int streamID, int varID, double *data, int *nmiss)
{
  int vlistID, fileID;
  int levID, nlevs, gridID, gridsize;
  off_t recpos, currentfilepos;
  int header[8];
  int tsid;
  int recID;
  int i;
  int status;
  double missval;
  srvrec_t *srvp;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  srvp     = streamptr->record->srvp;
  vlistID  = streamInqVlist(streamID);
  fileID   = streamInqFileID(streamID);
  nlevs    = streamptr->vars[varID].nlevs;
  missval  = vlistInqVarMissval(vlistID, varID);
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  tsid     = streamptr->curTsID;

  if ( CDI_Debug )
    Message("nlevs = %d gridID = %d gridsize = %d", nlevs, gridID, gridsize);

  currentfilepos = fileGetPos(fileID);

  for (levID = 0; levID < nlevs; levID++)
    {
      recID = streamptr->vars[varID].level[levID];
      recpos = streamptr->tsteps[tsid].records[recID].position;
      fileSetPos(fileID, recpos, SEEK_SET);
      status = srvRead(fileID, srvp);
      srvInqHeader(srvp, header);
      srvInqDataDP(srvp, &data[levID*gridsize]);
    }
  fileSetPos(fileID, currentfilepos, SEEK_SET);

  *nmiss = 0;
  for ( i = 0; i < nlevs*gridsize; i++ )
    if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
      {
	data[i] = missval;
	(*nmiss)++;
      }
}


void srvReadVarSliceDP(int streamID, int varID, int levID, double *data, int *nmiss)
{
  int vlistID, fileID;
  int nlevs, gridID, gridsize;
  off_t recpos, currentfilepos;
  int header[8];
  int tsid;
  int recID;
  int i;
  int status;
  double missval;
  srvrec_t *srvp;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  srvp     = streamptr->record->srvp;
  vlistID  = streamInqVlist(streamID);
  fileID   = streamInqFileID(streamID);
  nlevs    = streamptr->vars[varID].nlevs;
  missval  = vlistInqVarMissval(vlistID, varID);
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  tsid     = streamptr->curTsID;

  if ( CDI_Debug )
    Message("nlevs = %d gridID = %d gridsize = %d",
	     nlevs, gridID, gridsize);

  currentfilepos = fileGetPos(fileID);

  recID = streamptr->vars[varID].level[levID];
  recpos = streamptr->tsteps[tsid].records[recID].position;
  fileSetPos(fileID, recpos, SEEK_SET);
  status = srvRead(fileID, srvp);
  srvInqHeader(srvp, header);
  srvInqDataDP(srvp, data);

  fileSetPos(fileID, currentfilepos, SEEK_SET);

  *nmiss = 0;
  for ( i = 0; i < gridsize; i++ )
    if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
      {
	data[i] = missval;
	(*nmiss)++;
      }
}


void srvWriteVarDP(int streamID, int varID, const double *data)
{
  int fileID;
  int levID, nlevs, gridID, gridsize;
  int zaxisID;
  double level;
  int header[8];
  int xsize, ysize;
  int datatype;
  int tsID;
  int vlistID;
  int pdis, pcat, pnum;
  srvrec_t *srvp;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( CDI_Debug )
    Message("streamID = %d  varID = %d", streamID, varID);

  srvp     = streamptr->record->srvp;
  vlistID  = streamInqVlist(streamID);
  fileID   = streamInqFileID(streamID);
  tsID     = streamptr->curTsID;
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  zaxisID  = vlistInqVarZaxis(vlistID, varID);
  nlevs    = zaxisInqSize(zaxisID);

  if ( CDI_Debug )
    Message("nlevs = %d gridID = %d gridsize = %d", nlevs, gridID, gridsize);

  cdiDecodeParam(vlistInqVarParam(vlistID, varID), &pnum, &pcat, &pdis);

  header[0] = pnum;
  header[2] = streamptr->tsteps[tsID].taxis.vdate;
  header[3] = streamptr->tsteps[tsID].taxis.vtime;

  xsize = gridInqXsize(gridID);
  ysize = gridInqYsize(gridID);
  if ( xsize == 0 || ysize == 0 )
    {
      xsize = gridInqSize(gridID);
      ysize = 1;
    }
  if ( gridInqType(gridID) == GRID_UNSTRUCTURED ) ysize = 1;
  if ( gridInqSize(gridID) != xsize*ysize )
    Error("Internal problem with gridsize!");

  header[4] = xsize;
  header[5] = ysize;
  header[6] = 0;
  header[7] = 0;

  datatype = vlistInqVarDatatype(vlistID, varID);

  srvp->dprec = srvDefDatatype(datatype);

  for ( levID = 0; levID < nlevs; levID++ )
    {
      level = zaxisInqLevel(zaxisID, levID);

      header[1] = (int) level;
      srvDefHeader(srvp, header);
      srvDefDataDP(srvp, &data[levID*gridsize]);
      srvWrite(fileID, srvp);
    }
}


void srvWriteVarSliceDP(int streamID, int varID, int levID, const double *data)
{
  int fileID;
  int gridID;
  int zaxisID;
  double level;
  int header[8];
  int xsize, ysize;
  int datatype;
  int tsID;
  int vlistID;
  int pdis, pcat, pnum;
  srvrec_t *srvp;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  srvp     = streamptr->record->srvp;
  vlistID  = streamInqVlist(streamID);
  fileID   = streamInqFileID(streamID);
  tsID     = streamptr->curTsID;
  gridID   = vlistInqVarGrid(vlistID, varID);
  zaxisID  = vlistInqVarZaxis(vlistID, varID);
  level    = zaxisInqLevel(zaxisID, levID);

  if ( CDI_Debug )
    Message("gridID = %d zaxisID = %d", gridID, zaxisID);

  cdiDecodeParam(vlistInqVarParam(vlistID, varID), &pnum, &pcat, &pdis);

  header[0] = pnum;
  header[1] = (int) level;
  header[2] = streamptr->tsteps[tsID].taxis.vdate;
  header[3] = streamptr->tsteps[tsID].taxis.vtime;

  xsize = gridInqXsize(gridID);
  ysize = gridInqYsize(gridID);
  if ( xsize == 0 || ysize == 0 )
    {
      xsize = gridInqSize(gridID);
      ysize = 1;
    }
  if ( gridInqType(gridID) == GRID_UNSTRUCTURED ) ysize = 1;
  if ( gridInqSize(gridID) != xsize*ysize )
    Error("Internal problem with gridsize!");

  header[4] = xsize;
  header[5] = ysize;
  header[6] = 0;
  header[7] = 0;

  datatype = vlistInqVarDatatype(vlistID, varID);

  srvp->dprec = srvDefDatatype(datatype);

  srvDefHeader(srvp, header);
  srvDefDataDP(srvp, data);
  srvWrite(fileID, srvp);
}

#endif /* HAVE_LIBSERVICE */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
