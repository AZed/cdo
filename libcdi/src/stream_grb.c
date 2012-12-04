#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>

#include "dmemory.h"
#include "cdi.h"
#include "stream_int.h"
#include "stream_cgribex.h"
#include "stream_gribapi.h"
#include "file.h"
#include "cgribex.h"  /* gribZip gribGetZip gribGinfo */
#include "gribapi.h"



int grib1ltypeToZaxisType(int grib_ltype)
{
  int zaxistype = ZAXIS_GENERIC;

  switch ( grib_ltype )
    {
    case GRIB1_LTYPE_SURFACE:         { zaxistype = ZAXIS_SURFACE;           break; }
    case GRIB1_LTYPE_TOA:             { zaxistype = ZAXIS_TOA;               break; }
    case GRIB1_LTYPE_SEA_BOTTOM:      { zaxistype = ZAXIS_SEA_BOTTOM;        break; }
    case GRIB1_LTYPE_ATMOSPHERE:      { zaxistype = ZAXIS_ATMOSPHERE;        break; }
    case GRIB1_LTYPE_MEANSEA:         { zaxistype = ZAXIS_MEANSEA;           break; }
    case GRIB1_LTYPE_99:
    case GRIB1_LTYPE_ISOBARIC:        { zaxistype = ZAXIS_PRESSURE;          break; }
    case GRIB1_LTYPE_HEIGHT:          { zaxistype = ZAXIS_HEIGHT;            break; }
    case GRIB1_LTYPE_ALTITUDE:        { zaxistype = ZAXIS_ALTITUDE;	     break; }
    case GRIB1_LTYPE_SIGMA:
    case GRIB1_LTYPE_SIGMA_LAYER:     { zaxistype = ZAXIS_SIGMA;	     break; }
    case GRIB1_LTYPE_HYBRID:
    case GRIB1_LTYPE_HYBRID_LAYER:    { zaxistype = ZAXIS_HYBRID;	     break; }
    case GRIB1_LTYPE_LANDDEPTH:
    case GRIB1_LTYPE_LANDDEPTH_LAYER: { zaxistype = ZAXIS_DEPTH_BELOW_LAND;  break; }
    case GRIB1_LTYPE_ISENTROPIC:      { zaxistype = ZAXIS_ISENTROPIC;	     break; }
    case GRIB1_LTYPE_SEADEPTH:        { zaxistype = ZAXIS_DEPTH_BELOW_SEA;   break; }
    }

  return (zaxistype);
}


int grib2ltypeToZaxisType(int grib_ltype)
{
  int zaxistype = ZAXIS_GENERIC;

  switch ( grib_ltype )
    {
    case GRIB2_LTYPE_SURFACE:            { zaxistype = ZAXIS_SURFACE;           break; }
    case GRIB2_LTYPE_TOA:                { zaxistype = ZAXIS_TOA;               break; }
    case GRIB2_LTYPE_SEA_BOTTOM:         { zaxistype = ZAXIS_SEA_BOTTOM;        break; }
    case GRIB2_LTYPE_ATMOSPHERE:         { zaxistype = ZAXIS_ATMOSPHERE;        break; }
    case GRIB2_LTYPE_MEANSEA:            { zaxistype = ZAXIS_MEANSEA;           break; }
    case GRIB2_LTYPE_ISOBARIC:           { zaxistype = ZAXIS_PRESSURE;          break; }
    case GRIB2_LTYPE_HEIGHT:             { zaxistype = ZAXIS_HEIGHT;            break; }
    case GRIB2_LTYPE_ALTITUDE:           { zaxistype = ZAXIS_ALTITUDE;          break; }
    case GRIB2_LTYPE_SIGMA:              { zaxistype = ZAXIS_SIGMA;             break; }
    case GRIB2_LTYPE_HYBRID:
 /* case GRIB2_LTYPE_HYBRID_LAYER: */    { zaxistype = ZAXIS_HYBRID;            break; }
    case GRIB2_LTYPE_LANDDEPTH:
 /* case GRIB2_LTYPE_LANDDEPTH_LAYER: */ { zaxistype = ZAXIS_DEPTH_BELOW_LAND;  break; }
    case GRIB2_LTYPE_ISENTROPIC:         { zaxistype = ZAXIS_ISENTROPIC;        break; }
    case GRIB2_LTYPE_SEADEPTH:           { zaxistype = ZAXIS_DEPTH_BELOW_SEA;   break; }
    }

  return (zaxistype);
}


int grbBitsPerValue(int datatype)
{
  int bitsPerValue = 16;

  if ( datatype == DATATYPE_CPX32 || datatype == DATATYPE_CPX64 )
    Error("CDI/GRIB library does not support complex numbers!");

  if ( datatype != CDI_UNDEFID )
    {
      if ( datatype > 0 && datatype <= 32 )
	bitsPerValue = datatype;
      else if ( datatype == DATATYPE_FLT64 )
	bitsPerValue = 24;
      else
	bitsPerValue = 16;
    }

  return (bitsPerValue);
}


/*
int grbInqRecord(int streamID, int *varID, int *levelID)
{
  int status;

  status = cgribexInqRecord(streamID, varID, levelID);

  return (status);
}
*/

int grbDefRecord(int streamID)
{
  int fileID;
  int gridID;
  int status = 0;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  fileID = streamInqFileID(streamID);
  gridID = streamptr->record->gridID;

  return (status);
}

static
int grbDecode(int filetype, unsigned char *gribbuffer, int gribsize, double *data, int gridsize,
	      int unreduced, int *nmiss, int *zip, double missval)
{
  int status = 0;

#if  defined  (HAVE_LIBCGRIBEX)
  if ( filetype == FILETYPE_GRB )
    status = cgribexDecode(gribbuffer, gribsize, data, gridsize, unreduced, nmiss, zip, missval);
  else
#endif
    status = gribapiDecode(gribbuffer, gribsize, data, gridsize, unreduced, nmiss, zip, missval);
 
  return (status);
}


int grbReadRecord(int streamID, double *data, int *nmiss)
{
  int status = 0;
  unsigned char *gribbuffer;
  int fileID;
  int recID, vrecID, tsID, gridID, varID;
  long recsize;
  off_t recpos;
  int gridsize;
  int vlistID;
  int zip;
  int filetype;
  double missval;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  filetype = streamptr->filetype;

  gribbuffer = (unsigned char *) streamptr->record->buffer;

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);
  tsID    = streamptr->curTsID;
  vrecID  = streamptr->tsteps[tsID].curRecID;
  recID   = streamptr->tsteps[tsID].recIDs[vrecID];
  recpos  = streamptr->tsteps[tsID].records[recID].position;
  recsize = streamptr->tsteps[tsID].records[recID].size;
  varID   = streamptr->tsteps[tsID].records[recID].varID;

  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);

  streamptr->numvals += gridsize;

  fileSetPos(fileID, recpos, SEEK_SET);

  fileRead(fileID, gribbuffer, (size_t) recsize);

  missval = vlistInqVarMissval(vlistID, varID);

  grbDecode(filetype, gribbuffer, recsize, data, gridsize, streamptr->unreduced, nmiss, &zip, missval);

  streamptr->tsteps[tsID].records[recID].zip = zip;

  return (status);
}

static
int grbScanTimestep1(int streamID)
{
  int status;
  int filetype;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);
  filetype  = streamptr->filetype;

#if  defined  (HAVE_LIBCGRIBEX)
  if ( filetype == FILETYPE_GRB )
    {
      status = cgribexScanTimestep1(streamID);
    }
  else
#endif
    {
      status = gribapiScanTimestep1(streamID);
    }

  return (status);
}

static
int grbScanTimestep2(int streamID)
{
  int status;
  int filetype;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);
  filetype  = streamptr->filetype;

#if  defined  (HAVE_LIBCGRIBEX)
  if ( filetype == FILETYPE_GRB )
    {
      status = cgribexScanTimestep2(streamID);
    }
  else
#endif
    {
      status = gribapiScanTimestep2(streamID);
    }

  return (status);
}

static
int grbScanTimestep(int streamID)
{
  int status;
  int filetype;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);
  filetype  = streamptr->filetype;

#if  defined  (HAVE_LIBCGRIBEX)
  if ( filetype == FILETYPE_GRB )
    {
      status = cgribexScanTimestep(streamID);
    }
  else
#endif
    {
      status = gribapiScanTimestep(streamID);
    }

  return (status);
}


int grbInqContents(int streamID)
{
  int fileID;
  int status = 0;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  fileID = streamInqFileID(streamID);

  streamptr->curTsID = 0;

  status = grbScanTimestep1(streamID);
 
  if ( status == 0 && streamptr->ntsteps == -1 ) status = grbScanTimestep2(streamID);

  fileSetPos(fileID, 0, SEEK_SET);

  return (status);
}


int grbInqTimestep(int streamID, int tsID)
{
  int ntsteps, nrecs;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  if ( tsID == 0 && streamptr->rtsteps == 0 )
    Error("Call to cdiInqContents missing!");

  if ( CDI_Debug )
    Message("tsid = %d rtsteps = %d", tsID, streamptr->rtsteps);
  
  ntsteps = CDI_UNDEFID;
  while ( (tsID + 1) > streamptr->rtsteps && ntsteps == CDI_UNDEFID )
    {
      ntsteps = grbScanTimestep(streamID);
      if ( ntsteps == CDI_EUFSTRUCT )
	{
	  streamptr->ntsteps = streamptr->rtsteps;
	  break;
	}
    }

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


void grbReadVarDP(int streamID, int varID, double *data, int *nmiss)
{
  int fileID;
  int levelID, nlevs, gridID, gridsize;
  unsigned char *gribbuffer;
  int tsID, recID;
  long recsize;
  off_t recpos, currentfilepos;
  int imiss;
  int vlistID;
  int zip;
  int filetype;
  double missval;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  filetype = streamptr->filetype;

  gribbuffer = (unsigned char *) streamptr->record->buffer;

  vlistID  = streamInqVlist(streamID);
  fileID   = streamInqFileID(streamID);
  tsID     = streamptr->curTsID;

  nlevs    = streamptr->vars[varID].nlevs;
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);

  if ( CDI_Debug )
    Message("nlevs = %d gridID = %d gridsize = %d", nlevs, gridID, gridsize);

  currentfilepos = fileGetPos(fileID);

  *nmiss = 0;
  for ( levelID = 0; levelID < nlevs; levelID++ )
    {
      recID   = streamptr->vars[varID].level[levelID];
      recpos  = streamptr->tsteps[tsID].records[recID].position;
      recsize = streamptr->tsteps[tsID].records[recID].size;

      fileSetPos(fileID, recpos, SEEK_SET);

      fileRead(fileID, gribbuffer, recsize);

      missval = vlistInqVarMissval(vlistID, varID);

      grbDecode(filetype, gribbuffer, recsize, &data[levelID*gridsize], gridsize, 
		streamptr->unreduced, &imiss, &zip, missval);

      *nmiss += imiss;

      streamptr->tsteps[tsID].records[recID].zip = zip;
    }

  fileSetPos(fileID, currentfilepos, SEEK_SET);
}


void grbReadVarSliceDP(int streamID, int varID, int levelID, double *data, int *nmiss)
{
  int fileID;
  int gridID, gridsize;
  unsigned char *gribbuffer;
  long recsize;
  off_t recpos, currentfilepos;
  int tsID, recID;
  int vlistID;
  int zip;
  int filetype;
  double missval;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  filetype = streamptr->filetype;

  gribbuffer = (unsigned char *) streamptr->record->buffer;

  vlistID  = streamInqVlist(streamID);
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  tsID     = streamptr->curTsID;

  if ( CDI_Debug )
    Message("gridID = %d gridsize = %d", gridID, gridsize);

  fileID = streamInqFileID(streamID);

  currentfilepos = fileGetPos(fileID);

  recID   = streamptr->vars[varID].level[levelID];
  recpos  = streamptr->tsteps[tsID].records[recID].position;
  recsize = streamptr->tsteps[tsID].records[recID].size;

  if ( recsize == 0 )
    Error("Internal problem! Recordsize is zero for record %d at timestep %d",
	  recID+1, tsID+1);

  fileSetPos(fileID, recpos, SEEK_SET);

  fileRead(fileID, gribbuffer, recsize);

  missval = vlistInqVarMissval(vlistID, varID);

  grbDecode(filetype, gribbuffer, recsize, data, gridsize, streamptr->unreduced, nmiss, &zip, missval);

  fileSetPos(fileID, currentfilepos, SEEK_SET);

  streamptr->tsteps[tsID].records[recID].zip = zip;
}

static
size_t grbEncode(int filetype, int varID, int levelID, int vlistID, int gridID, int zaxisID,
		 int date, int time, int tsteptype, int numavg, 
		 long datasize, const double *data, int nmiss, unsigned char **gribbuffer,
		 int ljpeg, void *gribContainer)
{
  size_t nbytes;
  size_t gribbuffersize;

#if  defined  (HAVE_LIBCGRIBEX)
  if ( filetype == FILETYPE_GRB )
    {
      gribbuffersize = datasize*4+3000;
      *gribbuffer = (unsigned char *) malloc(gribbuffersize);

      nbytes = cgribexEncode(varID, levelID, vlistID, gridID, zaxisID,
			     date, time, tsteptype, numavg, 
			     datasize, data, nmiss, *gribbuffer, gribbuffersize);
    }
  else
#endif
    {
      nbytes = gribapiEncode(varID, levelID, vlistID, gridID, zaxisID,
			     date, time, tsteptype, numavg, 
			     datasize, data, nmiss, gribbuffer, &gribbuffersize,
			     ljpeg, gribContainer);
    }

  return (nbytes);
}

static
size_t grbSzip(int filetype, unsigned char *gribbuffer, size_t gribbuffersize)
{
  size_t nbytes = 0;
  unsigned char *buffer;
  size_t buffersize;
  static int lszip_warn = 1;

  buffersize = gribbuffersize + 1000; /* compressed record can be greater than source record */
  buffer = (unsigned char *) malloc(buffersize);

  /*  memcpy(buffer, gribbuffer, gribbuffersize); */
  
  if ( filetype == FILETYPE_GRB )
    {
      nbytes = gribZip(gribbuffer, (long) gribbuffersize, buffer, (long) buffersize);
    }
  else
    {
      if ( lszip_warn ) Warning("Szip compression of GRIB2 records not implemented!");
      lszip_warn = 0;
      nbytes = gribbuffersize;
    }
      
  free(buffer);

  return (nbytes);
}


int grbWriteVarSliceDP(int streamID, int varID, int levelID, const double *data, int nmiss)
{
  size_t nwrite;
  int fileID;
  int gridID;
  int zaxisID;
  unsigned char *gribbuffer = NULL;
  long datasize;
  int tsID;
  int vlistID;
  int date, time;
  int tsteptype;
  int numavg = 0;
  size_t nbytes;
  int filetype;
  stream_t *streamptr;
  int ljpeg = 0;
  int ljpeg_warn = 1;
  void *gc = NULL;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  filetype  = streamptr->filetype;

  fileID    = streamInqFileID(streamID);
  vlistID   = streamInqVlist(streamID);
  gridID    = vlistInqVarGrid(vlistID, varID);
  zaxisID   = vlistInqVarZaxis(vlistID, varID);
  tsteptype = vlistInqVarTsteptype(vlistID, varID);

  tsID      = streamptr->curTsID;
  date      = streamptr->tsteps[tsID].taxis.vdate;
  time      = streamptr->tsteps[tsID].taxis.vtime;
  if ( vlistInqVarTimave(vlistID, varID) )
    numavg = streamptr->tsteps[tsID].taxis.numavg;

  if ( CDI_Debug )
    Message("gridID = %d zaxisID = %d", gridID, zaxisID);

  datasize = gridInqSize(gridID);
  /*
  gribbuffersize = datasize*4+3000;
  gribbuffer = (unsigned char *) malloc(gribbuffersize);
  */
#if  defined  (HAVE_LIBCGRIBEX)
  if ( filetype == FILETYPE_GRB )
    {
    }
  else
#endif
    {
      gribContainer_t *gribContainers =  (gribContainer_t *) streamptr->gribContainers;
      gc = (void *) &gribContainers[varID];
    }

  if ( streamptr->comptype == COMPRESS_JPEG )
    {
      if ( filetype == FILETYPE_GRB2 )
	{
	  ljpeg = 1;
	}
      else
	{
	  if ( ljpeg_warn ) Warning("JPEG compression of GRIB1 records not available!");
	  ljpeg_warn = 0;
	}
    }

  nbytes = grbEncode(filetype, varID, levelID, vlistID, gridID, zaxisID, date, time, tsteptype, numavg, 
		     datasize, data, nmiss, &gribbuffer, ljpeg, gc);

  if ( streamptr->comptype == COMPRESS_SZIP )
    nbytes = grbSzip(filetype, gribbuffer, nbytes);

  nwrite = fileWrite(fileID, gribbuffer, nbytes);
  if ( nwrite != nbytes ) perror(__func__);

  if ( gribbuffer ) free(gribbuffer);

  return ((int)nwrite);
}


void grbWriteVarDP(int streamID, int varID, const double *data, int nmiss)
{
  int vlistID, gridID, zaxisID, levelID, nlevs;
  int gridsize;

  vlistID  = streamInqVlist(streamID);
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  zaxisID  = vlistInqVarZaxis(vlistID, varID);
  nlevs    = zaxisInqSize(zaxisID);

  for ( levelID = 0; levelID < nlevs; levelID++ )
    {
      grbWriteVarSliceDP(streamID, varID, levelID, data+levelID*gridsize, nmiss);
    }
}


int grbCopyRecord(int streamID2, int streamID1)
{
  int fileID1, fileID2;
  int tsID, recID, vrecID;
  long recsize;
  size_t gribbuffersize;
  off_t recpos;
  size_t nwrite;
  unsigned char *gribbuffer;
  int filetype;
  size_t nbytes;
  long unzipsize;
  int izip;
  stream_t *streamptr1;
  stream_t *streamptr2;

  streamptr1 = stream_to_pointer(streamID1);
  streamptr2 = stream_to_pointer(streamID2);

  stream_check_ptr(__func__, streamptr1);
  stream_check_ptr(__func__, streamptr2);

  filetype = streamptr1->filetype;

  fileID1 = streamInqFileID(streamID1);
  fileID2 = streamInqFileID(streamID2);

  tsID    = streamptr1->curTsID;
  vrecID  = streamptr1->tsteps[tsID].curRecID;
  recID   = streamptr1->tsteps[tsID].recIDs[vrecID];
  recpos  = streamptr1->tsteps[tsID].records[recID].position;
  recsize = streamptr1->tsteps[tsID].records[recID].size;

  fileSetPos(fileID1, recpos, SEEK_SET);

  gribbuffersize = recsize == (recsize>>3)<<3 ? recsize : (1+(recsize>>3))<<3;

  gribbuffer = (unsigned char *) malloc(gribbuffersize);

  fileRead(fileID1, gribbuffer, recsize);

  nbytes = recsize;

  izip = gribGetZip(recsize, gribbuffer, &unzipsize);
  
  if ( izip == 0 )
    if ( streamptr2->comptype == COMPRESS_SZIP )
      nbytes = grbSzip(filetype, gribbuffer, nbytes);

  while ( nbytes & 7 ) gribbuffer[nbytes++] = 0;

  nwrite = fileWrite(fileID2, gribbuffer, nbytes);
  if ( nwrite != nbytes ) perror(__func__);

  free(gribbuffer);

  return ((int)nwrite);
}


int grbWriteRecord(int streamID, const double *data, int nmiss)
{
  int status = 0;
  int varID, levelID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  varID   = streamptr->record->varID;
  levelID = streamptr->record->levelID;

  status = grbWriteVarSliceDP(streamID, varID, levelID, data, nmiss);

  return (status);
}


void streamInqGinfo(int streamID, int *intnum, float *fltnum)
{
  int recID, vrecID, tsID;
  int filetype;
  void *gribbuffer;
  long recsize;
  long gribbuffersize;
  off_t recpos;
  int zip;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  filetype = streamptr->filetype;

  if ( filetype == FILETYPE_GRB )
    {
      tsID    = streamptr->curTsID;
      vrecID  = streamptr->tsteps[tsID].curRecID;
      recID   = streamptr->tsteps[tsID].recIDs[vrecID];
      recpos  = streamptr->tsteps[tsID].records[recID].position;
      recsize = streamptr->tsteps[tsID].records[recID].size;
      zip     = streamptr->tsteps[tsID].records[recID].zip;

      gribbuffer = streamptr->record->buffer;
      gribbuffersize = streamptr->record->buffersize;

      if ( zip > 0 )
	Error("Compressed GRIB records unsupported!");
      else
	gribGinfo(recpos, gribbuffersize, (unsigned char *) gribbuffer, intnum, fltnum);
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
