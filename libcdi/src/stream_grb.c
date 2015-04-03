#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "stream_cgribex.h"
#include "stream_gribapi.h"
#include "file.h"
#include "cgribex.h"  /* gribZip gribGetZip gribGinfo */
#include "gribapi.h"
#include "namespace.h"


int grib1ltypeToZaxisType(int grib_ltype)
{
  int zaxistype = ZAXIS_GENERIC;

  switch ( grib_ltype )
    {
    case GRIB1_LTYPE_SURFACE:            { zaxistype = ZAXIS_SURFACE;                break; }
    case GRIB1_LTYPE_CLOUD_BASE:         { zaxistype = ZAXIS_CLOUD_BASE;             break; }
    case GRIB1_LTYPE_CLOUD_TOP:          { zaxistype = ZAXIS_CLOUD_TOP;              break; }
    case GRIB1_LTYPE_ISOTHERM0:          { zaxistype = ZAXIS_ISOTHERM_ZERO;          break; }
    case GRIB1_LTYPE_TOA:                { zaxistype = ZAXIS_TOA;                    break; }
    case GRIB1_LTYPE_SEA_BOTTOM:         { zaxistype = ZAXIS_SEA_BOTTOM;             break; }
    case GRIB1_LTYPE_ATMOSPHERE:         { zaxistype = ZAXIS_ATMOSPHERE;             break; }
    case GRIB1_LTYPE_MEANSEA:            { zaxistype = ZAXIS_MEANSEA;                break; }
    case GRIB1_LTYPE_99:
    case GRIB1_LTYPE_ISOBARIC:           { zaxistype = ZAXIS_PRESSURE;               break; }
    case GRIB1_LTYPE_HEIGHT:             { zaxistype = ZAXIS_HEIGHT;                 break; }
    case GRIB1_LTYPE_ALTITUDE:           { zaxistype = ZAXIS_ALTITUDE;	             break; }
    case GRIB1_LTYPE_SIGMA:
    case GRIB1_LTYPE_SIGMA_LAYER:        { zaxistype = ZAXIS_SIGMA;	             break; }
    case GRIB1_LTYPE_HYBRID:
    case GRIB1_LTYPE_HYBRID_LAYER:       { zaxistype = ZAXIS_HYBRID;	             break; }
    case GRIB1_LTYPE_LANDDEPTH:
    case GRIB1_LTYPE_LANDDEPTH_LAYER:    { zaxistype = ZAXIS_DEPTH_BELOW_LAND;       break; }
    case GRIB1_LTYPE_ISENTROPIC:         { zaxistype = ZAXIS_ISENTROPIC;             break; }
    case GRIB1_LTYPE_SEADEPTH:           { zaxistype = ZAXIS_DEPTH_BELOW_SEA;        break; }
    case GRIB1_LTYPE_LAKE_BOTTOM:        { zaxistype = ZAXIS_LAKE_BOTTOM;            break; }
    case GRIB1_LTYPE_SEDIMENT_BOTTOM:    { zaxistype = ZAXIS_SEDIMENT_BOTTOM;        break; }
    case GRIB1_LTYPE_SEDIMENT_BOTTOM_TA: { zaxistype = ZAXIS_SEDIMENT_BOTTOM_TA;     break; }
    case GRIB1_LTYPE_SEDIMENT_BOTTOM_TW: { zaxistype = ZAXIS_SEDIMENT_BOTTOM_TW;     break; }
    case GRIB1_LTYPE_MIX_LAYER:          { zaxistype = ZAXIS_MIX_LAYER;              break; }
    }

  return (zaxistype);
}


int grib2ltypeToZaxisType(int grib_ltype)
{
  int zaxistype = ZAXIS_GENERIC;

  switch ( grib_ltype )
    {
    case GRIB2_LTYPE_SURFACE:            { zaxistype = ZAXIS_SURFACE;                break; }
    case GRIB2_LTYPE_CLOUD_BASE:         { zaxistype = ZAXIS_CLOUD_BASE;             break; }
    case GRIB2_LTYPE_CLOUD_TOP:          { zaxistype = ZAXIS_CLOUD_TOP;              break; }
    case GRIB2_LTYPE_ISOTHERM0:          { zaxistype = ZAXIS_ISOTHERM_ZERO;          break; }
    case GRIB2_LTYPE_TOA:                { zaxistype = ZAXIS_TOA;                    break; }
    case GRIB2_LTYPE_SEA_BOTTOM:         { zaxistype = ZAXIS_SEA_BOTTOM;             break; }
    case GRIB2_LTYPE_ATMOSPHERE:         { zaxistype = ZAXIS_ATMOSPHERE;             break; }
    case GRIB2_LTYPE_MEANSEA:            { zaxistype = ZAXIS_MEANSEA;                break; }
    case GRIB2_LTYPE_ISOBARIC:           { zaxistype = ZAXIS_PRESSURE;               break; }
    case GRIB2_LTYPE_HEIGHT:             { zaxistype = ZAXIS_HEIGHT;                 break; }
    case GRIB2_LTYPE_ALTITUDE:           { zaxistype = ZAXIS_ALTITUDE;               break; }
    case GRIB2_LTYPE_SIGMA:              { zaxistype = ZAXIS_SIGMA;                  break; }
    case GRIB2_LTYPE_HYBRID:
 /* case GRIB2_LTYPE_HYBRID_LAYER: */    { zaxistype = ZAXIS_HYBRID;                 break; }
    case GRIB2_LTYPE_LANDDEPTH:
 /* case GRIB2_LTYPE_LANDDEPTH_LAYER: */ { zaxistype = ZAXIS_DEPTH_BELOW_LAND;       break; }
    case GRIB2_LTYPE_ISENTROPIC:         { zaxistype = ZAXIS_ISENTROPIC;             break; }
    case GRIB2_LTYPE_SNOW:               { zaxistype = ZAXIS_SNOW;                   break; }
    case GRIB2_LTYPE_SEADEPTH:           { zaxistype = ZAXIS_DEPTH_BELOW_SEA;        break; }
    case GRIB2_LTYPE_LAKE_BOTTOM:        { zaxistype = ZAXIS_LAKE_BOTTOM;            break; }
    case GRIB2_LTYPE_SEDIMENT_BOTTOM:    { zaxistype = ZAXIS_SEDIMENT_BOTTOM;        break; }
    case GRIB2_LTYPE_SEDIMENT_BOTTOM_TA: { zaxistype = ZAXIS_SEDIMENT_BOTTOM_TA;     break; }
    case GRIB2_LTYPE_SEDIMENT_BOTTOM_TW: { zaxistype = ZAXIS_SEDIMENT_BOTTOM_TW;     break; }
    case GRIB2_LTYPE_MIX_LAYER:          { zaxistype = ZAXIS_MIX_LAYER;              break; }
    case GRIB2_LTYPE_REFERENCE:          { zaxistype = ZAXIS_REFERENCE;              break; }
    }

  return (zaxistype);
}


int zaxisTypeToGrib1ltype(int zaxistype)
{
  int grib_ltype = -1;

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:               { grib_ltype = GRIB1_LTYPE_SURFACE;            break; }
    case ZAXIS_MEANSEA:               { grib_ltype = GRIB1_LTYPE_MEANSEA;            break; }
    case ZAXIS_HEIGHT:                { grib_ltype = GRIB1_LTYPE_HEIGHT;             break; }
    case ZAXIS_ALTITUDE:              { grib_ltype = GRIB1_LTYPE_ALTITUDE;           break; }
    case ZAXIS_SIGMA:                 { grib_ltype = GRIB1_LTYPE_SIGMA;              break; }
    case ZAXIS_DEPTH_BELOW_SEA:       { grib_ltype = GRIB1_LTYPE_SEADEPTH;           break; }
    case ZAXIS_ISENTROPIC:            { grib_ltype = GRIB1_LTYPE_ISENTROPIC;         break; }
    case ZAXIS_CLOUD_BASE:            { grib_ltype = GRIB1_LTYPE_CLOUD_BASE;         break; }
    case ZAXIS_CLOUD_TOP:             { grib_ltype = GRIB1_LTYPE_CLOUD_TOP;          break; }
    case ZAXIS_ISOTHERM_ZERO:         { grib_ltype = GRIB1_LTYPE_ISOTHERM0;          break; }
    case ZAXIS_TOA:                   { grib_ltype = GRIB1_LTYPE_TOA;                break; }
    case ZAXIS_SEA_BOTTOM:            { grib_ltype = GRIB1_LTYPE_SEA_BOTTOM;         break; }
    case ZAXIS_LAKE_BOTTOM:           { grib_ltype = GRIB1_LTYPE_LAKE_BOTTOM;        break; }
    case ZAXIS_SEDIMENT_BOTTOM:       { grib_ltype = GRIB1_LTYPE_SEDIMENT_BOTTOM;    break; }
    case ZAXIS_SEDIMENT_BOTTOM_TA:    { grib_ltype = GRIB1_LTYPE_SEDIMENT_BOTTOM_TA; break; }
    case ZAXIS_SEDIMENT_BOTTOM_TW:    { grib_ltype = GRIB1_LTYPE_SEDIMENT_BOTTOM_TW; break; }
    case ZAXIS_MIX_LAYER:             { grib_ltype = GRIB1_LTYPE_MIX_LAYER;          break; }
    case ZAXIS_ATMOSPHERE:            { grib_ltype = GRIB1_LTYPE_ATMOSPHERE;         break; }
    }

  return (grib_ltype);
}


int zaxisTypeToGrib2ltype(int zaxistype)
{
  int grib_ltype = -1;

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:               { grib_ltype = GRIB2_LTYPE_SURFACE;            break; }
    case ZAXIS_MEANSEA:               { grib_ltype = GRIB2_LTYPE_MEANSEA;            break; }
    case ZAXIS_HEIGHT:                { grib_ltype = GRIB2_LTYPE_HEIGHT;             break; }
    case ZAXIS_ALTITUDE:              { grib_ltype = GRIB2_LTYPE_ALTITUDE;           break; }
    case ZAXIS_SIGMA:                 { grib_ltype = GRIB2_LTYPE_SIGMA;              break; }
    case ZAXIS_DEPTH_BELOW_SEA:       { grib_ltype = GRIB2_LTYPE_SEADEPTH;           break; }
    case ZAXIS_ISENTROPIC:            { grib_ltype = GRIB2_LTYPE_ISENTROPIC;         break; }
    case ZAXIS_CLOUD_BASE:            { grib_ltype = GRIB2_LTYPE_CLOUD_BASE;         break; }
    case ZAXIS_CLOUD_TOP:             { grib_ltype = GRIB2_LTYPE_CLOUD_TOP;          break; }
    case ZAXIS_ISOTHERM_ZERO:         { grib_ltype = GRIB2_LTYPE_ISOTHERM0;          break; }
    case ZAXIS_TOA:                   { grib_ltype = GRIB2_LTYPE_TOA;                break; }
    case ZAXIS_SEA_BOTTOM:            { grib_ltype = GRIB2_LTYPE_SEA_BOTTOM;         break; }
    case ZAXIS_LAKE_BOTTOM:           { grib_ltype = GRIB2_LTYPE_LAKE_BOTTOM;        break; }
    case ZAXIS_SEDIMENT_BOTTOM:       { grib_ltype = GRIB2_LTYPE_SEDIMENT_BOTTOM;    break; }
    case ZAXIS_SEDIMENT_BOTTOM_TA:    { grib_ltype = GRIB2_LTYPE_SEDIMENT_BOTTOM_TA; break; }
    case ZAXIS_SEDIMENT_BOTTOM_TW:    { grib_ltype = GRIB2_LTYPE_SEDIMENT_BOTTOM_TW; break; }
    case ZAXIS_MIX_LAYER:             { grib_ltype = GRIB2_LTYPE_MIX_LAYER;          break; }
    case ZAXIS_ATMOSPHERE:            { grib_ltype = GRIB2_LTYPE_ATMOSPHERE;         break; }
    }

  return (grib_ltype);
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
int grbInqRecord(stream_t * streamptr, int *varID, int *levelID)
{
  int status;

  status = cgribexInqRecord(streamptr, varID, levelID);

  return (status);
}
*/

int grbDefRecord(stream_t * streamptr)
{
  int status = 0;

  return (status);
}

static
int grbDecode(int filetype, unsigned char *gribbuffer, int gribsize, double *data, int gridsize,
	      int unreduced, int *nmiss, int *zip, double missval, int vlistID, int varID)
{
  int status = 0;

#if  defined  (HAVE_LIBCGRIBEX)
  if ( filetype == FILETYPE_GRB )
    {
#if  defined  (HAVE_LIBGRIB_API)
      extern int cdiNAdditionalGRIBKeys;
      if ( cdiNAdditionalGRIBKeys > 0 )
	Error("CGRIBEX decode does not support reading of additional GRIB keys!");
#endif
      status = cgribexDecode(gribbuffer, gribsize, data, gridsize, unreduced, nmiss, zip, missval);
    }
  else
#endif
    {
      status = gribapiDecode(gribbuffer, gribsize, data, gridsize, unreduced, nmiss, zip, missval, vlistID, varID);
    }

  return (status);
}


int grbReadRecord(stream_t * streamptr, double *data, int *nmiss)
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

  filetype = streamptr->filetype;

  gribbuffer = (unsigned char *) streamptr->record->buffer;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;
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

  grbDecode(filetype, gribbuffer, recsize, data, gridsize, streamptr->unreduced, nmiss, &zip, missval, vlistID, varID);

  streamptr->tsteps[tsID].records[recID].zip = zip;

  return (status);
}

static
int grbScanTimestep1(stream_t * streamptr)
{
  int status;
  int filetype;

  filetype  = streamptr->filetype;

#if  defined  (HAVE_LIBCGRIBEX)
  if ( filetype == FILETYPE_GRB )
    {
      status = cgribexScanTimestep1(streamptr);
    }
  else
#endif
    {
      status = gribapiScanTimestep1(streamptr);
    }

  return (status);
}

static
int grbScanTimestep2(stream_t * streamptr)
{
  int status;
  int filetype;

  filetype  = streamptr->filetype;

#if  defined  (HAVE_LIBCGRIBEX)
  if ( filetype == FILETYPE_GRB )
    {
      status = cgribexScanTimestep2(streamptr);
    }
  else
#endif
    {
      status = gribapiScanTimestep2(streamptr);
    }

  return (status);
}

static
int grbScanTimestep(stream_t * streamptr)
{
  int status;
  int filetype;

  filetype  = streamptr->filetype;

#if  defined  (HAVE_LIBCGRIBEX)
  if ( filetype == FILETYPE_GRB )
    {
      status = cgribexScanTimestep(streamptr);
    }
  else
#endif
    {
      status = gribapiScanTimestep(streamptr);
    }

  return (status);
}


int grbInqContents(stream_t * streamptr)
{
  int fileID;
  int status = 0;

  fileID = streamptr->fileID;

  streamptr->curTsID = 0;

  status = grbScanTimestep1(streamptr);

  if ( status == 0 && streamptr->ntsteps == -1 ) status = grbScanTimestep2(streamptr);

  fileSetPos(fileID, 0, SEEK_SET);

  return (status);
}


int grbInqTimestep(stream_t * streamptr, int tsID)
{
  int ntsteps, nrecs;

  if ( tsID == 0 && streamptr->rtsteps == 0 )
    Error("Call to cdiInqContents missing!");

  if ( CDI_Debug )
    Message("tsid = %d rtsteps = %d", tsID, streamptr->rtsteps);

  ntsteps = CDI_UNDEFID;
  while ( (tsID + 1) > streamptr->rtsteps && ntsteps == CDI_UNDEFID )
    {
      ntsteps = grbScanTimestep(streamptr);
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


void grbReadVarDP(stream_t * streamptr, int varID, double *data, int *nmiss)
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

  filetype = streamptr->filetype;

  gribbuffer = (unsigned char *) streamptr->record->buffer;

  vlistID  = streamptr->vlistID;
  fileID   = streamptr->fileID;
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
		streamptr->unreduced, &imiss, &zip, missval, vlistID, varID);

      *nmiss += imiss;

      streamptr->tsteps[tsID].records[recID].zip = zip;
    }

  fileSetPos(fileID, currentfilepos, SEEK_SET);
}


void grbReadVarSliceDP(stream_t * streamptr, int varID, int levelID, double *data, int *nmiss)
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

  filetype = streamptr->filetype;

  gribbuffer = (unsigned char *) streamptr->record->buffer;

  vlistID  = streamptr->vlistID;
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  tsID     = streamptr->curTsID;

  if ( CDI_Debug )
    Message("gridID = %d gridsize = %d", gridID, gridsize);

  fileID = streamptr->fileID;

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

  grbDecode(filetype, gribbuffer, recsize, data, gridsize, streamptr->unreduced, nmiss, &zip, missval, vlistID, varID);

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


int grb_write_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, const void *data, int nmiss)
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
  int ljpeg = 0;
  int ljpeg_warn = 1;
  void *gc = NULL;

  if ( memtype == MEMTYPE_FLOAT ) Error("grb_write_var_slice not implemented for memtype float!");

  filetype  = streamptr->filetype;
  fileID    = streamptr->fileID;
  vlistID   = streamptr->vlistID;
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
#if defined (GRIBCONTAINER2D)
      gribContainer_t **gribContainers =  (gribContainer_t **) streamptr->gribContainers;
      gc = (void *) &gribContainers[varID][levelID];
#else
      gribContainer_t *gribContainers =  (gribContainer_t *) streamptr->gribContainers;
      gc = (void *) &gribContainers[varID];
#endif
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
		     datasize, (const double*) data, nmiss, &gribbuffer, ljpeg, gc);

  if ( streamptr->comptype == COMPRESS_SZIP )
    nbytes = grbSzip(filetype, gribbuffer, nbytes);

  {
    size_t (*myFileWrite)(int fileID, const void *restrict buffer,
                          size_t len, int tsID)
      = (size_t (*)(int, const void *restrict, size_t, int))
      namespaceSwitchGet(NSSWITCH_FILE_WRITE).func;
    nwrite = myFileWrite(fileID, gribbuffer, nbytes, tsID);
  }

  if ( nwrite != nbytes ) perror(__func__);

  if ( gribbuffer ) free(gribbuffer);

  return ((int)nwrite);
}


void grb_write_var(stream_t *streamptr, int varID, int memtype, const void *data, int nmiss)
{
  int vlistID, gridID, zaxisID, levelID, nlevs;
  int gridsize;

  vlistID  = streamptr->vlistID;
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  zaxisID  = vlistInqVarZaxis(vlistID, varID);
  nlevs    = zaxisInqSize(zaxisID);

  for ( levelID = 0; levelID < nlevs; levelID++ )
    {
      if ( memtype == MEMTYPE_FLOAT )
        grb_write_var_slice(streamptr, varID, levelID, memtype, ((float*)data)+levelID*gridsize, nmiss);
      else
        grb_write_var_slice(streamptr, varID, levelID, memtype, ((double*)data)+levelID*gridsize, nmiss);
    }
}


int grbCopyRecord(stream_t * streamptr2, stream_t * streamptr1)
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

  filetype = streamptr1->filetype;

  fileID1 = streamptr1->fileID;
  fileID2 = streamptr2->fileID;

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


int grb_write_record(stream_t * streamptr, int memtype, const void *data, int nmiss)
{
  int status = 0;
  int varID, levelID;

  varID   = streamptr->record->varID;
  levelID = streamptr->record->levelID;

  status = grb_write_var_slice(streamptr, varID, levelID, memtype, data, nmiss);

  return (status);
}


void streamInqGinfo(int streamID, int *intnum, float *fltnum, off_t *bignum)
{
  int recID, vrecID, tsID;
  int filetype;
  void *gribbuffer;
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
      zip     = streamptr->tsteps[tsID].records[recID].zip;

      gribbuffer = streamptr->record->buffer;
      gribbuffersize = streamptr->record->buffersize;

      if ( zip > 0 )
	Error("Compressed GRIB records unsupported!");
      else
	gribGinfo(recpos, gribbuffersize, (unsigned char *) gribbuffer, intnum, fltnum, bignum);
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
