#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <ctype.h>

#include "dmemory.h"
#include "cdi.h"
#include "stream_int.h"
#include "cdf.h"
#include "stream_grb.h"
#include "stream_cdf.h"
#include "stream_srv.h"
#include "stream_ext.h"
#include "stream_ieg.h"
#include "file.h"
#include "cgribex.h"
#include "gribapi.h"
#include "cdf.h"
#include "service.h"
#include "extra.h"
#include "ieg.h"
#include "vlist.h"

#define  MAX_FNAMES  3

FILE *popen(const char *command, const char *type);
int pclose(FILE *stream);

void cdiPrintDefaults(void)
{
  fprintf (stderr, "default instID     :  %d\n", cdiDefaultInstID); 
  fprintf (stderr, "default modelID    :  %d\n", cdiDefaultModelID); 
  fprintf (stderr, "default tableID    :  %d\n", cdiDefaultTableID); 
  fprintf (stderr, "default missval    :  %g\n", cdiDefaultMissval); 
}


void cdiDebug(int level)
{
  if ( level == 1 || level &  2 ) CDI_Debug = 1;

  if ( CDI_Debug ) Message("debug level %d", level);

  if ( level == 1 || level &  4 ) memDebug(1);

  if ( level == 1 || level &  8 ) fileDebug(1);

  if ( level == 1 || level & 16 )
    {
#if  defined  (HAVE_LIBGRIB)
      gribSetDebug(1);
#endif
#if  defined  (HAVE_LIBNETCDF)
      cdfDebug(1);
#endif
#if  defined  (HAVE_LIBSERVICE)
      srvDebug(1);
#endif
#if  defined  (HAVE_LIBEXTRA)
      extDebug(1);
#endif
#if  defined  (HAVE_LIBIEG)
      iegDebug(1);
#endif
    }

  if ( CDI_Debug )
    {
      cdiPrintDefaults();
      cdiPrintDatatypes();
    }
}


#undef  IsBigendian
#define IsBigendian()  ( u_byteorder.c[sizeof(long) - 1] )


static int getByteorder(int byteswap)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  int byteorder = -1;

  if ( IsBigendian() )
    {
      if ( byteswap ) byteorder = CDI_LITTLEENDIAN;
      else            byteorder = CDI_BIGENDIAN;
    }
  else
    {
      if ( byteswap ) byteorder = CDI_BIGENDIAN;
      else            byteorder = CDI_LITTLEENDIAN;
    }

  return (byteorder);
}


static int getFiletype(const char *filename, int *byteorder)
{
  int filetype = CDI_EUFTYPE;
  int fileID;
  int swap = 0;
  int version;
  long recpos;
  char buffer[8];

  fileID = fileOpen(filename, "r");

  if ( fileID == CDI_UNDEFID )
    {
      if ( memcmp(filename, "http:", 5) == 0 )
	return (FILETYPE_NC);
      else
	return (CDI_ESYSTEM);
    }

  if ( fileRead(fileID, buffer, 8) != 8 ) return (CDI_EUFTYPE);

  fileRewind(fileID);

  if ( memcmp(buffer, "GRIB", 4) == 0 )
    {
      version = buffer[7];
      if ( version <= 1 )
	{
	  filetype = FILETYPE_GRB;
	  if ( CDI_Debug ) Message("found GRIB file = %s, version %d", filename, version);
	}
      else if ( version == 2 )
	{
	  filetype = FILETYPE_GRB2;
	  if ( CDI_Debug ) Message("found GRIB2 file = %s", filename);
	}
    }
  else if ( memcmp(buffer, "CDF\001", 4) == 0 )
    {
      filetype = FILETYPE_NC;
      if ( CDI_Debug ) Message("found CDF1 file = %s", filename);
    }
  else if ( memcmp(buffer, "CDF\002", 4) == 0 )
    {
      filetype = FILETYPE_NC2;
      if ( CDI_Debug ) Message("found CDF2 file = %s", filename);
    }
  else if ( memcmp(buffer+1, "HDF", 3) == 0 )
    {
      filetype = FILETYPE_NC4;
      if ( CDI_Debug ) Message("found HDF file = %s", filename);
    }
#if  defined  (HAVE_LIBSERVICE)
  else if ( srvCheckFiletype(fileID, &swap) )
    {
      filetype = FILETYPE_SRV;
      if ( CDI_Debug ) Message("found SRV file = %s", filename);
    }
#endif
#if  defined  (HAVE_LIBEXTRA)
  else if ( extCheckFiletype(fileID, &swap) )
    {
      filetype = FILETYPE_EXT;
      if ( CDI_Debug ) Message("found EXT file = %s", filename);
    }
#endif
#if  defined  (HAVE_LIBIEG)
  else if ( iegCheckFiletype(fileID, &swap) )
    {
      filetype = FILETYPE_IEG;
      if ( CDI_Debug ) Message("found IEG file = %s", filename);
    }
#endif
  else if ( gribCheckSeek(fileID, &recpos, &version) == 0 )
    {
      if ( version <= 1 )
	{
	  filetype = FILETYPE_GRB;
	  if ( CDI_Debug ) Message("found seeked GRIB file = %s", filename);
	}
      else if ( version == 2 )
	{
	  filetype = FILETYPE_GRB2;
	  if ( CDI_Debug ) Message("found seeked GRIB2 file = %s", filename);
	}
    }

  fileClose(fileID);

  *byteorder = getByteorder(swap);

  return (filetype);
}


int _readline_(FILE *fp, char *line, int len)
{
  int ichar, ipos = 0;

  while ( (ichar = fgetc(fp)) != EOF )
    {
      if ( ichar == '\n' ) break;
      line[ipos++] = ichar;
      if ( ipos >= len )
        {
          fprintf(stderr, "readline Warning: end of line not found (maxlen = %d)!\n", len);
          break;
        }
    }
  line[ipos] = 0;

  if ( feof(fp) && ipos == 0 ) return (0);

  return (1);
}

#define  MAX_LINE  4096

int get_fnames(const char *argument, char *fnames[], int max_fnames)
{
  int num_fnames = 0;
  int len;
  int nfiles = 0;
  int i, j;
  const char *pch;
  char line[MAX_LINE];

  len = (int) strlen(argument);
  for ( i = 0; i < len; ++i )
    if ( argument[i] == ':' ) break;

  if ( i < len )
    {
      pch = &argument[i+1];
      len -= (i+1);
      if ( len && ( memcmp(argument, "filelist:", i) == 0 || 
		    memcmp(argument, "flist:", i) == 0 ) )
	{
	  for ( i = 0; i < len; ++i ) if ( pch[i] == ',' ) nfiles++;

	  if ( nfiles == 0 )
	    {
	      FILE *fp;
	      fp = fopen(pch, "r");
	      if ( fp == NULL ) Error("Open failed on %s", pch);

	      if ( CDI_Debug )
		Message("Reading file names from %s", pch);

	      rewind(fp);

	      nfiles = 0;
	      while ( _readline_(fp, line, MAX_LINE) )
		{
		  if ( line[0] == '#' || line[0] == '\0' ||
		       line[0] == ' ' ) continue;
		  
		  if ( nfiles >= max_fnames )
		    {
		      Warning("Too many input files (limit: %d)", max_fnames);
		      break;
		    }
		  fnames[nfiles] = strdupx(line);
		  nfiles++;
		}
	      
	      fclose(fp);

	      if ( nfiles == 0 ) Error("No input file found in %s", pch);
	    }
	  else
	    {
	      char xline[65536];
	      	      
	      strcpy(xline, pch);
	      for ( i = 0; i < len; i++ ) if ( xline[i] == ',' ) xline[i] = 0;
	      
	      nfiles++;
	      if ( nfiles >= max_fnames )
		{
		  Warning("Too many input files (limit: %d)", max_fnames);
		  nfiles = max_fnames;
		}

	      i = 0;
	      for ( j = 0; j < nfiles; j++ )
		{
		  fnames[j] = strdupx(&xline[i]);
		  i += strlen(&xline[i]) + 1;
		}
	    }
	}
      else if ( len && memcmp(argument, "ls:", i) == 0 )
	{
	  char command[4096];
	  FILE *pfp;
	  
	  strcpy(command, "ls ");
	  strcat(command, pch);

	  pfp = popen(command, "r");
	  if ( pfp == NULL )
	    SysError("popen %s failed", command);
	  
	  nfiles = 0;
	  while ( _readline_(pfp, line, MAX_LINE) )
	    {
	      if ( nfiles >= max_fnames )
		{
		  Warning("Too many input files (limit: %d)", max_fnames);
		  break;
		}
	      fnames[nfiles++] = strdupx(line);
	    }

	  pclose(pfp);
	  /*
	  for ( j = 0; j < nfiles; j++ )
	    fnames[j] = fnames[j];
	  */
	}
    }

  num_fnames = nfiles;
  
  return (num_fnames);
}

/*
@Function  streamInqFiletype
@Title     Get the filetype

@Prototype int streamInqFiletype(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.

@Description
The function @func{streamInqFiletype} returns the filetype of a stream.

@Result
@func{streamInqFiletype} returns the type of the file format,
one of the set of predefined CDI file format types.
The valid CDI file format types are @func{FILETYPE_GRB}, @func{FILETYPE_GRB2}, @func{FILETYPE_NC}, @func{FILETYPE_NC2},
@func{FILETYPE_NC4}, @func{FILETYPE_NC4C}, @func{FILETYPE_SRV}, @func{FILETYPE_EXT} and @func{FILETYPE_IEG}.

@EndFunction
*/
int streamInqFiletype(int streamID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->filetype);
}


int getByteswap(int byteorder)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  int byteswap = 0;

  if ( IsBigendian() )
    {
      if ( byteorder == CDI_LITTLEENDIAN ) byteswap = TRUE;
    }
  else
    {
      if ( byteorder == CDI_BIGENDIAN ) byteswap = TRUE;
    }

  return (byteswap);
}

/*
@Function  streamDefByteorder
@Title     Define the byte order

@Prototype void streamDefByteorder(int streamID, int byteorder)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.
    @Item  byteorder The byte order of a dataset, one of the CDI constants @func{CDI_BIGENDIAN} and
                     @func{CDI_LITTLEENDIAN}.

@Description
The function @func{streamDefByteorder} defines the byte order of a binary dataset
with the file format type @func{FILETYPE_SRV}, @func{FILETYPE_EXT} or @func{FILETYPE_IEG}.

@EndFunction
*/
void streamDefByteorder(int streamID, int byteorder)
{
  int filetype, fileID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  streamptr->byteorder = byteorder;
  filetype = streamptr->filetype;
  fileID   = streamptr->fileID;

  switch (filetype)
    {
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
	srvrec_t *srvp = streamptr->record->srvp;
	srvp->byteswap = getByteswap(byteorder);

	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
	extrec_t *extp = streamptr->record->extp;
	extp->byteswap = getByteswap(byteorder);

	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
	iegrec_t *iegp = streamptr->record->iegp;
	iegp->byteswap = getByteswap(byteorder);

	break;
      }
#endif
    }
}

/*
@Function  streamInqByteorder
@Title     Get the byte order

@Prototype int streamInqByteorder(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.

@Description
The function @func{streamInqByteorder} returns the byte order of a binary dataset
with the file format type @func{FILETYPE_SRV}, @func{FILETYPE_EXT} or @func{FILETYPE_IEG}.

@Result
@func{streamInqByteorder} returns the type of the byte order.
The valid CDI byte order types are @func{CDI_BIGENDIAN} and @func{CDI_LITTLEENDIAN}

@EndFunction
*/
int streamInqByteorder(int streamID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->byteorder);
}


char *streamFilesuffix(int filetype)
{
  static char *fileSuffix[] = {"", ".grb", ".g2", ".nc", ".nc2", ".nc4", ".nc4", ".srv", ".ext", ".ieg", ".h5"};
  int size = (int) (sizeof(fileSuffix)/sizeof(char *));

  if ( filetype > 0 && filetype < size )
    return (fileSuffix[filetype]);
  else
    return (fileSuffix[0]);
}


char *streamFilename(int streamID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->filename);
}


int cdiInqTimeSize(int streamID)
{
  int ntsteps;
  int tsID = 0, nrecs;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  ntsteps = streamptr->ntsteps;

  if ( ntsteps == CDI_UNDEFID )
    while ( (nrecs = streamInqTimestep(streamID, tsID++)) )

  ntsteps = streamptr->ntsteps;

  return (ntsteps);
}


int cdiInqContents(int streamID)
{
  int filetype;
  int vlistID;
  int taxisID;
  int status = 0;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  filetype = streamptr->filetype;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        status = grbInqContents(streamID);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        status = srvInqContents(streamID);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        status = extInqContents(streamID);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        status = iegInqContents(streamID);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
        status = cdfInqContents(streamID);
	break;
      }
#endif
    default:
      {
	if ( CDI_Debug )
	  Message("%s support not compiled in!", strfiletype(filetype));

	status = CDI_ELIBNAVAIL;
      }
    }

  if ( status == 0 )
    {
      vlistID = streamInqVlist(streamID);
      taxisID = vlistInqTaxis(vlistID);
      if ( taxisID != -1 )
	ptaxisCopy(taxisPtr(taxisID), &streamptr->tsteps[0].taxis);
    }

  return (status);
}


int streamOpen(const char *filename, const char *filemode, int filetype)
{
  int fileID = CDI_UNDEFID;
  int streamID = CDI_ESYSTEM;
  int status;
  Record *record = NULL;
  stream_t *streamptr = NULL;

  if ( CDI_Debug )
    Message("Open %s mode %c file %s", strfiletype(filetype), (int) *filemode, filename);

  if ( ! filename || ! filemode || filetype < 0 ) return (CDI_EINVAL);

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        fileID = gribOpen(filename, filemode);
	record = (Record *) malloc(sizeof(Record));
	record->buffer = NULL;
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        fileID = fileOpen(filename, filemode);
	record = (Record *) malloc(sizeof(Record));
	record->buffer = NULL;
	record->srvp   = srvNew();
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        fileID = fileOpen(filename, filemode);
	record = (Record *) malloc(sizeof(Record));
	record->buffer = NULL;
	record->extp   = extNew();
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        fileID = fileOpen(filename, filemode);
	record = (Record *) malloc(sizeof(Record));
	record->buffer = NULL;
	record->iegp   = iegNew();
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
      {
	fileID = cdfOpen(filename, filemode);
	break;
      }
    case FILETYPE_NC2:
      {
	fileID = cdfOpen64(filename, filemode);
	break;
      }
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
	fileID = cdf4Open(filename, filemode, &filetype);
	break;
      }
 #endif
    default:
      {
	if ( CDI_Debug ) Message("%s support not compiled in!", strfiletype(filetype));
	return (CDI_ELIBNAVAIL);
      }
    }

  if ( fileID < 0 )
    {
      streamID = fileID;
    }
  else
    {
      streamptr = stream_new_entry();
      streamID  = streamptr->self;

      if ( streamID < 0 ) return(CDI_ELIMIT);

      streamptr->record   = record;
      streamptr->filetype = filetype;
      streamptr->filemode = tolower(*filemode);
      streamptr->filename = strdupx(filename);
      streamptr->fileID   = fileID;

      if ( streamptr->filemode == 'r' )
	{
	  vlist_t *vlistptr;
	  int vlistID;
	  vlistID = vlistCreate();
	  if ( vlistID < 0 ) return(CDI_ELIMIT);

	  streamptr->vlistID = vlistID;
	  /* cdiReadByteorder(streamID); */
	  status = cdiInqContents(streamID);
	  if ( status < 0 ) return (status);
	  vlistptr = vlist_to_pointer(streamptr->vlistID);
	  vlistptr->ntsteps = streamNtsteps(streamID);
	}
    }
 
  return (streamID);
}


int streamOpenA(const char *filename, const char *filemode, int filetype)
{
  int fileID = CDI_UNDEFID;
  int streamID = CDI_ESYSTEM;
  int status;
  Record *record = NULL;
  stream_t *streamptr = NULL;

  if ( CDI_Debug )
    Message("Open %s mode %c file %s", strfiletype(filetype), (int) *filemode, filename);

  if ( ! filename || ! filemode || filetype < 0 ) return (CDI_EINVAL);

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        fileID = gribOpen(filename, "r");
	record = (Record *) malloc(sizeof(Record));
	record->buffer = NULL;
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        fileID = fileOpen(filename, "r");
	record = (Record *) malloc(sizeof(Record));
	record->buffer = NULL;
	record->srvp   = srvNew();
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        fileID = fileOpen(filename, "r");
	record = (Record *) malloc(sizeof(Record));
	record->buffer = NULL;
	record->extp   = extNew();
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        fileID = fileOpen(filename, "r");
	record = (Record *) malloc(sizeof(Record));
	record->buffer = NULL;
	record->iegp   = iegNew();
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
      {
	fileID = cdfOpen(filename, "r");
	break;
      }
    case FILETYPE_NC2:
      {
	fileID = cdfOpen64(filename, "r");
	break;
      }
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
	fileID = cdf4Open(filename, "r", &filetype);
	break;
      }
#endif
    default:
      {
	if ( CDI_Debug ) Message("%s support not compiled in!", strfiletype(filetype));
	return (CDI_ELIBNAVAIL);
      }
    }

  if ( fileID == CDI_UNDEFID || fileID == CDI_ELIBNAVAIL )
    {
      streamID = fileID;
      return (streamID);
    }
  else
    {
      vlist_t *vlistptr;
      streamptr = stream_new_entry();
      streamID = streamptr->self;

      streamptr->record   = record;
      streamptr->filetype = filetype;
      streamptr->filemode = tolower(*filemode);
      streamptr->filename = strdupx(filename);
      streamptr->fileID   = fileID;

      streamptr->vlistID = vlistCreate();
      /* cdiReadByteorder(streamID); */
      status = cdiInqContents(streamID);
      if ( status < 0 ) return (status);
      vlistptr = vlist_to_pointer(streamptr->vlistID);
      vlistptr->ntsteps = cdiInqTimeSize(streamID);
    }
 
  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
	gribClose(fileID);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
	fileClose(fileID);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
	fileClose(fileID);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
	fileClose(fileID);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
	cdfClose(fileID);
	break;
      }
#endif
    default:
      {
	if ( CDI_Debug ) Message("%s support not compiled in!", strfiletype(filetype));
	return (CDI_ELIBNAVAIL);
      }
    }

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        fileID = gribOpen(filename, filemode);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        fileID = fileOpen(filename, filemode);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        fileID = fileOpen(filename, filemode);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        fileID = fileOpen(filename, filemode);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
      {
	fileID = cdfOpen(filename, filemode);
	streamptr->ncmode = 2;
	break;
      }
    case FILETYPE_NC2:
      {
	fileID = cdfOpen64(filename, filemode);
	streamptr->ncmode = 2;
	break;
      }
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
	fileID = cdf4Open(filename, filemode, &filetype);
	streamptr->ncmode = 2;
	break;
      }
#endif
    default:
      {
	if ( CDI_Debug ) Message("%s support not compiled in!", strfiletype(filetype));
	return (CDI_ELIBNAVAIL);
      }
    }

  if ( fileID == CDI_UNDEFID )
    streamID = CDI_UNDEFID;
  else
    streamptr->fileID   = fileID;

  return (streamID);
}

/*
@Function  streamOpenRead
@Title     Open a dataset for reading

@Prototype int streamOpenRead(const char *path)
@Parameter
    @Item  path  The name of the dataset to be read.

@Description
The function @func{streamOpenRead} opens an existing dataset for reading.

@Result
Upon successful completion @func{streamOpenRead} returns an identifier to the
open stream. Otherwise, a negative number with the error status is returned.

@Errors
@List
   @Item  CDI_ESYSTEM     Operating system error.
   @Item  CDI_EINVAL      Invalid argument.
   @Item  CDI_EUFILETYPE  Unsupported file type.
   @Item  CDI_ELIBNAVAIL  Library support not compiled in.
@EndList

@Example
Here is an example using @func{streamOpenRead} to open an existing netCDF
file named @func{foo.nc} for reading:

@Source
#include "cdi.h"
   ...
int streamID;
   ...
streamID = streamOpenRead("foo.nc");
if ( streamID < 0 ) handle_error(streamID);
   ...
@EndSource
@EndFunction
*/
int streamOpenRead(const char *filenames)
{
  int filetype, byteorder;
  int streamID;
  int num_fnames = 0;
  char *fnames[MAX_FNAMES];
  const char *filename;
  stream_t *streamptr = NULL;

  //num_fnames = get_fnames(filenames, fnames, MAX_FNAMES);

  if ( num_fnames == 0 )
    filename = filenames;
  else
    {
      int i; 
      for ( i = 0; i < num_fnames; ++i ) printf("fnames: %d %s\n", i, fnames[i]);
      filename = fnames[0];
    }

  filetype = getFiletype(filename, &byteorder);

  if ( filetype < 0 ) return (filetype);

  streamID = streamOpen(filename, "r", filetype);

  if ( streamID >= 0 )
    {
      streamptr = stream_to_pointer(streamID);
      streamptr->byteorder = byteorder;

      if ( num_fnames > 0 )
	{
	  int i;
	  streamptr->nfiles = num_fnames;
	  streamptr->fnames = (char **) malloc(num_fnames*sizeof(char *));
	  for ( i = 0; i < num_fnames; ++i )
	    streamptr->fnames[i] = fnames[i];
	}
    }

  return (streamID);
}


int streamOpenAppend(const char *filename)
{
  int filetype, byteorder;
  int streamID;
  stream_t *streamptr;

  filetype = getFiletype(filename, &byteorder);

  if ( filetype < 0 ) return (filetype);

  streamID = streamOpenA(filename, "a", filetype);

  if ( streamID >= 0 )
    {
      streamptr = stream_to_pointer(streamID);
      streamptr->byteorder = byteorder;
    }

  return (streamID);
}

/*
@Function  streamOpenWrite
@Title     Create a new dataset

@Prototype int streamOpenWrite(const char *path, int filetype)
@Parameter
    @Item  path      The name of the new dataset.
    @Item  filetype  The type of the file format, one of the set of predefined CDI file format types.
                     The valid CDI file format types are @func{FILETYPE_GRB}, @func{FILETYPE_GRB2}, @func{FILETYPE_NC},
                     @func{FILETYPE_NC2}, @func{FILETYPE_NC4}, @func{FILETYPE_NC4C}, @func{FILETYPE_SRV},
                     @func{FILETYPE_EXT} and @func{FILETYPE_IEG}.

@Description
The function @func{streamOpenWrite} creates a new datset.
@Result
Upon successful completion @func{streamOpenWrite} returns an identifier to the
open stream. Otherwise, a negative number with the error status is returned.

@Errors
@List
   @Item  CDI_ESYSTEM     Operating system error.
   @Item  CDI_EINVAL      Invalid argument.
   @Item  CDI_EUFILETYPE  Unsupported file type.
   @Item  CDI_ELIBNAVAIL  Library support not compiled in.
@EndList

@Example
Here is an example using @func{streamOpenWrite} to create a new netCDF file 
named @func{foo.nc} for writing:

@Source
#include "cdi.h"
   ...
int streamID;
   ...
streamID = streamOpenWrite("foo.nc", FILETYPE_NC);
if ( streamID < 0 ) handle_error(streamID);
   ...
@EndSource
@EndFunction
*/
int streamOpenWrite(const char *filename, int filetype)
{
  return (streamOpen(filename, "w", filetype));
}

/*
@Function  streamClose
@Title     Close an open dataset

@Prototype  void streamClose(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.

@Description
The function @func{streamClose} closes an open dataset.

@EndFunction
*/
void streamClose(int streamID)
{
  int filetype;
  int fileID;
  int index;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  if ( CDI_Debug )
    Message("fileID = %d filename = %s", streamID, streamptr->filename);

  fileID   = streamptr->fileID;
  filetype = streamptr->filetype;
  vlistID  = streamptr->vlistID;

  if ( fileID == CDI_UNDEFID )
    Warning("File %s not open!", streamptr->filename);
  else
    switch (filetype)
      {
#if  defined  (HAVE_LIBGRIB)
      case FILETYPE_GRB:
      case FILETYPE_GRB2:
	{
	  gribClose(fileID);
	  gribContainersDelete(streamID);
	  break;
	}
#endif
#if  defined  (HAVE_LIBSERVICE)
      case FILETYPE_SRV:
	{
	  fileClose(fileID);
	  srvDelete(streamptr->record->srvp);
	  break;
	}
#endif
#if  defined  (HAVE_LIBEXTRA)
      case FILETYPE_EXT:
	{
	  fileClose(fileID);
	  extDelete(streamptr->record->extp);
	  break;
	}
#endif
#if  defined  (HAVE_LIBIEG)
      case FILETYPE_IEG:
	{
	  fileClose(fileID);
	  iegDelete(streamptr->record->iegp);
	  break;
	}
#endif
#if  defined  (HAVE_LIBNETCDF)
      case FILETYPE_NC:
      case FILETYPE_NC2:
      case FILETYPE_NC4:
      case FILETYPE_NC4C:
	{
	  cdfClose(fileID);
	  break;
	}
#endif
      default:
	{
	  Error("%s support not compiled in!", strfiletype(filetype));
	  break;
	}
      }

  if ( streamptr->record )
    {
      if ( streamptr->record->buffer )
	free(streamptr->record->buffer);

      free(streamptr->record);
    }  

  streamptr->filetype = 0;
  if ( streamptr->filename ) free(streamptr->filename);

  for ( index = 0; index < streamptr->nvars; index++ )
    {
      if ( streamptr->vars[index].level )
	free(streamptr->vars[index].level);
      if ( streamptr->vars[index].lindex )
	free(streamptr->vars[index].lindex);
    }
  free(streamptr->vars);

  for ( index = 0; index < streamptr->ntsteps; ++index )
    {
      if ( streamptr->tsteps[index].records )
	free(streamptr->tsteps[index].records);
      if ( streamptr->tsteps[index].recIDs )
	free(streamptr->tsteps[index].recIDs);
    }
    
  if ( streamptr->tsteps ) free(streamptr->tsteps);

  if ( streamptr->nfiles > 0 )
    {
      for ( index = 0; index < streamptr->nfiles; ++index )
	free(streamptr->fnames[index]);

      free(streamptr->fnames);
    }

  if ( vlistID != -1 )
    {
      if ( streamptr->filemode != 'w' )
	if ( vlistInqTaxis(vlistID) != -1 )
	  {
	    taxisDestroy(vlistInqTaxis(vlistID));
	  }

      vlistDestroy(vlistID);
    }

  stream_delete_entry(streamptr);
}

/*
@Function  streamSync
@Title     Synchronize an Open Dataset to Disk

@Prototype  void streamSync(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.

@Description
The function @func{streamSync} offers a way to synchronize the disk copy of a dataset with in-memory buffers.

@EndFunction
*/
void streamSync(int streamID)
{
  int filetype;
  int fileID;
  int vlistID;
  int nvars;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  fileID   = streamptr->fileID;
  filetype = streamptr->filetype;
  vlistID  = streamInqVlist(streamID);
  nvars    = vlistNvars(vlistID);

  if ( fileID == CDI_UNDEFID )
    Warning("File %s not open!", streamptr->filename);
  else if ( vlistID == CDI_UNDEFID )
    Warning("Vlist undefined for file %s!", streamptr->filename);
  else if ( nvars == 0 )
    Warning("No variables defined!");
  else
    {
      if ( streamptr->filemode == 'w' || streamptr->filemode == 'a' )
	{
	  switch (filetype)
	    {
#if  defined  (HAVE_LIBNETCDF)
	    case FILETYPE_NC:
	    case FILETYPE_NC2:
	    case FILETYPE_NC4:
	    case FILETYPE_NC4C:
	      {
		void cdf_sync(int ncid);
		if ( streamptr->ncmode == 2 ) cdf_sync(fileID);
		break;
	      }
#endif
	    default:
	      {
		fileFlush(fileID);
		break;
	      }
	    }
	}
    }
}

/*
@Function  streamDefTimestep
@Title     Define time step

@Prototype int streamDefTimestep(int streamID, int tsID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.
    @Item  tsID      Timestep identifier.

@Description
The function @func{streamDefTimestep} defines the time step of a stream.

@Result
@func{streamDefTimestep} returns the number of records of the time step.

@EndFunction
*/
int streamDefTimestep(int streamID, int tsID)
{
  int newtsID;
  int taxisID;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  if ( CDI_Debug )
    Message("streamID = %d  tsID = %d", streamID, tsID);

  stream_check_ptr(__func__, streamptr);

  vlistID = streamInqVlist(streamID);
  taxisID = vlistInqTaxis(vlistID);
  if ( taxisID == CDI_UNDEFID )
    {
      Warning("taxisID undefined for fileID = %d! Using absolute time axis.", streamID);
      taxisID = taxisCreate(TAXIS_ABSOLUTE);
      vlistDefTaxis(vlistID, taxisID);
    }

  newtsID = tstepsNewEntry(streamID);

  if ( tsID != newtsID )
    Error("Internal problem: tsID = %d newtsID = %d", tsID, newtsID);

  streamptr->curTsID = tsID;

  ptaxisCopy(&streamptr->tsteps[tsID].taxis, taxisPtr(taxisID));

  streamptr->ntsteps = tsID + 1;

  if ( (streamptr->filetype == FILETYPE_NC  ||
	streamptr->filetype == FILETYPE_NC2 ||
	streamptr->filetype == FILETYPE_NC4 ||
	streamptr->filetype == FILETYPE_NC4C)
       && vlistHasTime(vlistID) )
    cdfDefTimestep(streamID, tsID);

  cdiCreateRecords(streamID, tsID);

  return (streamptr->ntsteps);
}

/*
@Function  streamInqTimestep
@Title     Get time step

@Prototype int streamInqTimestep(int streamID, int tsID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.
    @Item  tsID      Timestep identifier.

@Description
The function @func{streamInqTimestep} returns the time step of a stream.

@Result
@func{streamInqTimestep} returns the number of records of the time step.

@EndFunction
*/
int streamInqTimestep(int streamID, int tsID)
{
  int filetype;
  int nrecs = 0;
  int taxisID;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  vlistID = streamInqVlist(streamID);

  if ( tsID < streamptr->rtsteps )
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
      streamptr->tsteps[tsID].curRecID = CDI_UNDEFID;
      taxisID = vlistInqTaxis(vlistID);
      if ( taxisID == -1 )
	Error("Timestep undefined for fileID = %d", streamID);
      ptaxisCopy(taxisPtr(taxisID), &streamptr->tsteps[tsID].taxis);

      return (nrecs);
    }

  if ( tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID )
    {
      return (0);
    }

  filetype = streamptr->filetype;

  if ( CDI_Debug )
    Message("streamID = %d  tsID = %d  filetype = %d", streamID, tsID, filetype);

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        nrecs = grbInqTimestep(streamID, tsID);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        nrecs = srvInqTimestep(streamID, tsID);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        nrecs = extInqTimestep(streamID, tsID);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        nrecs = iegInqTimestep(streamID, tsID);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
        nrecs = cdfInqTimestep(streamID, tsID);
	break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }

  taxisID = vlistInqTaxis(vlistID);
  if ( taxisID == -1 )
    Error("Timestep undefined for fileID = %d", streamID);

  ptaxisCopy(taxisPtr(taxisID), &streamptr->tsteps[tsID].taxis);

  return (nrecs);
}

/*
@Function  streamReadVar
@Title     Read a variable

@Prototype void streamReadVar(int streamID, int varID, double *data, int *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVar reads all the values of one time step of a variable
from an open dataset.
@EndFunction
*/
void streamReadVar(int streamID, int varID, double *data, int *nmiss)
{
  int filetype;
  stream_t *streamptr;

  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamID, varID);

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
        grbReadVarDP(streamID, varID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        srvReadVarDP(streamID, varID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        extReadVarDP(streamID, varID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        iegReadVarDP(streamID, varID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
        cdfReadVarDP(streamID, varID, data, nmiss);
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

/*
@Function  streamWriteVar
@Title     Write a variable

@Prototype void streamWriteVar(int streamID, int varID, const double *data, int nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to a block of data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVar writes the values of one time step of a variable 
to an open dataset.
@EndFunction
*/
void streamWriteVar(int streamID, int varID, const double *data, int nmiss)
{
  int filetype;
  stream_t *streamptr;

  if ( CDI_Debug ) Message("streamID = %d varID = %d", streamID, varID);

  check_parg(data);

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  streamDefineTaxis(streamID);

  filetype = streamptr->filetype;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        grbWriteVarDP(streamID, varID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        srvWriteVarDP(streamID, varID, data);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        extWriteVarDP(streamID, varID, data);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        iegWriteVarDP(streamID, varID, data);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
	if ( streamptr->accessmode == 0 ) cdfEndDef(streamID);
        cdfWriteVarDP(streamID, varID, data, nmiss);
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

/*
@Function  streamReadVarSlice
@Title     Read a horizontal slice of a variable

@Prototype void streamReadVarSlice(int streamID, int varID, int levelID, double *data, int *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVar reads all the values of a horizontal slice of a variable
from an open dataset.
@EndFunction
*/
void streamReadVarSlice(int streamID, int varID, int levelID, double *data, int *nmiss)
{
  int filetype;
  int ierr = 0;
  stream_t *streamptr;

  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamID, varID);

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
        grbReadVarSliceDP(streamID, varID, levelID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        srvReadVarSliceDP(streamID, varID, levelID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        extReadVarSliceDP(streamID, varID, levelID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        iegReadVarSliceDP(streamID, varID, levelID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
        ierr = cdfReadVarSliceDP(streamID, varID, levelID, data, nmiss);
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

/*
@Function  streamWriteVarSlice
@Title     Write a horizontal slice of a variable

@Prototype void streamWriteVarSlice(int streamID, int varID, int levelID, const double *data, int nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to a block of data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVarSlice writes the values of a horizontal slice of a 
variable to an open dataset.
@EndFunction
*/
void streamWriteVarSlice(int streamID, int varID, int levelID, const double *data, int nmiss)
{
  int filetype;
  int ierr = 0;
  stream_t *streamptr;

  if ( CDI_Debug ) Message("streamID = %d varID = %d", streamID, varID);

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
        grbWriteVarSliceDP(streamID, varID, levelID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        srvWriteVarSliceDP(streamID, varID, levelID, data);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        extWriteVarSliceDP(streamID, varID, levelID, data);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        iegWriteVarSliceDP(streamID, varID, levelID, data);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
	if ( streamptr->accessmode == 0 ) cdfEndDef(streamID);
        ierr = cdfWriteVarSliceDP(streamID, varID, levelID, data, nmiss);
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


void streamWriteContents(int streamID, char *cname)
{
  FILE *cnp;
  int tsID, recID, varID, levelID;
  long recsize;
  int nrecs, nvars;
  int code, gridID, zaxisID, timeID, datatype;
  int ngrids, nzaxis;
  int filetype, gridtype;
  int xsize, ysize;
  int date, time;
  int i;
  off_t recpos, position;
  int vlistID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  vlistID = streamInqVlist(streamID);

  cnp = fopen(cname, "w");

  if ( cnp == NULL ) SysError(cname);

  fprintf(cnp, "#CDI library version %s\n", cdiLibraryVersion());
  fprintf(cnp, "#\n");

  fprintf(cnp, "filename: %s\n", streamptr->filename);
  filetype = streamptr->filetype;
  fprintf(cnp, "filetype: %s\n", strfiletype(filetype));

  fprintf(cnp, "#\n");
  fprintf(cnp, "#grids:\n");

  ngrids = vlistNgrids(vlistID);
  for ( i = 0; i < ngrids; i++ )
    {
      gridID   = vlistGrid(vlistID, i);
      gridtype = gridInqType(gridID);
      xsize    = gridInqXsize(gridID);
      ysize    = gridInqYsize(gridID);
      fprintf(cnp, "%4d:%4d:%4d:%4d\n", i+1, gridtype, xsize, ysize);
    }
  nzaxis = vlistNzaxis(vlistID);

  fprintf(cnp, "#\n");

  fprintf(cnp, "varID:code:gridID:zaxisID:timeID:datatype\n");

  nvars = vlistNvars(vlistID);
  for ( varID = 0; varID < nvars; varID++ )
    {
      code     = vlistInqVarCode(vlistID, varID);
      gridID   = vlistInqVarGrid(vlistID, varID);
      zaxisID  = vlistInqVarZaxis(vlistID, varID);
      timeID   = vlistInqVarTime(vlistID, varID);
      datatype = vlistInqVarDatatype(vlistID, varID);
      fprintf(cnp, "%4d:%4d:%4d:%4d:%4d:%4d:\n",
	      varID+1, code, gridID, zaxisID, timeID, datatype);
    }

  fprintf(cnp, "#\n");

  fprintf(cnp, "tsID:nrecs:date:time\n");
  
  tsID = 0;
  while (1)
    {
      nrecs = streamptr->tsteps[tsID].nallrecs;
      date  = streamptr->tsteps[tsID].taxis.vdate;
      time  = streamptr->tsteps[tsID].taxis.vtime;
      position = streamptr->tsteps[tsID].position;

      fprintf(cnp, "%4d:%4d:%4d:%4d:%ld\n",
	      tsID, nrecs, date, time, (long) position);

      if ( streamptr->tsteps[tsID].next )
	tsID++;
      else
	break;
    }

  fprintf(cnp, "#\n");

  fprintf(cnp, "tsID:recID:varID:levID:size:pos\n");

  tsID = 0;
  while (1)
    {
      nrecs = streamptr->tsteps[tsID].nallrecs;
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  varID   = streamptr->tsteps[tsID].records[recID].varID;
	  levelID = streamptr->tsteps[tsID].records[recID].levelID;
	  recpos  = streamptr->tsteps[tsID].records[recID].position;
	  recsize = (long)streamptr->tsteps[tsID].records[recID].size;
	  fprintf(cnp, "%4d:%4d:%4d:%4d:%4ld:%ld\n",
		  tsID, recID, varID, levelID, recsize, (long) recpos);
	}

      if ( streamptr->tsteps[tsID].next )
	tsID++;
      else
	break;
    }

  fclose(cnp);
}


void cdiDefTableID(int tableID)
{
  int modelID, instID;

  cdiDefaultTableID = tableID;

  modelID = tableInqModel(tableID);
  cdiDefaultModelID = modelID;

  instID = modelInqInstitut(modelID);
  cdiDefaultInstID = instID;
}


void cdiPrintVersion(void)
{
  fprintf(stderr, "     CDI library version : %s\n", cdiLibraryVersion());
#if  defined  (HAVE_LIBCGRIBEX)
  fprintf(stderr, " CGRIBEX library version : %s\n", cgribexLibraryVersion());
#endif
#if  defined  (HAVE_LIBGRIB_API)
  fprintf(stderr, "GRIB_API library version : %s\n", gribapiLibraryVersion());
#endif
#if  defined  (HAVE_LIBNETCDF)
  fprintf(stderr, "  netCDF library version : %s\n", cdfLibraryVersion());
#endif
#if  defined  (HAVE_LIBHDF5)
  fprintf(stderr, "    HDF5 library version : %s\n", hdfLibraryVersion());
#endif
#if  defined  (HAVE_LIBSERVICE)
  fprintf(stderr, " SERVICE library version : %s\n", srvLibraryVersion());
#endif
#if  defined  (HAVE_LIBEXTRA)
  fprintf(stderr, "   EXTRA library version : %s\n", extLibraryVersion());
#endif
#if  defined  (HAVE_LIBIEG)
  fprintf(stderr, "     IEG library version : %s\n", iegLibraryVersion());
#endif
  fprintf(stderr, "    FILE library version : %s\n", fileLibraryVersion());
}


int streamNtsteps(int streamID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->ntsteps);
}


off_t   streamNvals(int streamID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->numvals);
}

/*
@Function  streamDefVlist
@Title     Define the variable list

@Prototype void streamDefVlist(int streamID, int vlistID)
@Parameter
    @Item  streamID Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.

@Description
The function @func{streamDefVlist} defines the variable list of a stream.

@EndFunction
*/
void streamDefVlist(int streamID, int vlistID)
{
  int nvars, varID;
  int gridID, zaxisID;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  if ( streamptr->vlistID == CDI_UNDEFID )
    {
      streamptr->vlistID = vlistDuplicate(vlistID);

      nvars = vlistNvars(vlistID);
      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridID  = vlistInqVarGrid(vlistID, varID);
	  zaxisID = vlistInqVarZaxis(vlistID, varID);
	  streamNewVar(streamID, gridID, zaxisID);
	  if ( streamptr->have_missval )
	    vlistDefVarMissval(streamptr->vlistID, varID, vlistInqVarMissval(vlistID, varID));
	}

      if ( streamptr->filemode == 'w' )
	{
	  if ( streamptr->filetype == FILETYPE_NC  ||
	       streamptr->filetype == FILETYPE_NC2 ||
	       streamptr->filetype == FILETYPE_NC4 ||
	       streamptr->filetype == FILETYPE_NC4C )
	    {
	      cdfDefVars(streamID);
	    }
	  else if ( streamptr->filetype == FILETYPE_GRB  ||
		    streamptr->filetype == FILETYPE_GRB2 )
	    {
	      gribContainersNew(streamID);
	    }
	}
    }
  else
    {
      Warning("vlist already defined for %s!", streamptr->filename);
    }
}

/*
@Function  streamInqVlist
@Title     Get the variable list

@Prototype int streamInqVlist(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.

@Description
The function @func{streamInqVlist} returns the variable list of a stream.

@Result
@func{streamInqVlist} returns an identifier to the variable list.

@EndFunction
*/
int streamInqVlist(int streamID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);
  
  stream_check_ptr(__func__, streamptr);

  return (streamptr->vlistID);
}


void streamDefCompType(int streamID, int comptype)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  streamptr->comptype = comptype;
}


void streamDefCompLevel(int streamID, int complevel)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  streamptr->complevel = complevel;
}


int streamInqCompType(int streamID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->comptype);
}


int streamInqCompLevel(int streamID)
{
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->complevel);
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