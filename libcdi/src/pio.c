#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#ifdef USE_MPI

#include <ctype.h>
#include <yaxt.h>

#include "file.h"
#include "cdi_int.h"
#include "namespace.h"

#include "pio.h"
#include "cdi.h"
#include "pio_comm.h"
#include "pio_impl.h"
#include "pio_interface.h"
#include "pio_rpc.h"
#include "pio_util.h"

char * command2charP[6] = {"IO_Open_file", "IO_Close_file",
                           "IO_Get_fp","IO_Set_fp",
                           "IO_Send_buffer", "IO_Finalize"};

long initial_buffersize = 16 * 1024 * 1024;
/*  4 KB <= x < 256 MB */
/* 16 * 1024 * 1024; */
/* 16 * 1024; */
/* 4 * 1024; */

enum {
  tagKey = 100,
};

double accumProbe   = 0.0;
double accumRecv    = 0.0;
double accumSend    = 0.0;
double accumSuspend = 0.0;
double accumWait    = 0.0;
double accumWrite   = 0.0;

char *token = "%";

/***************************************************************/

int encodeFileOpTag(int ID, int sc)
{
  return ID * tagKey + sc;
}

/***************************************************************/

struct fileOpTag decodeFileOpTag(int tag)
{
  struct fileOpTag rtag;

  rtag.id = tag / tagKey;
  rtag.command = tag % tagKey;

  return rtag;
}

/***************************************************************/

size_t
cdiPioFileWrite(int fileID, const void *restrict buffer, size_t len, int tsID)
{
  size_t iret = CDI_UNDEFID;

  switch ( commInqIOMode ())
    {
    case PIO_MPI:
      iret = fwMPINONB ( fileID, tsID, buffer, len );
      break;
#ifndef _SX
    case PIO_ASYNCH:
#endif
    case PIO_WRITER:
      iret = pioSendWrite(fileID, tsID, buffer, len);
      break;
    case PIO_FPGUARD:
      iret = fwPOSIXFPGUARDSENDRECV ( fileID, tsID, buffer, len );
      break;
    }

  return iret;
}

/***************************************************************/

int pioFileClose ( int id )
{
  int iret = CDI_UNDEFID;
  switch ( commInqIOMode ())
    {
    case PIO_MPI:
      iret = fcMPINONB ( id );
      break;
#ifndef _SX
    case PIO_ASYNCH:
#endif
    case PIO_WRITER:
      iret = pioSendClose(id);
      break;
    case PIO_FPGUARD:
      iret = fcPOSIXFPGUARDSENDRECV ( id );
      break;
    }

  return iret;
}

/***************************************************************/

int pioFileOpen(const char *filename, const char *mode)
{
  int iret = CDI_UNDEFID;

  if ((mode[0] != 'w' && mode[0] != 'W') || mode[0] == 0 || mode[1] != 0)
    xabort("Unsupported mode \"%s\" in parallel file open.", mode);

  switch ( commInqIOMode ())
    {
    case PIO_MPI:
      iret = fowMPINONB ( filename );
      break;
#ifndef _SX
    case PIO_ASYNCH:
#endif
    case PIO_WRITER:
      iret = pioSendOpen(filename);
      break;
    case PIO_FPGUARD:
      iret = fowPOSIXFPGUARDSENDRECV ( filename );
      break;
    }

  return iret;
}

/***************************************************************/

void backendInit ( void )
{
  int IOMode = commInqIOMode ();

  commDefCommNode ();

  xassert ( IOMode != PIO_NONE  || commInqSizeNode () == 1 );

  switch ( IOMode )
    {
    case PIO_NONE:
      commDefCommColl ( 1 );
      commSendNodeInfo ();
      commRecvNodeMap ();
      commDefCommsIO ();
      break;
    case PIO_MPI:
      initMPINONB ();
      break;
#ifndef _SX
    case PIO_ASYNCH:
#endif
    case PIO_WRITER:
      pioSendInitialize();
      break;
    case PIO_FPGUARD:
      initPOSIXFPGUARDSENDRECV ();
      break;
    }
}

/***************************************************************/

void backendCleanup ( void )
{
  int IOMode = commInqIOMode ();
  switch ( IOMode )
    {
    case PIO_NONE:
      break;
    case PIO_MPI:
      finalizeMPINONB ();
      break;
#ifndef _SX
    case PIO_ASYNCH:
#endif
    case PIO_WRITER:
      pioSendFinalize();
      break;
    case PIO_FPGUARD:
      finalizePOSIXFPGUARDSENDRECV ();
      break;
    default:
      xdebug("%s", " BACKENDCLEANUP FUNCTION NOT IMPLEMENTED YET.");
    }
}

/***************************************************************/

void backendFinalize ( void )
{
  commDestroy ();
  MPI_Finalize ();
  exit ( EXIT_SUCCESS );
}
#endif
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
