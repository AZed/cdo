#ifndef _PIO_H
#define _PIO_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef USE_MPI
#include <stdlib.h>
#include <mpi.h>

#include "cdi_int.h"

void   backendCleanup  ( void );
void   backendInit     ( void );
void   backendFinalize ( void );
int pioFileOpen(const char *filename, const char *mode);
int    pioFileClose    ( int );
size_t cdiPioFileWrite(int fileID, const void *restrict buffer, size_t len,
                       int tsID);
#else
typedef int MPI_Comm;
#endif

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
