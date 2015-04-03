#ifndef PIO_INTERFACE_
#define PIO_INTERFACE_

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef USE_MPI

#include <mpi.h>
#include <yaxt.h>

#include "pio_rpc.h"

void
pioBufferPartData(int streamID, int varID, const double *data,
                  int nmiss, Xt_idxlist partDesc);
void pioBufferData (int, int, const double *, int );
void pioBufferFuncCall(union winHeaderEntry header,
                       const void *data, size_t data_len);

extern float cdiPIOpartInflate_;

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
