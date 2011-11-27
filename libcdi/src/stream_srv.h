#ifndef _STREAM_SRV_H
#define _STREAM_SRV_H

#ifndef _SERVICE_H
#  include "service.h"
#endif

int    srvInqContents(int streamID);
int    srvInqTimestep(int streamID, int tsID);

int    srvInqRecord(int streamID, int *varID, int *levelID);
int    srvDefRecord(int streamID);
int    srvCopyRecord(int streamIDdest, int streamIDsrc);
int    srvReadRecord(int streamID, double *data, int *nmiss);
int    srvWriteRecord(int streamID, const double *data);

void   srvReadVarDP (int streamID, int varID,       double *data, int *nmiss);
void   srvWriteVarDP(int streamID, int varID, const double *data);

void   srvReadVarSliceDP (int streamID, int varID, int levelID,       double *data, int *nmiss);
void   srvWriteVarSliceDP(int streamID, int varID, int levelID, const double *data);

#endif  /* _STREAM_SRV_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
