#ifndef _STREAM_IEG_H
#define _STREAM_IEG_H

#ifndef _IEG_H
#  include "ieg.h"
#endif

int    iegInqContents(int streamID);
int    iegInqTimestep(int streamID, int tsID);

int    iegInqRecord(int streamID, int *varID, int *levelID);
int    iegDefRecord(int streamID);
int    iegCopyRecord(int streamIDdest, int streamIDsrc);
int    iegReadRecord(int streamID, double *data, int *nmiss);
int    iegWriteRecord(int streamID, const double *data);

void   iegReadVarDP (int streamID, int varID,       double *data, int *nmiss);
void   iegWriteVarDP(int streamID, int varID, const double *data);

void   iegReadVarSliceDP (int streamID, int varID, int levelID,       double *data, int *nmiss);
void   iegWriteVarSliceDP(int streamID, int varID, int levelID, const double *data);

#endif  /* _STREAM_IEG_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
