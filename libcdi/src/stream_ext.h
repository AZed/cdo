#ifndef _STREAM_EXT_H
#define _STREAM_EXT_H

#ifndef _EXTRA_H
#  include "extra.h"
#endif

int    extInqContents(int streamID);
int    extInqTimestep(int streamID, int tsID);

int    extInqRecord(int streamID, int *varID, int *levelID);
int    extDefRecord(int streamID);
int    extCopyRecord(int streamIDdest, int streamIDsrc);
int    extReadRecord(int streamID, double *data, int *nmiss);
int    extWriteRecord(int streamID, const double *data);

void   extReadVarDP (int streamID, int varID,       double *data, int *nmiss);
void   extWriteVarDP(int streamID, int varID, const double *data);

void   extReadVarSliceDP (int streamID, int varID, int levelID,       double *data, int *nmiss);
void   extWriteVarSliceDP(int streamID, int varID, int levelID, const double *data);

#endif  /* _STREAM_EXT_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
