#ifndef _STREAM_CDF_H
#define _STREAM_CDF_H

void   cdfDefVars(int streamID);
void   cdfDefTimestep(int streamID, int tsID);
int    cdfInqTimestep(int streamID, int tsID);
int    cdfInqContents(int streamID);
void   cdfDefHistory(int streamID, int size, char *history);
int    cdfInqHistorySize(int streamID);
void   cdfInqHistoryString(int streamID, char *history);

void   cdfEndDef(int streamID);
int    cdfDefRecord(int streamID);

int    cdfCopyRecord(int streamIDdest, int streamIDsrc);

int    cdfReadRecord(int streamID, double *data, int *nmiss);
void   cdf_write_record(int streamID, int memtype, const void *data, int nmiss);

void   cdfReadVarDP(int streamID, int varID, double *data, int *nmiss);
void   cdf_write_var(int streamID, int varID, int memtype, const void *data, int nmiss);

int    cdfReadVarSliceDP(int streamID, int varID, int levelID, double *data, int *nmiss);
int    cdf_write_var_slice(int streamID, int varID, int levelID, int memtype, const void *data, int nmiss);

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
