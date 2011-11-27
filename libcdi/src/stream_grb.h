#ifndef _STREAM_GRB_H
#define _STREAM_GRB_H

int   grbBitsPerValue(int datatype);

int   grbInqContents(int streamID);
int   grbInqTimestep(int streamID, int tsID);

int   grbInqRecord(int streamID, int *varID, int *levelID);
int   grbDefRecord(int streamID);
int   grbWriteRecord(int streamID, const double *data, int nmiss);
int   grbReadRecord(int streamID, double *data, int *nmiss);
int   grbCopyRecord(int streamIDdest, int streamIDsrc);

void  grbReadVarDP(int streamID, int varID, double *data, int *nmiss);
void  grbWriteVarDP(int streamID, int varID, const double *data, int nmiss);

void  grbReadVarSliceDP(int streamID, int varID, int levelID, double *data, int *nmiss);
int   grbWriteVarSliceDP(int streamID, int varID, int levelID, const double *data, int nmiss);

int   grib1ltypeToZaxisType(int grib_ltype);
int   grib2ltypeToZaxisType(int grib_ltype);

#endif  /* _STREAM_GRB_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
