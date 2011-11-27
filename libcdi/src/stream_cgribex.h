#ifndef _STREAM_CGRIBEX_H
#define _STREAM_CGRIBEX_H

int cgribexScanTimestep1(int streamID);
int cgribexScanTimestep2(int streamID);
int cgribexScanTimestep(int streamID);

int cgribexDecode(unsigned char *gribbuffer, int gribsize, double *data, int gridsize,
		  int unreduced, int *nmiss, int *zip, double missval);

size_t cgribexEncode(int varID, int levelID, int vlistID, int gridID, int zaxisID,
		     int vdate, int vtime, int tsteptype, int numavg, 
		     long datasize, const double *data, int nmiss, unsigned char *gribbuffer, size_t gribbuffersize);

#endif  /* _STREAM_CGRIBEX_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
