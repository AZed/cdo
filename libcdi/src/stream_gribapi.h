#ifndef _STREAM_GRIBAPI_H
#define _STREAM_GRIBAPI_H

int gribapiScanTimestep1(int streamID);
int gribapiScanTimestep2(int streamID);
int gribapiScanTimestep(int streamID);

int gribapiDecode(unsigned char *gribbuffer, int gribsize, double *data, int gridsize,
		  int unreduced, int *nmiss, int *zip, double missval);

size_t gribapiEncode(int varID, int levelID, int vlistID, int gridID, int zaxisID,
		     int vdate, int vtime, int tsteptype, int numavg, 
		     long datasize, const double *data, int nmiss, unsigned char **gribbuffer, size_t *gribbuffersize,
		     int ljpeg, void *gribContainer);

#endif  /* _STREAM_GRIBAPI_H */
