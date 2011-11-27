#include <stdio.h>
#include <stdlib.h>


#include "cdi.h"

int main(int argc, char *argv[])
{
  char fname[] = "test.grb";
  int filetype = FILETYPE_GRB;
  int nlat = 18, nlon = 2*nlat;
  double *data = NULL;
  double missval;
  int nlevel;
  int varID;
  int datasize = 0;
  int streamID1, streamID2;
  int gridID, zaxisID, code, vdate, vtime;
  int nrecs, nvars;
  int gridtype, gridsize = 0;
  int tsID;
  int timeID = 0;
  int level, levelID;
  int offset;
  int vlistID, taxisID;
  int nmiss;

  
  datasize = nlon * nlat;
  data = (double *) malloc(datasize*sizeof(double));

  gridID = gridCreate(GRID_GAUSSIAN, datasize);
  gridDefXsize(gridID, nlon);
  gridDefYsize(gridID, nlat);

  zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  vlistID = vlistCreate();
  vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);

  taxisID = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID, taxisID);

  streamID1 = streamOpenWrite(fname, filetype);
  if ( streamID1 < 0 )
    {
      fprintf(stderr, "Open failed on %s\n", fname);
      fprintf(stderr, "%s\n", cdiStringError(streamID1));
      return (-1);
    }

  streamDefVlist(streamID1, vlistID);

  (void) streamDefTimestep(streamID1, 0);

  streamWriteVar(streamID1, 0, data, 0);

  return (0);

  vlistID = streamInqVlist(streamID1);

  filetype = streamInqFiletype(streamID1);

  streamID2 = streamOpenWrite(fname, filetype);
  if ( streamID2 < 0 )
    {
      fprintf(stderr, "Open failed on %s\n", fname);
      fprintf(stderr, "%s\n", cdiStringError(streamID2));
      return (-1);
    }

  streamDefVlist(streamID2, vlistID);

  nvars = vlistNvars(vlistID);

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID, varID);
      zaxisID  = vlistInqVarZaxis(vlistID, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);
      if ( gridsize*nlevel > datasize ) datasize = gridsize*nlevel;
    }

  data = (double *) malloc(datasize*sizeof(double));

  taxisID = vlistInqTaxis(vlistID);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);

      streamDefTimestep(streamID2, tsID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  streamReadVar(streamID1, varID, data, &nmiss);

	  code     = vlistInqVarCode(vlistID, varID);
	  gridID   = vlistInqVarGrid(vlistID, varID);
	  zaxisID  = vlistInqVarZaxis(vlistID, varID);
	  gridtype = gridInqType(gridID);
	  gridsize = gridInqSize(gridID);
	  nlevel   = zaxisInqSize(zaxisID);
	  missval  = vlistInqVarMissval(vlistID, varID);

	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      level  = (int) zaxisInqLevel(zaxisID, levelID);
	      offset = gridsize*levelID;
	    }

	  streamWriteVar(streamID2, varID, data, nmiss);
	}
      tsID++;
    }

  free(data);

  streamClose(streamID2);
  streamClose(streamID1);

  return (0);
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
