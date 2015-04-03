#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "cdi.h"

int main(int argc, const char **argv)
{
  const char *fname = "test.nc";
  int countMissingValues = 1;
  /* todo: handle optional arguments here to increase test coverage */
  if (argc)
    fname = argv[1];

  int streamID = streamOpenRead(fname);
  if (streamID < 0)
    {
      fprintf(stderr, "Open failed for file %s: %s\n",
              fname, cdiStringError(streamID));
      return EXIT_FAILURE;
    }
  int vlistID = streamInqVlist(streamID);
  size_t nVars = (size_t)vlistNvars(vlistID);

  double *buf = NULL;
  size_t bufSize = 0;
  size_t allNmissSum = 0;

  for (int tsID = 0; streamInqTimestep(streamID, tsID); ++tsID)
    {
      for (size_t varID = 0; varID < nVars; ++varID)
        {
          size_t memSize = (size_t)vlistInqVarSize(vlistID, varID)
            * sizeof (double);
          int nmiss;
          if (memSize > bufSize)
            {
              double *temp = realloc(buf, memSize);
              if (!temp)
                {
                  perror("read buffer reallocation failed");
                  return EXIT_FAILURE;
                }
              buf = temp;
            }
          streamReadVar(streamID, (int)varID, buf, &nmiss);
          allNmissSum += (size_t)nmiss;
        }
      ++tsID;
    }
  if (countMissingValues)
    printf("missing values count = %zu\n", allNmissSum);
  streamClose(streamID);
  return EXIT_SUCCESS;
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
