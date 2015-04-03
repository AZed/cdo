#include <stdio.h>
#include "cdi.h"

int nlon = 12; // Number of longitudes
int nlat =  6; // Number of latitudes
int nlev =  5; // Number of levels
int nts  =  3; // Number of time steps

int main(void)
{
  int taxisID, vlistID, varID1, varID2, streamID, tsID;
  int nmiss, vdate, vtime;
  double var1[nlon*nlat];
  double var2[nlon*nlat*nlev];


  // Open the dataset 
  streamID = streamOpenRead("example.nc");
  if ( streamID < 0 )
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      return(1);
    }

  // Get the variable list of the dataset 
  vlistID = streamInqVlist(streamID);

  // Set the variable IDs 
  varID1 = 0;
  varID2 = 1;

  // Get the Time axis from the variable list 
  taxisID = vlistInqTaxis(vlistID);

  // Loop over the number of time steps 
  for ( tsID = 0; tsID < nts; tsID++ )
    {
      // Inquire the time step 
      streamInqTimestep(streamID, tsID);

      // Get the verification date and time 
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);

      // Read var1 and var2 
      streamReadVar(streamID, varID1, var1, &nmiss);
      streamReadVar(streamID, varID2, var2, &nmiss);
    }

  // Close the input stream 
  streamClose(streamID);

  return 0;
}
