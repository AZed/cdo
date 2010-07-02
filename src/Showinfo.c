/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2010 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

      Showinfo   showparam       Show parameters
      Showinfo   showcode        Show code numbers
      Showinfo   showname        Show variable names
      Showinfo   showstdname     Show variable standard names
      Showinfo   showlevel       Show levels
      Showinfo   showyear        Show years
      Showinfo   showmon         Show months
      Showinfo   showdate        Show dates
      Showinfo   showtime        Show timesteps
      Showinfo   showltype       Show level types
      Showinfo   showformat      Show file format
*/


#include <stdio.h>
#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Showinfo(void *argument)
{
  int SHOWYEAR, SHOWMON, SHOWDATE, SHOWTIME, SHOWTIMESTAMP, SHOWCODE;
  int SHOWPARAM, SHOWNAME, SHOWSTDNAME, SHOWLEVEL, SHOWLTYPE, SHOWFORMAT;
  int operatorID;
  int varID, zaxisID;
  int vdate, vtime;
  int nrecs, nvars, nout, ntsteps;
  int nlevs, levelID;
  int ltype;
  int index, nzaxis;
  int tsID, ndate, date0 = 0;
  int taxisID;
  int streamID;
  int vlistID;
  int year, month, day;
  int month0 = 0, nmonth, year0 = 0, nyear;
  char varname[256];
  char stdname[256];
  char vdatestr[32], vtimestr[32];

  cdoInitialize(argument);

  SHOWYEAR      = cdoOperatorAdd("showyear",      0, 0, NULL);
  SHOWMON       = cdoOperatorAdd("showmon",       0, 0, NULL);
  SHOWDATE      = cdoOperatorAdd("showdate",      0, 0, NULL);
  SHOWTIME      = cdoOperatorAdd("showtime",      0, 0, NULL);
  SHOWTIMESTAMP = cdoOperatorAdd("showtimestamp", 0, 0, NULL);
  SHOWCODE      = cdoOperatorAdd("showcode",      0, 0, NULL);
  SHOWPARAM     = cdoOperatorAdd("showparam",     0, 0, NULL);
  SHOWNAME      = cdoOperatorAdd("showname",      0, 0, NULL);
  SHOWSTDNAME   = cdoOperatorAdd("showstdname",   0, 0, NULL);
  SHOWLEVEL     = cdoOperatorAdd("showlevel",     0, 0, NULL);
  SHOWLTYPE     = cdoOperatorAdd("showltype",     0, 0, NULL);
  SHOWFORMAT    = cdoOperatorAdd("showformat",    0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID = streamOpenRead(cdoStreamName(0));
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));

  vlistID = streamInqVlist(streamID);

  nvars   = vlistNvars(vlistID);
  taxisID = vlistInqTaxis(vlistID);
  ntsteps = vlistNtsteps(vlistID);

  if ( operatorID == SHOWYEAR )
    {
      nyear = 0;
      tsID = 0;
      if ( ntsteps != 0 )
	while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	  {
	    vdate = taxisInqVdate(taxisID);

	    cdiDecodeDate(vdate, &year, &month, &day);
	 
	    if ( tsID == 0 || year0 != year )
	      {
		/* if ( nyear == 10 ) { nyear = 0; fprintf(stdout, "\n"); } */
		year0 = year;
		fprintf(stdout, " %4d", year0);
		nyear++;
	      }

	    tsID++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWMON )
    {
      nmonth = 0;
      tsID = 0;
      if ( ntsteps != 0 )
	while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	  {
	    vdate = taxisInqVdate(taxisID);

	    cdiDecodeDate(vdate, &year, &month, &day);
	 
	    if ( tsID == 0 || month0 != month )
	      {
		/* if ( nmonth == 12 ) { nmonth = 0; fprintf(stdout, "\n"); } */
		month0 = month;
		fprintf(stdout, " %2d", month0);
		nmonth++;
	      }

	    tsID++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWDATE )
    {
      ndate = 0;
      tsID  = 0;
      if ( ntsteps != 0 )
	while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	  {
	    vdate = taxisInqVdate(taxisID);
	 
	    date2str(vdate, vdatestr, sizeof(vdatestr));

	    if ( tsID == 0 || date0 != vdate )
	      {
		/* if ( ndate == 10 ) { ndate = 0; fprintf(stdout, "\n"); } */
		date0 = vdate;
		fprintf(stdout, " %s", vdatestr);
		ndate++;
	      }

	    tsID++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWTIME )
    {
      nout = 0;
      tsID = 0;
      if ( ntsteps != 0 )
	while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	  {
	    /* if ( nout == 4 ) { nout = 0; fprintf(stdout, "\n"); } */
	    vtime = taxisInqVtime(taxisID);

	    time2str(vtime, vtimestr, sizeof(vtimestr));
	    fprintf(stdout, " %s", vtimestr);

	    tsID++;
	    nout++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWTIMESTAMP )
    {
      nout = 0;
      tsID = 0;
      if ( ntsteps != 0 )
	while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	  {
	    /* if ( nout == 4 ) { nout = 0; fprintf(stdout, "\n"); } */
	    vdate = taxisInqVdate(taxisID);
	    vtime = taxisInqVtime(taxisID);

	    date2str(vdate, vdatestr, sizeof(vdatestr));
	    time2str(vtime, vtimestr, sizeof(vtimestr));
	    fprintf(stdout, " %sT%s", vdatestr, vtimestr);

	    tsID++;
	    nout++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWCODE )
    {
      nout = 0;
      for ( varID = 0; varID < nvars; varID++ )
	{
	  /* if ( nout == 20 ) { nout = 0; fprintf(stdout, "\n"); } */
	  fprintf(stdout, " %d", vlistInqVarCode(vlistID, varID));
	  nout++;
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWPARAM )
    {
      int param;
      char paramstr[32];
      
      for ( varID = 0; varID < nvars; varID++ )
	{
	  param   = vlistInqVarParam(vlistID, varID);
	  cdiParamToString(param, paramstr, sizeof(paramstr));

	  fprintf(stdout, " %s", paramstr);
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWNAME )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  vlistInqVarName(vlistID, varID, varname);
	  fprintf(stdout, " %s", varname);
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWSTDNAME )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  vlistInqVarStdname(vlistID, varID, stdname);
	  if ( stdname[0] != 0 )
	    fprintf(stdout, " %s", stdname);
	  else
	    fprintf(stdout, " unknown");
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWLEVEL )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  zaxisID = vlistInqVarZaxis(vlistID, varID);
	  nlevs = zaxisInqSize(zaxisID);
	  for ( levelID = 0; levelID < nlevs; levelID++ )
	    fprintf(stdout, " %.9g", zaxisInqLevel(zaxisID, levelID));
	  fprintf(stdout, "\n");
	}
    }
  else if ( operatorID == SHOWLTYPE )
    {
      nzaxis = vlistNzaxis(vlistID);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID = vlistZaxis(vlistID, index);

	  ltype = zaxis2ltype(zaxisID);

	  if ( ltype != -1 ) fprintf(stdout, " %d", ltype);
	}
      fprintf(stdout, "\n"); 
    }
  else if ( operatorID == SHOWFORMAT )
    {
      printFiletype(streamID, vlistID);
    }

  streamClose(streamID);

  cdoFinish();

  return (0);
}
