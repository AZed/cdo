/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Sinfo      sinfo           Short dataset information
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#include "printinfo.h"

#define MAXCHARS 82

void *Sinfo(void *argument)
{
  int SINFO, SINFOP, SINFON, SINFOC;
  int operatorID;
  int indf;
  int varID;
  int gridsize = 0;
  int gridID, zaxisID, code, tabnum, param;
  int zaxistype, ltype;
  int vdate, vtime;
  int nrecs, nvars, nzaxis, ntsteps;
  int levelID, levelsize;
  int tsID, ntimeout;
  int tsteptype, taxisID;
  int nbyte, nbyte0;
  int index;
  char varname[CDI_MAX_NAME];
  char longname[CDI_MAX_NAME];
  char units[CDI_MAX_NAME];
  char paramstr[32];
  char vdatestr[32], vtimestr[32];
  double level;
  char *modelptr, *instptr;
  int streamID = 0;
  int vlistID;
  int datatype;
  char pstr[4];

  cdoInitialize(argument);

  SINFO  = cdoOperatorAdd("sinfo",  0, 0, NULL);
  SINFOP = cdoOperatorAdd("sinfop", 0, 0, NULL);
  SINFON = cdoOperatorAdd("sinfon", 0, 0, NULL);
  SINFOC = cdoOperatorAdd("sinfoc", 0, 0, NULL);

  operatorID = cdoOperatorID();

  for ( indf = 0; indf < cdoStreamCnt(); indf++ )
    {
      streamID = streamOpenRead(cdoStreamName(indf));

      vlistID = streamInqVlist(streamID);

      printf("   File format: ");
      printFiletype(streamID, vlistID);

      if ( operatorID == SINFON )
	fprintf(stdout,
		"%6d : Institut Source   Name        Ttype   Dtype  Gridsize Num  Levels Num\n",  -(indf+1));
      else if ( operatorID == SINFOC )
	fprintf(stdout,
		"%6d : Institut Source  Table Code   Ttype   Dtype  Gridsize Num  Levels Num\n",  -(indf+1));
      else
	fprintf(stdout,
		"%6d : Institut Source   Param       Ttype   Dtype  Gridsize Num  Levels Num\n",  -(indf+1));

      nvars = vlistNvars(vlistID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  param   = vlistInqVarParam(vlistID, varID);
	  code    = vlistInqVarCode(vlistID, varID);
	  tabnum  = tableInqNum(vlistInqVarTable(vlistID, varID));
	  gridID  = vlistInqVarGrid(vlistID, varID);
	  zaxisID = vlistInqVarZaxis(vlistID, varID);

	  cdiParamToString(param, paramstr, sizeof(paramstr));

	  if ( operatorID == SINFON ) vlistInqVarName(vlistID, varID, varname);

	  gridsize = gridInqSize(gridID);

	  fprintf(stdout, "%6d : ", varID + 1);

	  instptr = institutInqNamePtr(vlistInqVarInstitut(vlistID, varID));
	  if ( instptr )
	    fprintf(stdout, "%-8s ", instptr);
	  else
	    fprintf(stdout, "unknown  ");

	  modelptr = modelInqNamePtr(vlistInqVarModel(vlistID, varID));
	  if ( modelptr )
	    fprintf(stdout, "%-8s ", modelptr);
	  else
	    fprintf(stdout, "unknown  ");

	  if ( operatorID == SINFON )
	    fprintf(stdout, "%-11s ", varname);
	  else if ( operatorID == SINFOC )
	    fprintf(stdout, "%4d %4d   ", tabnum, code);
	  else
	    fprintf(stdout, "%-11s ", paramstr);

	  tsteptype = vlistInqVarTsteptype(vlistID, varID);
	  if      ( tsteptype == TSTEP_CONSTANT ) fprintf(stdout, "%-8s", "constant");
	  else if ( tsteptype == TSTEP_INSTANT  ) fprintf(stdout, "%-8s", "instant");
	  else if ( tsteptype == TSTEP_MIN      ) fprintf(stdout, "%-8s", "min");
	  else if ( tsteptype == TSTEP_MAX      ) fprintf(stdout, "%-8s", "max");
	  else if ( tsteptype == TSTEP_ACCUM    ) fprintf(stdout, "%-8s", "accum");
	  else                                    fprintf(stdout, "%-8s", "unknown");

	  datatype = vlistInqVarDatatype(vlistID, varID);

	  if      ( datatype == DATATYPE_PACK   ) strcpy(pstr, "P0");
	  else if ( datatype > 0 && datatype <= 32  ) sprintf(pstr, "P%d", datatype);
	  else if ( datatype == DATATYPE_CPX32  ) strcpy(pstr, "C32");
	  else if ( datatype == DATATYPE_CPX64  ) strcpy(pstr, "C64");
	  else if ( datatype == DATATYPE_FLT32  ) strcpy(pstr, "F32");
	  else if ( datatype == DATATYPE_FLT64  ) strcpy(pstr, "F64");
	  else if ( datatype == DATATYPE_INT8   ) strcpy(pstr, "I8");
	  else if ( datatype == DATATYPE_INT16  ) strcpy(pstr, "I16");
	  else if ( datatype == DATATYPE_INT32  ) strcpy(pstr, "I32");
	  else if ( datatype == DATATYPE_UINT8  ) strcpy(pstr, "U8");
	  else if ( datatype == DATATYPE_UINT16 ) strcpy(pstr, "U16");
	  else if ( datatype == DATATYPE_UINT32 ) strcpy(pstr, "U32");
	  else                                    strcpy(pstr, "-1");

	  fprintf(stdout, " %-3s", pstr);

	  if ( vlistInqVarCompType(vlistID, varID) == COMPRESS_NONE )
	    fprintf(stdout, " ");
	  else
	    fprintf(stdout, "z");

	  fprintf(stdout, "%9d", gridsize);

	  fprintf(stdout, " %3d ", vlistGridIndex(vlistID, gridID) + 1);

	  levelsize = zaxisInqSize(zaxisID);
	  fprintf(stdout, " %6d", levelsize);
	  fprintf(stdout, " %3d", vlistZaxisIndex(vlistID, zaxisID) + 1);

	  fprintf(stdout, "\n");
	}

      fprintf(stdout, "   Horizontal grids :\n");
      printGridInfo(vlistID);

      nzaxis = vlistNzaxis(vlistID);
      fprintf(stdout, "   Vertical grids :\n");
      for ( index = 0; index < nzaxis; index++)
	{
	  zaxisID   = vlistZaxis(vlistID, index);
	  zaxistype = zaxisInqType(zaxisID);
	  ltype     = zaxisInqLtype(zaxisID);
	  levelsize = zaxisInqSize(zaxisID);
	  /* zaxisInqLongname(zaxisID, longname); */
	  zaxisName(zaxistype, longname);
	  longname[17] = 0;
	  zaxisInqUnits(zaxisID, units);
	  units[12] = 0;
	  if ( zaxistype == ZAXIS_GENERIC && ltype != 0 )
	    nbyte0    = fprintf(stdout, "  %4d : %-11s  (ltype=%3d) : ", vlistZaxisIndex(vlistID, zaxisID)+1, longname, ltype);
	  else
	    nbyte0    = fprintf(stdout, "  %4d : %-17s  %5s : ", vlistZaxisIndex(vlistID, zaxisID)+1, longname, units);
	  nbyte = nbyte0;
	  for ( levelID = 0; levelID < levelsize; levelID++ )
	    {
	      if ( nbyte > MAXCHARS )
		{
		  fprintf(stdout, "\n");
		  fprintf(stdout, "%*s", nbyte0, "");
		  nbyte = nbyte0;
		}
	      level = zaxisInqLevel(zaxisID, levelID);
	      nbyte += fprintf(stdout, "%.9g ", level);
	    }
	  fprintf(stdout, "\n");
	  if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	    {
	      double level1, level2;
	      nbyte = nbyte0;
	      nbyte0 = fprintf(stdout, "%33s : ", "bounds");
	      for ( levelID = 0; levelID < levelsize; levelID++ )
		{
		  if ( nbyte > MAXCHARS )
		    {
		      fprintf(stdout, "\n");
		      fprintf(stdout, "%*s", nbyte0, "");
		      nbyte = nbyte0;
		    }
		  level1 = zaxisInqLbound(zaxisID, levelID);
		  level2 = zaxisInqUbound(zaxisID, levelID);
		  nbyte += fprintf(stdout, "%.9g-%.9g ", level1, level2);
		}
	      fprintf(stdout, "\n");
	    }
	}

      taxisID = vlistInqTaxis(vlistID);
      ntsteps = vlistNtsteps(vlistID);

      if ( ntsteps != 0 )
	{
	  if ( ntsteps == CDI_UNDEFID )
	    fprintf(stdout, "   Time axis :  unlimited steps\n");
	  else
	    fprintf(stdout, "   Time axis :  %d step%s\n", ntsteps, ntsteps == 1 ? "" : "s");

	  if ( taxisID != CDI_UNDEFID )
	    {
	      int calendar, tunits;

	      if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
		{
		  vdate = taxisInqRdate(taxisID);
		  vtime = taxisInqRtime(taxisID);

		  date2str(vdate, vdatestr, sizeof(vdatestr));
		  time2str(vtime, vtimestr, sizeof(vtimestr));

		  fprintf(stdout, "     RefTime = %s %s", vdatestr, vtimestr);
		      
		  tunits = taxisInqTunit(taxisID);
		  if ( tunits != CDI_UNDEFID )
		    {
		      if ( tunits == TUNIT_YEAR )
			fprintf(stdout, "  Units = years");
		      else if ( tunits == TUNIT_MONTH )
			fprintf(stdout, "  Units = months");
		      else if ( tunits == TUNIT_DAY )
			fprintf(stdout, "  Units = days");
		      else if ( tunits == TUNIT_12HOURS )
			fprintf(stdout, "  Units = 12hours");
		      else if ( tunits == TUNIT_6HOURS )
			fprintf(stdout, "  Units = 6hours");
		      else if ( tunits == TUNIT_3HOURS )
			fprintf(stdout, "  Units = 3hours");
		      else if ( tunits == TUNIT_HOUR )
			fprintf(stdout, "  Units = hours");
		      else if ( tunits == TUNIT_MINUTE )
			fprintf(stdout, "  Units = minutes");
		      else if ( tunits == TUNIT_SECOND )
			fprintf(stdout, "  Units = seconds");
		      else
			fprintf(stdout, "  Units = unknown");
		    }
	      
		  calendar = taxisInqCalendar(taxisID);
		  if ( calendar != CDI_UNDEFID )
		    {
		      if      ( calendar == CALENDAR_STANDARD )
			fprintf(stdout, "  Calendar = STANDARD");
		      else if ( calendar == CALENDAR_PROLEPTIC )
			fprintf(stdout, "  Calendar = PROLEPTIC");
		      else if ( calendar == CALENDAR_360DAYS )
			fprintf(stdout, "  Calendar = 360DAYS");
		      else if ( calendar == CALENDAR_365DAYS )
			fprintf(stdout, "  Calendar = 365DAYS");
		      else if ( calendar == CALENDAR_366DAYS )
			fprintf(stdout, "  Calendar = 366DAYS");
		      else
			fprintf(stdout, "  Calendar = unknown");
		    }

		  if ( taxisHasBounds(taxisID) )
		    fprintf(stdout, "  Bounds = true");

		  fprintf(stdout, "\n");
		}
	    }

	  fprintf(stdout, "  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss\n");

	  ntimeout = 0;
	  tsID = 0;
	  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	    {
	      if ( ntimeout == 4 )
		{
		  ntimeout = 0;
		  fprintf(stdout, "\n");
		}

	      vdate = taxisInqVdate(taxisID);
	      vtime = taxisInqVtime(taxisID);

	      date2str(vdate, vdatestr, sizeof(vdatestr));
	      time2str(vtime, vtimestr, sizeof(vtimestr));

	      fprintf(stdout, " %s %s", vdatestr, vtimestr);

	      ntimeout++;
	      tsID++;
	    }
	  fprintf(stdout, "\n");
	}

      streamClose(streamID);
    }

  cdoFinish();

  return (0);
}
