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

      Set        setpartab       Set parameter table
      Set        setcode         Set code number
      Set        setparam        Set parameter identifier
      Set        setname         Set variable name
      Set        setlevel        Set level
      Set        setltype        Set GRIB level type
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "namelist.h"


int stringToParam(const char *paramstr)
{
  int param = 0;
  int pnum = -1, pcat = 255, pdis = 255;
  size_t len;

  len = strlen(paramstr);

  sscanf(paramstr, "%d.%d.%d", &pnum, &pcat, &pdis);
  
  if ( cdoVerbose ) cdoPrint("pnum, pcat, pdis: %d.%d.%d", pnum, pcat, pdis);

  param = cdiEncodeParam(pnum, pcat, pdis);

  return (param);
}


void *Set(void *argument)
{
  int SETPARTAB, SETPARTABV, SETCODE, SETPARAM, SETNAME, SETLEVEL, SETLTYPE, SETTABNUM;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs, nvars, newval = -1, tabnum = 0;
  int tsID1, recID, varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int nmiss;
  int gridsize;
  int index, zaxisID1, zaxisID2, nzaxis, nlevs;
  int tableID = -1;
  int tableformat = 0;
  int zaxistype;
  int newparam = 0;
  char *newname = NULL, *partab = NULL;
  double newlevel = 0;
  double *levels = NULL;
  double *array = NULL;

  cdoInitialize(argument);

  SETPARTAB  = cdoOperatorAdd("setpartab",  0, 0, "parameter table");
  SETPARTABV = cdoOperatorAdd("setpartabv", 0, 0, "parameter table");
  SETCODE    = cdoOperatorAdd("setcode",    0, 0, "code number");
  SETPARAM   = cdoOperatorAdd("setparam",   0, 0, "parameter identifier (format: code[.tabnum] or num[.cat[.dis]])");
  SETNAME    = cdoOperatorAdd("setname",    0, 0, "variable name");
  SETLEVEL   = cdoOperatorAdd("setlevel",   0, 0, "level");
  SETLTYPE   = cdoOperatorAdd("setltype",   0, 0, "GRIB level type");
  SETTABNUM  = cdoOperatorAdd("settabnum",  0, 0, "GRIB table number");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));
  if ( operatorID == SETCODE || operatorID == SETLTYPE )
    {
      newval = atoi(operatorArgv()[0]);
    }
  else if ( operatorID == SETPARAM )
    {
      newparam = stringToParam(operatorArgv()[0]);
    }
  else if ( operatorID == SETNAME )
    {
      newname = operatorArgv()[0];
    }
  else if ( operatorID == SETTABNUM )
    {
      tabnum = atoi(operatorArgv()[0]);
    }
  else if ( operatorID == SETPARTAB )
    {
      FILE *fp;
      size_t fsize;
      char *parbuf = NULL;
      size_t nbytes;

      partab = operatorArgv()[0];
      fp = fopen(partab, "r");
      if ( fp != NULL )
	{
	  fseek(fp, 0L, SEEK_END);
	  fsize = (size_t) ftell(fp);
	  parbuf = (char *) malloc(fsize+1);
	  fseek(fp, 0L, SEEK_SET);
	  nbytes = fread(parbuf, fsize, 1, fp);
	  parbuf[fsize] = 0;
	  fseek(fp, 0L, SEEK_SET);

	  if ( atoi(parbuf) == 0 ) tableformat = 1;

	  fclose(fp);
	  free(parbuf);
	}

      if ( tableformat == 0 ) tableID = defineTable(partab);
    }
  else if ( operatorID == SETPARTABV )
    {
      tableformat = 1;
    }
  else if ( operatorID == SETLEVEL )
    {
      newlevel = atof(operatorArgv()[0]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);
  /* vlistPrint(vlistID2);*/

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operatorID == SETCODE )
    {
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarCode(vlistID2, varID, newval);
    }
  else if ( operatorID == SETPARAM )
    {
      vlistDefVarParam(vlistID2, 0, newparam);
    }
  else if ( operatorID == SETNAME )
    {
      vlistDefVarName(vlistID2, 0, newname);
    }
  else if ( operatorID == SETTABNUM )
    {
      int tableID;
      tableID = tableDef(-1, tabnum, NULL);
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarTable(vlistID2, varID, tableID);
    }
  else if ( operatorID == SETPARTAB || operatorID == SETPARTABV )
    {
      nvars = vlistNvars(vlistID2);

      if ( tableformat == 0 )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    vlistDefVarTable(vlistID2, varID, tableID);
	}
      else
	{
	  FILE *fp;
	  namelist_t *nml;
	  int nml_code, nml_new_code, nml_table, nml_datatype, nml_name, nml_new_name, nml_stdname;
	  int nml_longname, nml_units, nml_ltype;
	  int locc, i;
	  int code, new_code, table, ltype;
	  int nml_index = 0;
	  int codenum, tabnum, levtype;
	  char *datatype = NULL;
	  char *name = NULL, *new_name = NULL, *stdname = NULL, longname[CDI_MAX_NAME] = "", units[CDI_MAX_NAME] = "";
	  char varname[CDI_MAX_NAME];

	  partab = operatorArgv()[0];
	  fp = fopen(partab, "r");
	  if ( fp == NULL ) cdoAbort("Internal problem! Parameter table %s not available", partab);

	  nml = namelistNew("parameter");
	  nml->dis = 0;

	  nml_code     = namelistAdd(nml, "code",          NML_INT,  0, &code, 1);
	  nml_new_code = namelistAdd(nml, "new_code",      NML_INT,  0, &new_code, 1);
	  nml_table    = namelistAdd(nml, "table",         NML_INT,  0, &table, 1);
	  nml_ltype    = namelistAdd(nml, "ltype",         NML_INT,  0, &ltype, 1);
	  nml_datatype = namelistAdd(nml, "datatype",      NML_WORD, 0, &datatype, 1);
	  nml_name     = namelistAdd(nml, "name",          NML_WORD, 0, &name, 1);
	  nml_new_name = namelistAdd(nml, "new_name",      NML_WORD, 0, &new_name, 1);
	  nml_stdname  = namelistAdd(nml, "standard_name", NML_WORD, 0, &stdname, 1);
	  nml_longname = namelistAdd(nml, "long_name",     NML_TEXT, 0, longname, sizeof(longname));
	  nml_units    = namelistAdd(nml, "units",         NML_TEXT, 0, units, sizeof(units));
	      
	  while ( ! feof(fp) )
	    {
	      namelistReset(nml);

	      namelistRead(fp, nml);

	      locc = FALSE;
	      for ( i = 0; i < nml->size; i++ )
		{
		  if ( nml->entry[i]->occ ) { locc = TRUE; break; }
		}

	      if ( locc )
		{
		  /* namelistPrint(nml); */

		  nml_index++;

		  if ( operatorID == SETPARTAB )
		    {
		      if ( nml->entry[nml_code]->occ == 0 )
			{
			  cdoPrint("Parameter %d skipped, code number not found!", nml_index);
			  continue;
			}
		    }
		  else
		    {
		      if ( nml->entry[nml_name]->occ == 0 )
			{
			  cdoWarning("Parameter %d skipped, variable name not found!", nml_index);
			  continue;
			}
		    }

		  for ( varID = 0; varID < nvars; varID++ )
		    {
		      if ( operatorID == SETPARTAB )
			{
			  codenum = vlistInqVarCode(vlistID2, varID);
			  tableID = vlistInqVarTable(vlistID2, varID);
			  tabnum  = tableInqNum(tableID);
			  levtype = zaxisInqLtype(vlistInqVarZaxis(vlistID2, varID));
			  /*
			  printf("code = %d  tabnum = %d  ltype = %d\n", codenum, tabnum, levtype);
			  */
			  if ( nml->entry[nml_table]->occ == 0 ) table = tabnum;
			  if ( nml->entry[nml_ltype]->occ == 0 ) ltype = levtype;

			  if ( codenum == code && tabnum == table && levtype == ltype ) break;
			}
		      else
			{
			  vlistInqVarName(vlistID2, varID, varname);
			  if ( strcmp(varname, name) == 0 ) break;
			}
		    }

		  if ( varID < nvars )
		    {
		      if ( nml->entry[nml_code]->occ     ) vlistDefVarCode(vlistID2, varID, code);
		      if ( nml->entry[nml_new_code]->occ ) vlistDefVarCode(vlistID2, varID, new_code);
		      if ( nml->entry[nml_name]->occ     ) vlistDefVarName(vlistID2, varID, name);
		      if ( nml->entry[nml_new_name]->occ ) vlistDefVarName(vlistID2, varID, new_name);
		      if ( nml->entry[nml_stdname]->occ  ) vlistDefVarStdname(vlistID2, varID, stdname);
		      if ( nml->entry[nml_longname]->occ ) vlistDefVarLongname(vlistID2, varID, longname);
		      if ( nml->entry[nml_units]->occ    ) vlistDefVarUnits(vlistID2, varID, units);
		    }
		  else
		    {
		      if ( cdoVerbose )
			{
			  if ( operatorID == SETPARTAB )
			    {
			      if ( nml->entry[nml_table]->occ == 0 )
				cdoPrint("Code %d not found!", code);
			      else
				cdoPrint("Code %d and table %d not found!", code, table);
			    }
			  else
			    cdoPrint("Variable %s not found!", name);
			}
		    }
		}
	      else
		break;
	    }
	  
	  namelistDelete(nml);

	  fclose(fp);
	}
    }
  else if ( operatorID == SETLEVEL )
    {
      nzaxis = vlistNzaxis(vlistID2);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID2, index);
	  zaxisID2 = zaxisDuplicate(zaxisID1);
	  nlevs = zaxisInqSize(zaxisID2);
	  levels = (double *) malloc(nlevs*sizeof(double));
	  zaxisInqLevels(zaxisID2, levels);
	  levels[0] = newlevel;
	  zaxisDefLevels(zaxisID2, levels);
	  vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
	  free(levels);
	}
    }
  else if ( operatorID == SETLTYPE )
    {
      nzaxis = vlistNzaxis(vlistID2);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID2, index);
	  zaxisID2 = zaxisDuplicate(zaxisID1);

	  zaxistype = ltype2ztype(newval);

	  zaxisChangeType(zaxisID2, zaxistype);
	  if ( zaxistype == ZAXIS_GENERIC ) zaxisDefLtype(zaxisID2, newval);
	  vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
	}
    }

  /* vlistPrint(vlistID2);*/
  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  array = (double *) malloc(gridsize*sizeof(double));

  tsID1 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID1);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  streamReadRecord(streamID1, array, &nmiss);
	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID1++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}
