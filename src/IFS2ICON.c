/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2011 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "namelist.h"

#define  MAX_PARAM  512

typedef struct {
  char *name;
  int   pos;
  int   itype;
} param_t;

typedef struct {
  char nml_name[32];
  int debug_level;
  int force_remap;
} nml_control_t;


static
void nml_error(int status, const char *nml_name)
{
  switch (status)
    {
    case 1: 
      {
	cdoAbort("Namelist %s not found!", nml_name);
	break;
      }
    default:
      {
	cdoAbort("Namelist %s unknown error!", nml_name);
      }
    }
}

static
int pos_nml(FILE *fp_nml, const char *nml_name)
{
  int status = 1;
  fpos_t fpos;
  char line[MAX_LINE_LEN];
  

  rewind(fp_nml);
  
  fgetpos(fp_nml, &fpos);
  while ( readline(fp_nml, line, MAX_LINE_LEN) )
    {
      if ( line[0] == '#' ) continue;
      if ( line[0] == '\0' ) continue;
      if ( line[0] == '&' && strlen(line) > 2 )
	{
	  strtolower(line+1);
	  if ( strncmp(line+1, nml_name, strlen(nml_name)) == 0 )
	    {
	      status = 0;
	      fsetpos(fp_nml, &fpos);
	      break;
	    }
	}
      fgetpos(fp_nml, &fpos);
    }

  return status;
}

static
void nml_control_init(nml_control_t *nml_control)
{
  strcpy(nml_control->nml_name, "control");
  nml_control->debug_level = 0;
  nml_control->force_remap = 0;
}

static
void nml_control_print(nml_control_t nml_control)
{ 
  fprintf(stdout, "Namelist: %s\n", nml_control.nml_name);
  fprintf(stdout, "  debug_level = %d\n", nml_control.debug_level);
  fprintf(stdout, "  force_remap = %d\n", nml_control.force_remap);
}

static
int nml_control_read(FILE *fp_nml, nml_control_t *nml_control)
{
  int status = 0;

  status = pos_nml(fp_nml, nml_control->nml_name);

  if ( status == 0 )
    {
      namelist_t *nml;
      NML_DEF_INT(debug_level, 1, 0);
      NML_DEF_INT(force_remap, 1, 0);

      nml = namelistNew(nml_control->nml_name);

      NML_ADD_INT(nml, debug_level);
      NML_ADD_INT(nml, force_remap);
      
      namelistRead(fp_nml, nml);

      if ( (NML_NUM(nml, debug_level)) ) nml_control->debug_level = seldebug_level[0];
      if ( (NML_NUM(nml, force_remap)) ) nml_control->force_remap = selforce_remap[0];

      namelistDelete(nml);
    }

  return status;
}

static
void read_param(const char *paramfile, param_t *param, int maxparam, int *nparam)
{
  FILE *fp;
  namelist_t *nml;
  int nml_itype, nml_name, nml_pos;
  int locc, i;
  int code, new_code, table, pos;
  int nml_index = 0;
  int codenum, tabnum, levtype;
  int numparam = 0;
  char *datatype = NULL;
  char *name = NULL, *itype = NULL, *stdname = NULL, longname[CDI_MAX_NAME] = "", units[CDI_MAX_NAME] = "";
  char varname[CDI_MAX_NAME];

  fp = fopen(paramfile, "r");
  if ( fp == NULL ) cdoAbort("Open failed on %s!", paramfile);

  nml = namelistNew("parameter");
  nml->dis = 0;

  nml_name     = namelistAdd(nml, "name",          NML_WORD, 0, &name, 1);
  nml_pos      = namelistAdd(nml, "pos",           NML_INT,  0, &pos, 1);
  nml_itype    = namelistAdd(nml, "itype",         NML_WORD, 0, &itype, 1);
	      
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

	  if ( nml->entry[nml_name]->occ == 0 )
	    {
	      cdoWarning("Parameter %d skipped, variable name not found!", nml_index);
	      continue;
	    }

	  if ( numparam >= maxparam )
	    {
	      cdoWarning("Too many parameter (limit=%d)!", maxparam);
	      break;
	    }

	  param[numparam].name = strdup(name);

	  if ( nml->entry[nml_pos]->occ )
	    param[numparam].pos  = pos;
	  else
	    param[numparam].pos  = 0;

	  if ( nml->entry[nml_itype]->occ )
	    param[numparam].itype = (int) *itype;
	  else
	    param[numparam].itype = 'Q';
	
	  numparam++;
	  /*
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      vlistInqVarName(vlistID2, varID, varname);
	      if ( strcmp(varname, name) == 0 ) break;
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
		  cdoPrint("Variable %s not found!", name);
		}
	    }
	  */
	}
      else
	break;
    }
  
  namelistDelete(nml);

  fclose(fp);

  *nparam = numparam;
}


void *IFS2ICON(void *argument)
{
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
  int i;
  int status;
  int nparam;
  const char *paramfile;
  const char *nmlfile;
  double newlevel = 0;
  double *levels = NULL;
  double *array = NULL;
  param_t param[MAX_PARAM];
  nml_control_t nml_control;
  FILE *fp_nml;

  cdoInitialize(argument);

  operatorInputArg("namelist file and parameter file");
  operatorCheckArgc(2);

  /* read namelists */
  nmlfile = operatorArgv()[0];

  fp_nml = fopen(nmlfile, "r");
  if ( fp_nml == NULL ) cdoAbort("Open failed on %s", nmlfile);

  nml_control_init(&nml_control);
  status = nml_control_read(fp_nml, &nml_control);
  if ( status > 0 ) nml_error(status, nml_control.nml_name);

  fclose(fp_nml);

  if ( cdoVerbose )
    {
      nml_control_print(nml_control);
    }
  
  paramfile = operatorArgv()[1];
  read_param(paramfile, param, MAX_PARAM, &nparam);

  for ( i = 0; i < nparam; ++i )
    printf("%s %d %c\n", param[i].name, param[i].pos, param[i].itype);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);
  vlistPrint(vlistID2);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID2);
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
