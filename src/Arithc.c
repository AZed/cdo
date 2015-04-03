/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

      Arithc     addc            Add by constant
      Arithc     subc            Subtract by constant
      Arithc     mulc            Multiply by constant
      Arithc     divc            Divide by constant
      Arithc     mod             Modulo operator
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


int *fill_vars(int vlistID)
{
  int varID;
  int nvars = vlistNvars(vlistID);
  int *vars = (int*) malloc(nvars*sizeof(int));

  if ( cdoNumVarnames )
    {
      int nfound = 0;
      char varname[CDI_MAX_NAME];

      for ( varID = 0; varID < nvars; ++varID )
	{
	  vars[varID] = 0;

	  vlistInqVarName(vlistID, varID, varname);

	  for ( int i = 0; i < cdoNumVarnames; ++i )
	    if ( strcmp(varname, cdoVarnames[i]) == 0 )
	      {
		vars[varID] = 1;
		nfound++;
		break;
	      }
	}

      if ( nfound == 0 )
	cdoAbort("Variable -n %s%s not found!", cdoVarnames[0], cdoNumVarnames > 1 ? ",..." : "");
    }
  else
    {
      for ( varID = 0; varID < nvars; ++varID ) vars[varID] = 1;
    }

  return (vars);
}


void *Arithc(void *argument)
{
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int gridsize, i;
  int nrecs, recID;
  int tsID;
  int varID, levelID;
  int vlistID1, vlistID2;
  double rconst;
  field_t field;
  int taxisID1, taxisID2;
  int *vars = NULL;

  cdoInitialize(argument);

  cdoOperatorAdd("addc", func_add, 0, "constant value");
  cdoOperatorAdd("subc", func_sub, 0, "constant value");
  cdoOperatorAdd("mulc", func_mul, 0, "constant value");
  cdoOperatorAdd("divc", func_div, 0, "constant value");
  cdoOperatorAdd("mod",  func_mod, 0, "divisor");

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);

  operatorInputArg(cdoOperatorEnter(operatorID));
  rconst = atof(operatorArgv()[0]);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  vars = fill_vars(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);

  field_init(&field);
  field.ptr    = (double*) malloc(gridsize*sizeof(double));
  field.weight = NULL;

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field.ptr, &field.nmiss);

	  if ( vars[varID] )
	    {
	      field.grid    = vlistInqVarGrid(vlistID1, varID);
	      field.missval = vlistInqVarMissval(vlistID1, varID);

	      farcfun(&field, rconst, operfunc);

	      /* recalculate number of missing values */
	      gridsize = gridInqSize(field.grid);
	      field.nmiss = 0;
	      for ( i = 0; i < gridsize; ++i )
		if ( DBL_IS_EQUAL(field.ptr[i], field.missval) ) field.nmiss++;
	    }

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, field.ptr, field.nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field.ptr ) free(field.ptr);
  if ( vars ) free(vars);

  vlistDestroy(vlistID2);

  cdoFinish();

  return (0);
}
