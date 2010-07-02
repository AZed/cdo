#include <string.h>

#include "cdi.h"
#include "dmemory.h"
#include "field.h"


field_t **field_allocate(int vlistID, int ptype, int init)
{
  static const char *func = "field_allocate";
  int nvars, nlevel;
  int varID, zaxisID, levelID;
  int gridID, gridsize;
  double missval;
  field_t **field;

  nvars = vlistNvars(vlistID);

  field = (field_t **) malloc(nvars*sizeof(field_t *));

  for ( varID = 0; varID < nvars; ++varID )
    {
      gridID   = vlistInqVarGrid(vlistID, varID);
      gridsize = gridInqSize(gridID);
      zaxisID  = vlistInqVarZaxis(vlistID, varID);
      nlevel   = zaxisInqSize(zaxisID);
      missval  = vlistInqVarMissval(vlistID, varID);

      field[varID] = (field_t *)  malloc(nlevel*sizeof(field_t));
      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  field[varID][levelID].grid    = gridID;
	  field[varID][levelID].nsamp   = 0;
	  field[varID][levelID].nmiss   = 0;
	  field[varID][levelID].missval = missval;
	  field[varID][levelID].ptr     = NULL;
	  field[varID][levelID].weight  = NULL;

	  if ( ptype == FIELD_ALL || ptype == FIELD_PTR )
	    {
	      field[varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));
	      if ( init ) memset(field[varID][levelID].ptr, 0, gridsize*sizeof(double));
	    }

	  if ( ptype == FIELD_ALL || ptype == FIELD_PTR )
	    {
	      field[varID][levelID].weight = (double *) malloc(gridsize*sizeof(double));
	      if ( init ) memset(field[varID][levelID].weight, 0, gridsize*sizeof(double));
	    }    
	}
    }

  return (field);
}


field_t **field_malloc(int vlistID, int ptype)
{
  return (field_allocate(vlistID, ptype, 0));
}


field_t **field_calloc(int vlistID, int ptype)
{
  return (field_allocate(vlistID, ptype, 1));
}


void field_free(field_t **field, int vlistID)
{
  static const char *func = "field_free";
  int nvars, nlevel;
  int varID, levelID;

  nvars = vlistNvars(vlistID);

  for ( varID = 0; varID < nvars; ++varID )
    {
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  if ( field[varID][levelID].ptr )    free(field[varID][levelID].ptr);
       	  if ( field[varID][levelID].weight ) free(field[varID][levelID].weight);
	}

      free(field[varID]);
    }

  free(field);
}
