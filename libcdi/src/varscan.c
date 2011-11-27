#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <string.h>
#include <math.h>

#include "cdi.h"
#include "stream_int.h"
#include "dmemory.h"
#include "varscan.h"
#include "vlist.h"

#undef  UNDEFID
#define UNDEFID -1

static size_t Vctsize = 0;
static double *Vct = NULL;

typedef struct
{
  int      level1;
  int      level2;
  int      recID;
  int      lindex;
}
leveltable_t;

typedef struct
{
  int         param;
  int         prec;
  int         tsteptype;
  int         timave;
  int         timaccu;
  int         gridID;
  int         zaxistype;
  int         ltype;     /* GRIB level type */
  int         lbounds;
  int         zaxisID;
  int         nlevels;
  int         levelTableSize;
  leveltable_t *levelTable;
  int         instID;
  int         modelID;
  int         tableID;
  int         comptype;       // compression type
  int         complevel;      // compression level
  int         lmissval;
  double      missval;
  char       *name;
  char       *longname;
  char       *units;
}
vartable_t;


int vartableInit = 0;
vartable_t *vartable;
static int varTablesize = 0;
int nvars = 0;


static
void paramInitEntry(int varID, int param)
{
  vartable[varID].param          = param;
  vartable[varID].prec           = 0;
  vartable[varID].tsteptype      = TSTEP_INSTANT;
  vartable[varID].timave         = 0;
  vartable[varID].timaccu        = 0;
  vartable[varID].gridID         = UNDEFID;
  vartable[varID].zaxistype      = 0;
  vartable[varID].ltype          = 0;
  vartable[varID].levelTable     = NULL;
  vartable[varID].levelTableSize = 0;
  vartable[varID].nlevels        = 0;
  vartable[varID].instID         = UNDEFID;
  vartable[varID].modelID        = UNDEFID;
  vartable[varID].tableID        = UNDEFID;
  vartable[varID].comptype       = COMPRESS_NONE;
  vartable[varID].complevel      = 1;
  vartable[varID].lmissval       = 0;
  vartable[varID].missval        = 0;
  vartable[varID].name           = NULL;
  vartable[varID].longname       = NULL;
  vartable[varID].units          = NULL;
}


static
int varGetEntry(int param, int zaxistype, int ltype)
{
  int varID;

  for ( varID = 0; varID < varTablesize; varID++ )
    {
      if ( vartable[varID].param     == param     &&
	   vartable[varID].zaxistype == zaxistype &&
	   vartable[varID].ltype     == ltype )
	return (varID);
    }

  return (UNDEFID);
}


void varFree(void)
{
  int varID;

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vartable[varID].levelTable )
	free(vartable[varID].levelTable);

      if ( vartable[varID].name )     free(vartable[varID].name);
      if ( vartable[varID].longname ) free(vartable[varID].longname);
      if ( vartable[varID].units )    free(vartable[varID].units);
    }

  if ( vartable )
    free(vartable);

  vartable = NULL;
  varTablesize = 0;
  nvars = 0;

  if ( Vct )
    free(Vct);

  Vct = NULL;
  Vctsize = 0;
}

int levelNewEntry(int varID, int level1, int level2)
{
  int levelID = 0;
  int levelTableSize;
  leveltable_t *levelTable;

  levelTableSize = vartable[varID].levelTableSize;
  levelTable     = vartable[varID].levelTable;

  /*
    Look for a free slot in levelTable.
    (Create the table the first time through).
  */
  if ( ! levelTableSize )
    {
      int i;

      levelTableSize = 2;
      levelTable = (leveltable_t *) malloc(levelTableSize*sizeof(leveltable_t));
      if( levelTable == NULL )
	{
          Message("levelTableSize = %d", levelTableSize);
	  SysError("Allocation of leveltable failed!");
	}

      for( i = 0; i < levelTableSize; i++ )
	levelTable[i].recID = UNDEFID;
    }
  else
    {
      while( levelID < levelTableSize )
	{
	  if ( levelTable[levelID].recID == UNDEFID ) break;
	  levelID++;
	}
    }
  /*
    If the table overflows, double its size.
  */
  if( levelID == levelTableSize )
    {
      int i;

      levelTableSize = 2*levelTableSize;
      levelTable = (leveltable_t *) realloc(levelTable, levelTableSize*sizeof(leveltable_t));
      if( levelTable == NULL )
	{
          Message("levelTableSize = %d", levelTableSize);
	  SysError("Reallocation of leveltable failed");
	}
      levelID = levelTableSize/2;

      for( i = levelID; i < levelTableSize; i++ )
	levelTable[i].recID = UNDEFID;
    }

  levelTable[levelID].level1 = level1;
  levelTable[levelID].level2 = level2;
  levelTable[levelID].lindex = levelID;

  vartable[varID].nlevels = levelID+1;
  vartable[varID].levelTableSize = levelTableSize;
  vartable[varID].levelTable = levelTable;

  return (levelID);
}

#define  UNDEF_PARAM  -4711

int paramNewEntry (int param)
{
  int varID = 0;

  /*
    Look for a free slot in vartable.
    (Create the table the first time through).
  */
  if ( ! varTablesize )
    {
      int i;

      varTablesize = 2;
      vartable = (vartable_t *) malloc(varTablesize*sizeof(vartable_t));
      if( vartable == NULL )
	{
          Message("varTablesize = %d", varTablesize);
	  SysError("Allocation of vartable failed");
	}

      for( i = 0; i < varTablesize; i++ )
	vartable[i].param = UNDEF_PARAM;
    }
  else
    {
      while( varID < varTablesize )
	{
	  if ( vartable[varID].param == UNDEF_PARAM ) break;
	  varID++;
	}
    }
  /*
    If the table overflows, double its size.
  */
  if ( varID == varTablesize )
    {
      int i;

      varTablesize = 2*varTablesize;
      vartable = (vartable_t *) realloc(vartable, varTablesize*sizeof(vartable_t));
      if( vartable == NULL )
	{
          Message("varTablesize = %d", varTablesize);
	  SysError("Reallocation of vartable failed!");
	}
      varID = varTablesize/2;

      for( i = varID; i < varTablesize; i++ )
	vartable[i].param = UNDEF_PARAM;
    }

  paramInitEntry(varID, param);

  return (varID);
}


void varAddRecord(int recID, int param, int gridID, int zaxistype, int lbounds,
		  int level1, int level2, int prec,
		  int *pvarID, int *plevelID, int tsteptype, int numavg, int ltype,
		  const char *name, const char *longname, const char *units)
{
  int varID = UNDEFID;
  int levelID = -1;

  if ( ! (cdiSplitLtype105 == 1 && zaxistype == ZAXIS_HEIGHT) )
    varID = varGetEntry(param, zaxistype, ltype);

  if ( varID == UNDEFID )
    {
      nvars++;
      varID = paramNewEntry(param);
      vartable[varID].gridID    = gridID;
      vartable[varID].zaxistype = zaxistype;
      vartable[varID].ltype     = ltype;
      vartable[varID].lbounds   = lbounds;
      if ( tsteptype > 0 ) vartable[varID].tsteptype = tsteptype;
      if ( numavg ) vartable[varID].timave = 1;

      if ( name )     if ( name[0] )     vartable[varID].name     = strdup(name);
      if ( longname ) if ( longname[0] ) vartable[varID].longname = strdup(longname);
      if ( units )    if ( units[0] )    vartable[varID].units    = strdup(units);
    }
  else
    {
      if ( vartable[varID].gridID != gridID )
	{
	  char paramstr[32];
	  cdiParamToString(param, paramstr, sizeof(paramstr));
	  Message("param = %s gridID = %d", paramstr, gridID);
	  Error("horizontal grid must not change for same param!");
	}
      if ( vartable[varID].zaxistype != zaxistype )
	{
	  char paramstr[32];
	  cdiParamToString(param, paramstr, sizeof(paramstr));
	  Message("param = %s zaxistype = %d", paramstr, zaxistype);
	  Error("zaxistype must not change for same param!");
	}
    }

  if ( prec > vartable[varID].prec ) vartable[varID].prec = prec;

  levelID = levelNewEntry(varID, level1, level2);
  vartable[varID].levelTable[levelID].recID = recID;

  if ( CDI_Debug )
    Message("varID = %d  levelID = %d", varID, levelID);

  *pvarID   = varID;
  *plevelID = levelID;
}


int dblcmp(const void *s1, const void *s2)
{
  int cmp = 0;

  if      ( *((double *) s1) < *((double *) s2) ) cmp = -1;
  else if ( *((double *) s1) > *((double *) s2) ) cmp =  1;

  return (cmp);
}


int cmpLevelTable(const void *s1, const void *s2)
{
  int cmp = 0;
  leveltable_t *x = (leveltable_t *) s1;
  leveltable_t *y = (leveltable_t *) s2;
  /*
  printf("%g %g  %d %d\n", x->leve11, y->level1, x, y);
  */
  if      ( x->level1 < y->level1 ) cmp = -1;
  else if ( x->level1 > y->level1 ) cmp =  1;

  return (cmp);
}


typedef struct
{
  int      varid;
  int      param;
  int      ltype;
}
param_t;


static
int cmpparam(const void *s1, const void *s2)
{
  int cmp = 0;
  param_t *x = (param_t *) s1;
  param_t *y = (param_t *) s2;

  if      ( x->param > y->param ) cmp =  1;
  else if ( x->param < y->param ) cmp = -1;

  return (cmp);
}


static
int cmpltype(const void *s1, const void *s2)
{
  int cmp = 0;
  param_t *x = (param_t *) s1;
  param_t *y = (param_t *) s2;

  if      ( x->ltype > y->ltype ) cmp =  1;
  else if ( x->ltype < y->ltype ) cmp = -1;

  return (cmp);
}


void cdiGenVars(int streamID)
{
  int varID, gridID, zaxisID, levelID;
  int instID, modelID, tableID;
  int param, nlevels, zaxistype, lindex, ltype;
  int prec;
  int tsteptype;
  int timave, timaccu;
  int lbounds;
  int comptype;
  char name[256], longname[256], units[256];
  double *dlevels = NULL;
  double *dlevels1 = NULL;
  double *dlevels2 = NULL;
  int vlistID;
  int *varids, index, varid;
  stream_t *streamptr;

  streamptr = stream_to_pointer(streamID);

  vlistID =  streamInqVlist(streamID);

  varids = (int *) malloc(nvars*sizeof(int));
  for ( varID = 0; varID < nvars; varID++ ) varids[varID] = varID;

  if ( streamptr->sortname )
    {
      int index;
      param_t **varInfo;
      varInfo    = (param_t **) malloc(nvars*sizeof(param_t *));
      varInfo[0] = (param_t *)  malloc(nvars*sizeof(param_t));

      for ( index = 1; index < nvars; index++ )
	varInfo[index] = varInfo[0] + index;

      for ( varid = 0; varid < nvars; varid++ )
	{
	  varInfo[varid]->varid = varids[varid];
	  varInfo[varid]->param = vartable[varid].param;
	  varInfo[varid]->ltype = vartable[varid].ltype;
	}
      qsort(varInfo[0], nvars, sizeof(param_t), cmpltype);
      qsort(varInfo[0], nvars, sizeof(param_t), cmpparam);
      for ( varid = 0; varid < nvars; varid++ )
	{
	  varids[varid] = varInfo[varid]->varid;
	}
      free(varInfo[0]);
      free(varInfo);
    }

  for ( index = 0; index < nvars; index++ )
    {
      varid     = varids[index];

      gridID    = vartable[varid].gridID;
      param     = vartable[varid].param;
      nlevels   = vartable[varid].nlevels;
      ltype     = vartable[varid].ltype;
      zaxistype = vartable[varid].zaxistype;
      if ( ltype == 0 && zaxistype == ZAXIS_GENERIC && cdiDefaultLeveltype != -1 )
	zaxistype = cdiDefaultLeveltype;
      lbounds   = vartable[varid].lbounds;
      prec      = vartable[varid].prec;
      instID    = vartable[varid].instID;
      modelID   = vartable[varid].modelID;
      tableID   = vartable[varid].tableID;
      tsteptype = vartable[varid].tsteptype;
      timave    = vartable[varid].timave;
      timaccu   = vartable[varid].timaccu;
      comptype  = vartable[varid].comptype;

      zaxisID = UNDEFID;

      if ( ltype == 0 && zaxistype == ZAXIS_GENERIC && nlevels == 1 &&
	   ! (fabs(vartable[varid].levelTable[0].level1)>0) )
	zaxistype = ZAXIS_SURFACE;

      dlevels = (double *) malloc(nlevels*sizeof(double));

      if ( lbounds && zaxistype != ZAXIS_HYBRID && zaxistype != ZAXIS_HYBRID_HALF )
	for ( levelID = 0; levelID < nlevels; levelID++ )
	  dlevels[levelID] = (vartable[varid].levelTable[levelID].level1 +
	                      vartable[varid].levelTable[levelID].level2)/2;
      else
	for ( levelID = 0; levelID < nlevels; levelID++ )
	  dlevels[levelID] = vartable[varid].levelTable[levelID].level1;

      if ( nlevels > 1 )
	{
	  int linc = FALSE, ldec = FALSE;
	  /* check increasing of levels */
	  for ( levelID = 1; levelID < nlevels; levelID++ )
	    if ( dlevels[levelID] < dlevels[levelID-1] ) break;

	  if ( levelID == nlevels ) linc = TRUE;

	  if ( linc == FALSE )
	    {
	      /* check decreasing of levels */
	      for ( levelID = 1; levelID < nlevels; levelID++ )
		if ( dlevels[levelID] > dlevels[levelID-1] ) break;

	      if ( levelID == nlevels ) ldec = TRUE;

	      if ( ldec == FALSE ||
		   zaxistype == ZAXIS_HYBRID ||
		   zaxistype == ZAXIS_DEPTH_BELOW_LAND )
		{
		  /*
		  qsort(dlevels, nlevels, sizeof(double), dblcmp);
		  */
		  qsort(vartable[varid].levelTable, nlevels, 
			sizeof(leveltable_t), cmpLevelTable);

		  if ( lbounds && zaxistype != ZAXIS_HYBRID && zaxistype != ZAXIS_HYBRID_HALF )
		    for ( levelID = 0; levelID < nlevels; levelID++ )
		      dlevels[levelID] = (vartable[varid].levelTable[levelID].level1 +
					  vartable[varid].levelTable[levelID].level2)/2.;
		  else
		    for ( levelID = 0; levelID < nlevels; levelID++ )
		      dlevels[levelID] = vartable[varid].levelTable[levelID].level1;
		}
	    }
	}

      if ( lbounds )
	{
	  dlevels1 = (double *) malloc(nlevels*sizeof(double));
	  for ( levelID = 0; levelID < nlevels; levelID++ )
	    dlevels1[levelID] = vartable[varid].levelTable[levelID].level1;
	  dlevels2 = (double *) malloc(nlevels*sizeof(double));
	  for ( levelID = 0; levelID < nlevels; levelID++ )
	    dlevels2[levelID] = vartable[varid].levelTable[levelID].level2;
	}

      zaxisID = varDefZaxis(vlistID, zaxistype, nlevels, dlevels, lbounds, dlevels1, dlevels2,
			    Vctsize, Vct, NULL, NULL, NULL, 0, 0, ltype);

      if ( lbounds ) free(dlevels1);
      if ( lbounds ) free(dlevels2);
      free(dlevels);

      varID = streamNewVar(streamID, gridID, zaxisID);
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);

      vlistDefVarParam(vlistID, varID, param);
      vlistDefVarDatatype(vlistID, varID, prec);
      vlistDefVarTsteptype(vlistID, varID, tsteptype);
      vlistDefVarTimave(vlistID, varID, timave);
      vlistDefVarTimaccu(vlistID, varID, timaccu);
      vlistDefVarCompType(vlistID, varID, comptype);

      if ( vartable[varid].lmissval ) vlistDefVarMissval(vlistID, varID, vartable[varid].missval);

      if ( vartable[varid].name )     vlistDefVarName(vlistID, varID, vartable[varid].name);
      if ( vartable[varid].longname ) vlistDefVarLongname(vlistID, varID, vartable[varid].longname);
      if ( vartable[varid].units )    vlistDefVarUnits(vlistID, varID, vartable[varid].units);

      if ( cdiDefaultTableID != UNDEFID )
	{
	  int pdis, pcat, pnum;
	  cdiDecodeParam(param, &pnum, &pcat, &pdis);
	  if ( tableInqParNamePtr(cdiDefaultTableID, pnum) )
	    {
	      if ( tableID != UNDEFID )
		{
		  strcpy(name, tableInqParNamePtr(cdiDefaultTableID, pnum));
		  vlistDefVarName(vlistID, varID, name);
		  if ( tableInqParLongnamePtr(cdiDefaultTableID, pnum) )
		    {
		      strcpy(longname, tableInqParLongnamePtr(cdiDefaultTableID, pnum));
		      vlistDefVarLongname(vlistID, varID, longname);
		    }
		  if ( tableInqParUnitsPtr(cdiDefaultTableID, pnum) )
		    {
		      strcpy(units, tableInqParUnitsPtr(cdiDefaultTableID, pnum));
		      vlistDefVarUnits(vlistID, varID, units);
		    }
		}
	      else
		tableID = cdiDefaultTableID;
	    }
	  if ( cdiDefaultModelID != UNDEFID ) modelID = cdiDefaultModelID;
	  if ( cdiDefaultInstID  != UNDEFID )  instID = cdiDefaultInstID;
	}

      if ( instID  != UNDEFID ) vlistDefVarInstitut(vlistID, varID, instID);
      if ( modelID != UNDEFID ) vlistDefVarModel(vlistID, varID, modelID);
      if ( tableID != UNDEFID ) vlistDefVarTable(vlistID, varID, tableID);
    }

  for ( index = 0; index < nvars; index++ )
    {
      varID     = index;
      varid     = varids[index];

      nlevels   = vartable[varid].nlevels;
      /*
      for ( levelID = 0; levelID < nlevels; levelID++ )
	{
	  lindex = vartable[varid].levelTable[levelID].lindex;
	  printf("%d %d %d %d %d\n", varID, levelID, 
		 vartable[varid].levelTable[levelID].lindex,
		 vartable[varid].levelTable[levelID].recID,
		 vartable[varid].levelTable[levelID].level1);
	}
      */
      for ( levelID = 0; levelID < nlevels; levelID++ )
	{
	  streamptr->vars[varID].level[levelID] =
	    vartable[varid].levelTable[levelID].recID;
	  for ( lindex = 0; lindex < nlevels; lindex++ )
	    if ( levelID == vartable[varid].levelTable[lindex].lindex ) break;

	  if ( lindex == nlevels )
	    Error("Internal problem! lindex not found.");

	  streamptr->vars[varID].lindex[levelID] = lindex;
	}
    }

  free(varids);

  varFree();
}


void varDefVCT(size_t vctsize, double *vctptr)
{
  if ( Vct == NULL && vctptr != NULL && vctsize > 0 )
    {
      Vctsize = vctsize;
      Vct = (double *) malloc(vctsize*sizeof(double));
      memcpy(Vct, vctptr, vctsize*sizeof(double)); 
    }
}


int varDefGrid(int vlistID, grid_t grid, int mode)
{
  /*
    mode: 0 search in vlist and grid table
          1 search in grid table
   */
  int gridglobdefined = FALSE;
  int griddefined;
  int ngrids;
  int gridID = UNDEFID;
  int index;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  griddefined = FALSE;
  ngrids = vlistptr->ngrids;

  if ( mode == 0 )
    for ( index = 0; index < ngrids; index++ )
      {
	gridID = vlistptr->gridIDs[index];
	if ( gridID == UNDEFID )
	  Error("Internal problem: undefined gridID %d!", gridID);

	if ( gridCompare(gridID, grid) == 0 )
	  {
	    griddefined = TRUE;
	    break;
	  }
      }

  if ( ! griddefined )
    {
      ngrids = gridSize();
      for ( gridID = 0; gridID < ngrids; gridID++ )
	{
	  if ( gridCompare(gridID, grid) == 0 )
	    {
	      gridglobdefined = TRUE;
	      break;
	    }
	}

      ngrids = vlistptr->ngrids;
      if ( mode == 1 )
	for ( index = 0; index < ngrids; index++ )
	  if ( vlistptr->gridIDs[index] == gridID )
	    {
	      gridglobdefined = FALSE;
	      break;
	    }
    }

  if ( ! griddefined )
    {
      if ( ! gridglobdefined ) gridID = gridGenerate(grid);
      ngrids = vlistptr->ngrids;
      vlistptr->gridIDs[ngrids] = gridID;
      vlistptr->ngrids++;
    }

  return (gridID);
}


int zaxisCompare(int zaxisID, int zaxistype, int nlevels, int lbounds, double *levels, char *longname, char *units, int ltype)
{
  int differ = 1;
  int levelID;
  int zlbounds = 0;
  int ltype_is_equal = FALSE;

  if ( ltype == zaxisInqLtype(zaxisID) ) ltype_is_equal = TRUE;

  if ( ltype_is_equal && (zaxistype == zaxisInqType(zaxisID) || zaxistype == ZAXIS_GENERIC) )
    {
      if ( zaxisInqLbounds(zaxisID, NULL) > 0 ) zlbounds = 1;
      if ( nlevels == zaxisInqSize(zaxisID) && zlbounds == lbounds )
	{
	  const double *dlevels;
	  char zlongname[256];
	  char zunits[256];

	  dlevels = zaxisInqLevelsPtr(zaxisID);
	  for ( levelID = 0; levelID < nlevels; levelID++ )
	    {
	      if ( fabs(dlevels[levelID] - levels[levelID]) > 1.e-9 )
		break;
	    }

	  if ( levelID == nlevels ) differ = 0;

	  if ( ! differ )
	    {
	      zaxisInqLongname(zaxisID, zlongname);
	      zaxisInqUnits(zaxisID, zunits);
	      if ( longname && zlongname[0] )
		{
		  if ( strcmp(longname, zlongname) != 0 ) differ = 1;
		}
	      if ( units && zunits[0] )
		{
		  if ( strcmp(units, zunits) != 0 ) differ = 1;
		}
	    }
	}
    }

  return (differ);
}


int varDefZaxis(int vlistID, int zaxistype, int nlevels, double *levels, int lbounds,
		double *levels1, double *levels2, int vctsize, double *vct, char *name,
		char *longname, char *units, int prec, int mode, int ltype)
{
  /*
    mode: 0 search in vlist and zaxis table
          1 search in zaxis table
   */
  int zaxisdefined;
  int nzaxis;
  int zaxisID = UNDEFID;
  int index;
  int zaxisglobdefined = 0;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  zaxisdefined = 0;
  nzaxis = vlistptr->nzaxis;

  if ( mode == 0 )
    for ( index = 0; index < nzaxis; index++ )
      {
	zaxisID = vlistptr->zaxisIDs[index];

	if ( zaxisCompare(zaxisID, zaxistype, nlevels, lbounds, levels, longname, units, ltype) == 0 )
	  {
	    zaxisdefined = 1;
	    break;
	  }
      }

  if ( ! zaxisdefined )
    {
      nzaxis = zaxisSize();
      for ( zaxisID = 0; zaxisID < nzaxis; zaxisID++ )
	if ( zaxisCompare(zaxisID, zaxistype, nlevels, lbounds, levels, longname, units, ltype) == 0 )
	  {
	    zaxisglobdefined = 1;
	    break;
	  }

      nzaxis = vlistptr->nzaxis;
      if ( mode == 1 )
	for ( index = 0; index < nzaxis; index++ )
	  if ( vlistptr->zaxisIDs[index] == zaxisID )
	    {
	      zaxisglobdefined = FALSE;
	      break;
	    }
    }

  if ( ! zaxisdefined )
    {
      if ( ! zaxisglobdefined )
	{
	  zaxisID = zaxisCreate(zaxistype, nlevels);
	  zaxisDefLevels(zaxisID, levels);
	  if ( lbounds )
	    {
	      zaxisDefLbounds(zaxisID, levels1);
	      zaxisDefUbounds(zaxisID, levels2);
	    }

	  if ( zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF )
	    {
	      /* if ( vctsize > 0 && vctsize >= 2*(nlevels+1)) */
	      /* if ( vctsize > 0 && vctsize >= 2*(nlevels)) */
	      if ( vctsize > 0 )
		zaxisDefVct(zaxisID, vctsize, vct);
	      else
		Warning("VCT missing");
	    }

	  zaxisDefName(zaxisID, name);
	  zaxisDefLongname(zaxisID, longname);
	  zaxisDefUnits(zaxisID, units);
	  zaxisDefPrec(zaxisID, prec);
	  zaxisDefLtype(zaxisID, ltype);
	}

      nzaxis = vlistptr->nzaxis;
      vlistptr->zaxisIDs[nzaxis] = zaxisID;
      vlistptr->nzaxis++;
    }

  return (zaxisID);
}


void varDefMissval(int varID, double missval)
{
  vartable[varID].lmissval = 1;
  vartable[varID].missval = missval;
}


void varDefCompType(int varID, int comptype)
{
  if ( vartable[varID].comptype == COMPRESS_NONE )
    vartable[varID].comptype = comptype;
}


void varDefCompLevel(int varID, int complevel)
{
  vartable[varID].complevel = complevel;
}


int varInqInst(int varID)
{
  return (vartable[varID].instID);
}


void varDefInst(int varID, int instID)
{
  vartable[varID].instID = instID;
}


int varInqModel(int varID)
{
  return (vartable[varID].modelID);
}


void varDefModel(int varID, int modelID)
{
  vartable[varID].modelID = modelID;
}


int varInqTable(int varID)
{
  return (vartable[varID].tableID);
}


void varDefTable(int varID, int tableID)
{
  vartable[varID].tableID = tableID;
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
