#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "cdi.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "resource_handle.h"
#include "varscan.h"
#include "vlist.h"
#include "grid.h"
#include "zaxis.h"


#undef  UNDEFID
#define UNDEFID -1

static size_t Vctsize = 0;
static double *Vct = NULL;

static int numberOfVerticalLevels = 0;
static int numberOfVerticalGrid = 0;
static unsigned char uuidVGrid[CDI_UUID_SIZE];

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
  int           param;
  int           prec;
  int           tsteptype;
  int           timave;
  int           timaccu;
  int           gridID;
  int           zaxistype;
  int           ltype1;     /* GRIB first level type */
  int           ltype2;     /* GRIB second level type */
  int           lbounds;
  int           level_sf;
  int           level_unit;
  int           zaxisID;
  unsigned      nlevels;
  int           levelTableSize;
  leveltable_t *levelTable;
  int           instID;
  int           modelID;
  int           tableID;
  int           comptype;       // compression type
  int           complevel;      // compression level
  int           lmissval;
  double        missval;
  char         *name;
  char         *stdname;
  char         *longname;
  char         *units;
  ensinfo_t    *ensdata;
  int           typeOfGeneratingProcess;
  int           productDefinitionTemplate;
#if  defined  (HAVE_LIBGRIB_API)
  /* (Optional) list of keyword/double value pairs */
  int           opt_grib_dbl_nentries;
  char         *opt_grib_dbl_keyword[MAX_OPT_GRIB_ENTRIES];
  double        opt_grib_dbl_val[MAX_OPT_GRIB_ENTRIES];
  /* (Optional) list of keyword/integer value pairs */
  int           opt_grib_int_nentries;
  char         *opt_grib_int_keyword[MAX_OPT_GRIB_ENTRIES];
  int           opt_grib_int_val[MAX_OPT_GRIB_ENTRIES];
#endif
}
vartable_t;


static vartable_t *vartable;
static unsigned varTablesize = 0;
static unsigned nvars = 0;


static void
paramInitEntry(unsigned varID, int param)
{
  vartable[varID].param          = param;
  vartable[varID].prec           = 0;
  vartable[varID].tsteptype      = TSTEP_INSTANT;
  vartable[varID].timave         = 0;
  vartable[varID].timaccu        = 0;
  vartable[varID].gridID         = UNDEFID;
  vartable[varID].zaxistype      = 0;
  vartable[varID].ltype1         = 0;
  vartable[varID].ltype2         = -1;
  vartable[varID].lbounds        = 0;
  vartable[varID].level_sf       = 0;
  vartable[varID].level_unit     = 0;
  vartable[varID].levelTable     = NULL;
  vartable[varID].levelTableSize = 0;
  vartable[varID].nlevels        = 0;
  vartable[varID].instID         = UNDEFID;
  vartable[varID].modelID        = UNDEFID;
  vartable[varID].tableID        = UNDEFID;
  vartable[varID].typeOfGeneratingProcess   = UNDEFID;
  vartable[varID].productDefinitionTemplate = UNDEFID;
  vartable[varID].comptype       = COMPRESS_NONE;
  vartable[varID].complevel      = 1;
  vartable[varID].lmissval       = 0;
  vartable[varID].missval        = 0;
  vartable[varID].name           = NULL;
  vartable[varID].stdname        = NULL;
  vartable[varID].longname       = NULL;
  vartable[varID].units          = NULL;
  vartable[varID].ensdata        = NULL;
}

static unsigned
varGetEntry(int param, int zaxistype, int ltype1, int tsteptype, const char *name)
{
  for ( unsigned varID = 0; varID < varTablesize; varID++ )
    {
      if ( vartable[varID].param      == param       &&
	   vartable[varID].zaxistype  == zaxistype   &&
	   vartable[varID].ltype1     == ltype1      &&
	   vartable[varID].tsteptype  == tsteptype )
        {
          if ( name && name[0] && vartable[varID].name && vartable[varID].name[0] )
            {
              if ( strcmp(name, vartable[varID].name) == 0 ) return (varID);
            }
          else
            {
              return (varID);
            }
        }
    }

  return (unsigned)-1;
}

static
void varFree(void)
{
  for ( unsigned varID = 0; varID < nvars; varID++ )
    {
      if ( vartable[varID].levelTable )
	free(vartable[varID].levelTable);

      if ( vartable[varID].name )     free(vartable[varID].name);
      if ( vartable[varID].stdname )  free(vartable[varID].stdname);
      if ( vartable[varID].longname ) free(vartable[varID].longname);
      if ( vartable[varID].units )    free(vartable[varID].units);
      if ( vartable[varID].ensdata )  free(vartable[varID].ensdata);
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

static int
levelNewEntry(unsigned varID, int level1, int level2)
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
      levelTable = (leveltable_t *)xmalloc((size_t)levelTableSize
                                           * sizeof (leveltable_t));
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
      levelTable = (leveltable_t *)xrealloc(levelTable, (size_t)levelTableSize
                                            * sizeof (leveltable_t));
      if( levelTable == NULL )
	{
          Message("levelTableSize = %d", levelTableSize);
	  SysError("Reallocation of leveltable failed");
	}
      levelID = levelTableSize/2;

      for( i = levelID; i < levelTableSize; i++ )
	levelTable[i].recID = UNDEFID;
    }

  levelTable[levelID].level1   = level1;
  levelTable[levelID].level2   = level2;
  levelTable[levelID].lindex   = levelID;

  vartable[varID].nlevels = (unsigned)levelID+1;
  vartable[varID].levelTableSize = levelTableSize;
  vartable[varID].levelTable = levelTable;

  return (levelID);
}

#define  UNDEF_PARAM  -4711

static unsigned
paramNewEntry(int param)
{
  unsigned varID = 0;

  /*
    Look for a free slot in vartable.
    (Create the table the first time through).
  */
  if ( ! varTablesize )
    {
      varTablesize = 2;
      vartable = (vartable_t *)xmalloc((size_t)varTablesize
                                       * sizeof (vartable_t));
      if( vartable == NULL )
	{
          Message("varTablesize = %d", varTablesize);
	  SysError("Allocation of vartable failed");
	}

      for( unsigned i = 0; i < varTablesize; i++ )
	{
	  vartable[i].param = UNDEF_PARAM;
#if  defined  (HAVE_LIBGRIB_API)
	  vartable[i].opt_grib_int_nentries = 0;
	  vartable[i].opt_grib_dbl_nentries = 0;
#endif
	}
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

      varTablesize = 2 * varTablesize;
      vartable = (vartable_t *)xrealloc(vartable, (size_t)varTablesize
                                        * sizeof (vartable_t));
      if( vartable == NULL )
	{
          Message("varTablesize = %d", varTablesize);
	  SysError("Reallocation of vartable failed!");
	}
      varID = varTablesize/2;

      for( unsigned i = varID; i < varTablesize; i++ )
	{
	  vartable[i].param = UNDEF_PARAM;
#if  defined  (HAVE_LIBGRIB_API)
	  vartable[i].opt_grib_int_nentries = 0;
	  vartable[i].opt_grib_dbl_nentries = 0;
#endif
	}
    }

  paramInitEntry(varID, param);

  return (varID);
}


void varAddRecord(int recID, int param, int gridID, int zaxistype, int lbounds,
		  int level1, int level2, int level_sf, int level_unit, int prec,
		  int *pvarID, int *plevelID, int tsteptype, int numavg, int ltype1, int ltype2,
		  const char *name, const char *stdname, const char *longname, const char *units)
{
  unsigned varID = (cdiSplitLtype105 != 1 || zaxistype != ZAXIS_HEIGHT) ?
    varGetEntry(param, zaxistype, ltype1, tsteptype, name) : (unsigned)UNDEFID;

  if ( varID == (unsigned)UNDEFID )
    {
      nvars++;
      varID = paramNewEntry(param);
      vartable[varID].gridID     = gridID;
      vartable[varID].zaxistype  = zaxistype;
      vartable[varID].ltype1     = ltype1;
      vartable[varID].ltype2     = ltype2;
      vartable[varID].lbounds    = lbounds;
      vartable[varID].level_sf   = level_sf;
      vartable[varID].level_unit = level_unit;
      vartable[varID].tsteptype  = tsteptype;
      if ( numavg ) vartable[varID].timave = 1;

      if ( name )     if ( name[0] )     vartable[varID].name     = strdup(name);
      if ( stdname )  if ( stdname[0] )  vartable[varID].stdname  = strdup(stdname);
      if ( longname ) if ( longname[0] ) vartable[varID].longname = strdup(longname);
      if ( units )    if ( units[0] )    vartable[varID].units    = strdup(units);
    }
  else
    {
      char paramstr[32];
      cdiParamToString(param, paramstr, sizeof(paramstr));

      if ( vartable[varID].gridID != gridID )
	{
	  Message("param = %s gridID = %d", paramstr, gridID);
	  Error("horizontal grid must not change for same parameter!");
	}
      if ( vartable[varID].zaxistype != zaxistype )
	{
	  Message("param = %s zaxistype = %d", paramstr, zaxistype);
	  Error("zaxistype must not change for same parameter!");
	}
    }

  if ( prec > vartable[varID].prec ) vartable[varID].prec = prec;

  int levelID = levelNewEntry(varID, level1, level2);
  vartable[varID].levelTable[levelID].recID = recID;

  *pvarID   = (int)varID;
  *plevelID = levelID;
}
/*
static
int dblcmp(const void *s1, const void *s2)
{
  int cmp = 0;

  if      ( *((double *) s1) < *((double *) s2) ) cmp = -1;
  else if ( *((double *) s1) > *((double *) s2) ) cmp =  1;

  return (cmp);
}
*/
static
int cmpLevelTable(const void* s1, const void* s2)
{
  int cmp = 0;
  const leveltable_t* x = (const leveltable_t*) s1;
  const leveltable_t* y = (const leveltable_t*) s2;
  /*
  printf("%g %g  %d %d\n", x->leve11, y->level1, x, y);
  */
  if      ( x->level1 < y->level1 ) cmp = -1;
  else if ( x->level1 > y->level1 ) cmp =  1;

  return (cmp);
}

static
int cmpLevelTableInv(const void* s1, const void* s2)
{
  int cmp = 0;
  const leveltable_t* x = (const leveltable_t*) s1;
  const leveltable_t* y = (const leveltable_t*) s2;
  /*
  printf("%g %g  %d %d\n", x->leve11, y->level1, x, y);
  */
  if      ( x->level1 < y->level1 ) cmp =  1;
  else if ( x->level1 > y->level1 ) cmp = -1;

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
int cmpparam(const void* s1, const void* s2)
{
  const param_t* x = (const param_t*) s1;
  const param_t* y = (const param_t*) s2;

  int cmp = (( x->param > y->param ) - ( x->param < y->param )) * 2
           + ( x->ltype > y->ltype ) - ( x->ltype < y->ltype );

  return (cmp);
}


void cdi_generate_vars(stream_t *streamptr)
{
  int gridID, zaxisID;
  int instID, modelID, tableID;
  int param, zaxistype, ltype1, ltype2;
  int prec;
  int tsteptype;
  int timave, timaccu;
  int lbounds;
  int comptype;
  char name[CDI_MAX_NAME], longname[CDI_MAX_NAME], units[CDI_MAX_NAME];
  double *dlevels = NULL;
  double *dlevels1 = NULL;
  double *dlevels2 = NULL;
  double level_sf = 1;
  int vlistID = streamptr->vlistID;

  int *varids = (int *)xmalloc(nvars*sizeof(int));
  for ( unsigned varID = 0; varID < nvars; varID++ ) varids[varID] = (int)varID;

  if ( streamptr->sortname )
    {
      param_t *varInfo = (param_t *)xmalloc((size_t)nvars * sizeof (param_t));

      for ( unsigned varID = 0; varID < nvars; varID++ )
	{
	  varInfo[varID].varid = varids[varID];
	  varInfo[varID].param = vartable[varID].param;
	  varInfo[varID].ltype = vartable[varID].ltype1;
	}
      qsort(varInfo, (size_t)nvars, sizeof(param_t), cmpparam);
      for ( unsigned varID = 0; varID < nvars; varID++ )
	{
	  varids[varID] = varInfo[varID].varid;
	}
      free(varInfo);
    }

  for ( unsigned index = 0; index < nvars; index++ )
    {
      int varid      = varids[index];

      gridID     = vartable[varid].gridID;
      param      = vartable[varid].param;
      unsigned nlevels = vartable[varid].nlevels;
      ltype1     = vartable[varid].ltype1;
      ltype2     = vartable[varid].ltype2;
      zaxistype = vartable[varid].zaxistype;
      if ( ltype1 == 0 && zaxistype == ZAXIS_GENERIC && cdiDefaultLeveltype != -1 )
	zaxistype = cdiDefaultLeveltype;
      lbounds    = vartable[varid].lbounds;
      prec       = vartable[varid].prec;
      instID     = vartable[varid].instID;
      modelID    = vartable[varid].modelID;
      tableID    = vartable[varid].tableID;
      tsteptype  = vartable[varid].tsteptype;
      timave     = vartable[varid].timave;
      timaccu    = vartable[varid].timaccu;
      comptype   = vartable[varid].comptype;

      level_sf  = 1;
      if ( vartable[varid].level_sf != 0 ) level_sf = 1./vartable[varid].level_sf;

      zaxisID = UNDEFID;

      if ( ltype1 == 0 && zaxistype == ZAXIS_GENERIC && nlevels == 1 &&
	   vartable[varid].levelTable[0].level1 == 0 )
	zaxistype = ZAXIS_SURFACE;

      dlevels = (double *) malloc(nlevels*sizeof(double));

      if ( lbounds && zaxistype != ZAXIS_HYBRID && zaxistype != ZAXIS_HYBRID_HALF )
	for (unsigned levelID = 0; levelID < nlevels; levelID++ )
	  dlevels[levelID] = (level_sf*vartable[varid].levelTable[levelID].level1 +
	                      level_sf*vartable[varid].levelTable[levelID].level2)/2;
      else
	for (unsigned levelID = 0; levelID < nlevels; levelID++ )
	  dlevels[levelID] = level_sf*vartable[varid].levelTable[levelID].level1;

      if ( nlevels > 1 )
	{
          bool linc = true, ldec = true, lsort = false;
          for (unsigned levelID = 1; levelID < nlevels; levelID++ )
            {
              /* check increasing of levels */
              linc &= (dlevels[levelID] > dlevels[levelID-1]);
              /* check decreasing of levels */
              ldec &= (dlevels[levelID] < dlevels[levelID-1]);
            }
          /*
           * always sort pressure z-axis to ensure
           * vartable[varid].levelTable[levelID1].level1 < vartable[varid].levelTable[levelID2].level1 <=> levelID1 > levelID2
           * unless already sorted in decreasing order
           */
          if ( !ldec && zaxistype == ZAXIS_PRESSURE )
            {
              qsort(vartable[varid].levelTable, nlevels, sizeof(leveltable_t), cmpLevelTableInv);
              lsort = true;
            }
          /*
           * always sort hybrid and depth-below-land z-axis to ensure
           * vartable[varid].levelTable[levelID1].level1 < vartable[varid].levelTable[levelID2].level1 <=> levelID1 < levelID2
           * unless already sorted in increasing order
           */
          else if ( (!linc && !ldec) ||
                    zaxistype == ZAXIS_HYBRID ||
                    zaxistype == ZAXIS_DEPTH_BELOW_LAND )
            {
              qsort(vartable[varid].levelTable, nlevels, sizeof(leveltable_t), cmpLevelTable);
              lsort = true;
            }

          if ( lsort )
            {
              if ( lbounds && zaxistype != ZAXIS_HYBRID && zaxistype != ZAXIS_HYBRID_HALF )
                for (unsigned levelID = 0; levelID < nlevels; levelID++ )
                  dlevels[levelID] = (level_sf*vartable[varid].levelTable[levelID].level1 +
                                      level_sf*vartable[varid].levelTable[levelID].level2)/2.;
              else
                for (unsigned levelID = 0; levelID < nlevels; levelID++ )
                  dlevels[levelID] = level_sf*vartable[varid].levelTable[levelID].level1;
            }
	}

      if ( lbounds )
	{
	  dlevels1 = (double *) malloc(nlevels*sizeof(double));
	  for (unsigned levelID = 0; levelID < nlevels; levelID++)
	    dlevels1[levelID] = level_sf*vartable[varid].levelTable[levelID].level1;
	  dlevels2 = (double *) malloc(nlevels*sizeof(double));
	  for (unsigned levelID = 0; levelID < nlevels; levelID++)
	    dlevels2[levelID] = level_sf*vartable[varid].levelTable[levelID].level2;
        }

      char *unitptr = cdiUnitNamePtr(vartable[varid].level_unit);
      zaxisID = varDefZaxis(vlistID, zaxistype, (int)nlevels, dlevels, lbounds, dlevels1, dlevels2,
                            (int)Vctsize, Vct, NULL, NULL, unitptr, 0, 0, ltype1);

      if ( ltype1 != ltype2 && ltype2 != -1 )
        {
          zaxisDefLtype2(zaxisID, ltype2);
        }

      if ( zaxisInqType(zaxisID) == ZAXIS_REFERENCE )
        {
          if ( numberOfVerticalLevels > 0 ) zaxisDefNlevRef(zaxisID, numberOfVerticalLevels);
          if ( numberOfVerticalGrid > 0 ) zaxisDefNumber(zaxisID, numberOfVerticalGrid);
          if ( !cdiUUIDIsNull(uuidVGrid) ) zaxisDefUUID(zaxisID, uuidVGrid);
        }

      if ( lbounds ) free(dlevels1);
      if ( lbounds ) free(dlevels2);
      free(dlevels);

      int varID = stream_new_var(streamptr, gridID, zaxisID);
      varID = vlistDefVar(vlistID, gridID, zaxisID, tsteptype);

      vlistDefVarParam(vlistID, varID, param);
      vlistDefVarDatatype(vlistID, varID, prec);
      vlistDefVarTimave(vlistID, varID, timave);
      vlistDefVarTimaccu(vlistID, varID, timaccu);
      vlistDefVarCompType(vlistID, varID, comptype);

      if ( vartable[varid].typeOfGeneratingProcess != UNDEFID )
        vlistDefVarTypeOfGeneratingProcess(vlistID, varID, vartable[varid].typeOfGeneratingProcess);

      if ( vartable[varid].productDefinitionTemplate != UNDEFID )
        vlistDefVarProductDefinitionTemplate(vlistID, varID, vartable[varid].productDefinitionTemplate);

      if ( vartable[varid].lmissval ) vlistDefVarMissval(vlistID, varID, vartable[varid].missval);

      if ( vartable[varid].name )     vlistDefVarName(vlistID, varID, vartable[varid].name);
      if ( vartable[varid].stdname )  vlistDefVarStdname(vlistID, varID, vartable[varid].stdname);
      if ( vartable[varid].longname ) vlistDefVarLongname(vlistID, varID, vartable[varid].longname);
      if ( vartable[varid].units )    vlistDefVarUnits(vlistID, varID, vartable[varid].units);

      if ( vartable[varid].ensdata )  vlistDefVarEnsemble(vlistID, varID, vartable[varid].ensdata->ens_index,
	                                                  vartable[varid].ensdata->ens_count,
							  vartable[varid].ensdata->forecast_init_type);

#if  defined  (HAVE_LIBGRIB_API)
      /* ---------------------------------- */
      /* Local change: 2013-04-23, FP (DWD) */
      /* ---------------------------------- */

      int    i;
      vlist_t *vlistptr;
      vlistptr = vlist_to_pointer(vlistID);
      for (i=0; i<vartable[varid].opt_grib_int_nentries; i++)
        {
          int idx = vlistptr->vars[varID].opt_grib_int_nentries;
          vlistptr->vars[varID].opt_grib_int_nentries++;
          if ( idx >= MAX_OPT_GRIB_ENTRIES ) Error("Too many optional keyword/integer value pairs!");
          vlistptr->vars[varID].opt_grib_int_update[idx] = TRUE;
          vlistptr->vars[varID].opt_grib_int_val[idx] = vartable[varid].opt_grib_int_val[idx];
          vlistptr->vars[varID].opt_grib_int_keyword[idx] = strdupx(vartable[varid].opt_grib_int_keyword[idx]);
        }
      for (i=0; i<vartable[varid].opt_grib_dbl_nentries; i++)
        {
          int idx = vlistptr->vars[varID].opt_grib_dbl_nentries;
          vlistptr->vars[varID].opt_grib_dbl_nentries++;
          if ( idx >= MAX_OPT_GRIB_ENTRIES ) Error("Too many optional keyword/double value pairs!");
          vlistptr->vars[varID].opt_grib_dbl_update[idx] = TRUE;
          vlistptr->vars[varID].opt_grib_dbl_val[idx] = vartable[varid].opt_grib_dbl_val[idx];
          vlistptr->vars[varID].opt_grib_dbl_keyword[idx] = strdupx(vartable[varid].opt_grib_dbl_keyword[idx]);
        }
      /* note: if the key is not defined, we do not throw an error! */
#endif

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

  for ( unsigned index = 0; index < nvars; index++ )
    {
      int varID = (int)index;
      int varid = varids[index];

      unsigned nlevels = vartable[varid].nlevels;
      /*
      for ( levelID = 0; levelID < nlevels; levelID++ )
	{
	  printf("%d %d %d %d %d\n", varID, levelID,
		 vartable[varid].levelTable[levelID].lindex,
		 vartable[varid].levelTable[levelID].recID,
		 vartable[varid].levelTable[levelID].level1);
	}
      */
      for (unsigned levelID = 0; levelID < nlevels; levelID++)
	{
	  streamptr->vars[varID].level[levelID] = vartable[varid].levelTable[levelID].recID;
          unsigned lindex;
	  for (lindex = 0; lindex < nlevels; lindex++ )
	    if ( levelID == (unsigned)vartable[varid].levelTable[lindex].lindex ) break;

	  if ( lindex == nlevels )
	    Error("Internal problem! lindex not found.");

	  streamptr->vars[varID].lindex[levelID] = (int)lindex;
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


void varDefZAxisReference(int nhlev, int nvgrid, unsigned char uuid[CDI_UUID_SIZE])
{
  numberOfVerticalLevels = nhlev;
  numberOfVerticalGrid = nvgrid;
  memcpy(uuidVGrid, uuid, CDI_UUID_SIZE);
}

struct varDefGridSearchState
{
  int resIDValue;
  const grid_t *queryKey;
};

static enum cdiApplyRet
varDefGridSearch(int id, void *res, void *data)
{
  struct varDefGridSearchState *state = data;
  (void)res;
  if (gridCompare(id, state->queryKey) == 0)
    {
      state->resIDValue = id;
      return CDI_APPLY_STOP;
    }
  else
    return CDI_APPLY_GO_ON;
}

int varDefGrid(int vlistID, const grid_t *grid, int mode)
{
  /*
    mode: 0 search in vlist and grid table
          1 search in grid table
   */
  int gridglobdefined = FALSE;
  int griddefined;
  int gridID = CDI_UNDEFID;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  griddefined = FALSE;
  unsigned ngrids = (unsigned)vlistptr->ngrids;

  if ( mode == 0 )
    for (unsigned index = 0; index < ngrids; index++ )
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
      struct varDefGridSearchState query = { .queryKey = grid };
      if ((gridglobdefined
           = (cdiResHFilterApply(&gridOps, varDefGridSearch, &query)
              == CDI_APPLY_STOP)))
        gridID = query.resIDValue;

      if ( mode == 1 && gridglobdefined)
	for (unsigned index = 0; index < ngrids; index++ )
	  if ( vlistptr->gridIDs[index] == gridID )
	    {
	      gridglobdefined = FALSE;
	      break;
	    }
    }

  if ( ! griddefined )
    {
      if ( ! gridglobdefined ) gridID = gridGenerate(grid);
      ngrids = (unsigned)vlistptr->ngrids;
      vlistptr->gridIDs[ngrids] = gridID;
      vlistptr->ngrids++;
    }

  return (gridID);
}


int zaxisCompare(int zaxisID, int zaxistype, int nlevels, int lbounds, const double *levels, char *longname, char *units, int ltype1)
{
  int differ = 1;
  int levelID;
  int zlbounds = 0;
  int ltype_is_equal = FALSE;

  if ( ltype1 == zaxisInqLtype(zaxisID) ) ltype_is_equal = TRUE;

  if ( ltype_is_equal && (zaxistype == zaxisInqType(zaxisID) || zaxistype == ZAXIS_GENERIC) )
    {
      if ( zaxisInqLbounds(zaxisID, NULL) > 0 ) zlbounds = 1;
      if ( nlevels == zaxisInqSize(zaxisID) && zlbounds == lbounds )
	{
	  const double *dlevels;
	  char zlongname[CDI_MAX_NAME];
	  char zunits[CDI_MAX_NAME];

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

struct varDefZAxisSearchState
{
  int resIDValue;
  int zaxistype;
  int nlevels;
  double *levels;
  int lbounds;
  char *longname, *units;
  int ltype;
};

static enum cdiApplyRet
varDefZAxisSearch(int id, void *res, void *data)
{
  struct varDefZAxisSearchState *state = data;
  (void)res;
  if (zaxisCompare(id, state->zaxistype, state->nlevels, state->lbounds,
                   state->levels, state->longname, state->units, state->ltype)
      == 0)
    {
      state->resIDValue = id;
      return CDI_APPLY_STOP;
    }
  else
    return CDI_APPLY_GO_ON;
}


int varDefZaxis(int vlistID, int zaxistype, int nlevels, double *levels, int lbounds,
		double *levels1, double *levels2, int vctsize, double *vct, char *name,
		char *longname, char *units, int prec, int mode, int ltype1)
{
  /*
    mode: 0 search in vlist and zaxis table
          1 search in zaxis table
   */
  int zaxisdefined = 0;
  int nzaxis;
  int zaxisID = UNDEFID;
  int index;
  int zaxisglobdefined = 0;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  nzaxis = vlistptr->nzaxis;

  if ( mode == 0 )
    for ( index = 0; index < nzaxis; index++ )
      {
	zaxisID = vlistptr->zaxisIDs[index];

	if ( zaxisCompare(zaxisID, zaxistype, nlevels, lbounds, levels, longname, units, ltype1) == 0 )
	  {
	    zaxisdefined = 1;
	    break;
	  }
      }

  if ( ! zaxisdefined )
    {
      struct varDefZAxisSearchState query = {
        .zaxistype = zaxistype,
        .nlevels = nlevels,
        .levels = levels,
        .lbounds = lbounds,
        .longname = longname,
        .units = units,
        .ltype = ltype1,
      };
      if ((zaxisglobdefined
           = (cdiResHFilterApply(&zaxisOps, varDefZAxisSearch, &query)
              == CDI_APPLY_STOP)))
        zaxisID = query.resIDValue;

      if ( mode == 1 && zaxisglobdefined)
	for (int index = 0; index < nzaxis; index++ )
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
	  zaxisDefLtype(zaxisID, ltype1);
	}

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


void varDefEnsembleInfo(int varID, int ens_idx, int ens_count, int forecast_type)
{
  if ( vartable[varID].ensdata == NULL )
      vartable[varID].ensdata = (ensinfo_t *)xmalloc( sizeof( ensinfo_t ) );

  vartable[varID].ensdata->ens_index = ens_idx;
  vartable[varID].ensdata->ens_count = ens_count;
  vartable[varID].ensdata->forecast_init_type = forecast_type;
}


void varDefTypeOfGeneratingProcess(int varID, int typeOfGeneratingProcess)
{
  vartable[varID].typeOfGeneratingProcess = typeOfGeneratingProcess;
}


void varDefProductDefinitionTemplate(int varID, int productDefinitionTemplate)
{
  vartable[varID].productDefinitionTemplate = productDefinitionTemplate;
}


#if  defined  (HAVE_LIBGRIB_API)
void varDefOptGribInt(int varID, long lval, const char *keyword)
{
  int idx = vartable[varID].opt_grib_int_nentries;
  vartable[varID].opt_grib_int_nentries++;
  if ( idx >= MAX_OPT_GRIB_ENTRIES ) Error("Too many optional keyword/integer value pairs!");
  vartable[varID].opt_grib_int_val[idx] = (int) lval;
  vartable[varID].opt_grib_int_keyword[idx] = strdupx(keyword);
}
#endif


#if  defined  (HAVE_LIBGRIB_API)
void varDefOptGribDbl(int varID, double dval, const char *keyword)
{
  int idx = vartable[varID].opt_grib_dbl_nentries;
  vartable[varID].opt_grib_dbl_nentries++;
  if ( idx >= MAX_OPT_GRIB_ENTRIES ) Error("Too many optional keyword/double value pairs!");
  vartable[varID].opt_grib_dbl_val[idx] = dval;
  vartable[varID].opt_grib_dbl_keyword[idx] = strdupx(keyword);
}
#endif


#if  defined  (HAVE_LIBGRIB_API)
int varOptGribNentries(int varID)
{
  int nentries = 0;
  nentries = vartable[varID].opt_grib_int_nentries + vartable[varID].opt_grib_dbl_nentries;
  return (nentries);
}
#endif

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
