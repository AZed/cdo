#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "vlist.h"
#include "zaxis.h"
#include "varscan.h"
#include "namespace.h"
#include "pio_util.h"
#include "resource_handle.h"
#include "vlist_var.h"
#include "vlist_att.h"
#include "pio_rpc.h"

#include "resource_unpack.h"
#include "serialize.h"

#if  defined  (HAVE_LIBGRIB_API)
/* list of additional GRIB2 keywords which are read by the open process */
int    cdiNAdditionalGRIBKeys = 0;
char*  cdiAdditionalGRIBKeys[MAX_OPT_GRIB_ENTRIES];
#endif

extern void zaxisGetIndexList ( int, int * );

static int VLIST_Debug = 0;

static void vlist_initialize(void);

#if  defined  (HAVE_LIBPTHREAD)
#  include <pthread.h>

static pthread_once_t  _vlist_init_thread = PTHREAD_ONCE_INIT;

#  define VLIST_INIT()        \
  pthread_once(&_vlist_init_thread, vlist_initialize)

#else

static int vlistIsInitialized = 0;

#  define VLIST_INIT()               \
  if ( !vlistIsInitialized ) vlist_initialize()
#endif


/* FIXME: implementation incomplete, fix once leaf nodes are complete */
static int
vlist_compare(vlist_t *a, vlist_t *b)
{
  int diff;
  diff = (a->nvars != b->nvars) || (a->ngrids != b->ngrids)
    || (a->nzaxis != b->nzaxis) || (a->instID != b->instID)
    || (a->modelID != b->modelID) || (a->tableID != b->tableID)
    || (a->ntsteps != b->ntsteps);
  return diff;
}

static void
vlistPrintKernel(vlist_t *vlistptr, FILE * fp );
static void
vlist_delete(vlist_t *vlistptr);

static int  vlistGetSizeP ( void * vlistptr, void *context);
static void vlistPackP    ( void * vlistptr, void * buff, int size,
                            int *position, void *context);
static int  vlistTxCode   ( void );

resOps vlist_ops = {
  (valCompareFunc)vlist_compare,
  (valDestroyFunc)vlist_delete,
  (valPrintFunc)vlistPrintKernel
  , vlistGetSizeP,
  vlistPackP,
  vlistTxCode
};


vlist_t *vlist_to_pointer(int code)
{
  VLIST_INIT();
  return reshGetVal(code, &vlist_ops );
}

static
void vlist_init_entry(vlist_t *vlistptr)
{
  vlistptr->self           = reshPut(vlistptr, &vlist_ops);

  vlistptr->nvars          = 0;
  vlistptr->vars           = NULL;
  vlistptr->ngrids         = 0;
  vlistptr->nzaxis         = 0;
  vlistptr->taxisID        = CDI_UNDEFID;
  vlistptr->instID         = cdiDefaultInstID;
  vlistptr->modelID        = cdiDefaultModelID;
  vlistptr->tableID        = cdiDefaultTableID;
  vlistptr->varsAllocated  = 0;
  vlistptr->ntsteps        = CDI_UNDEFID;
  vlistptr->atts.nalloc    = MAX_ATTRIBUTES;
  vlistptr->atts.nelems    = 0;
}

static
vlist_t *vlist_new_entry(void)
{
  vlist_t *vlistptr;

  vlistptr = (vlist_t *)xmalloc(sizeof(vlist_t));

  vlist_init_entry(vlistptr);

  return (vlistptr);
}

static
void vlist_delete_entry(vlist_t *vlistptr)
{
  int idx;

  idx = vlistptr->self;

  reshRemove(idx, &vlist_ops );

  free(vlistptr);

  if ( VLIST_Debug )
    Message("Removed idx %d from vlist list", idx);
}

static
void vlist_initialize(void)
{
  char *env;

  env = getenv("VLIST_DEBUG");
  if ( env ) VLIST_Debug = atoi(env);
#ifndef HAVE_LIBPTHREAD
  vlistIsInitialized = TRUE;
#endif
}

static
void vlist_copy(vlist_t *vlistptr2, vlist_t *vlistptr1)
{
  int vlistID2;

  vlistID2 = vlistptr2->self;
  memcpy(vlistptr2, vlistptr1, sizeof(vlist_t));
  vlistptr2->atts.nelems = 0;
  vlistptr2->self = vlistID2;
}

static
void vlist_check_ptr(const char *caller, vlist_t *vlistptr)
{
  if ( vlistptr == NULL )
    Errorc("vlist undefined!");
}

/*
@Function  vlistCreate
@Title     Create a variable list

@Prototype int vlistCreate(void)

@Example
Here is an example using @func{vlistCreate} to create a variable list
and add a variable with @func{vlistDefVar}.

@Source
#include "cdi.h"
   ...
int vlistID, varID;
   ...
vlistID = vlistCreate();
varID = vlistDefVar(vlistID, gridID, zaxisID, TSTEP_INSTANT);
   ...
streamDefVlist(streamID, vlistID);
   ...
vlistDestroy(vlistID);
   ...
@EndSource
@EndFunction
*/
int vlistCreate(void)
{
  int vlistID = 0;
  vlist_t *vlistptr;

  cdiInitialize();

  VLIST_INIT();

  vlistptr = vlist_new_entry();

  vlistID = vlistptr->self;

  return (vlistID);
}

static void
vlist_delete(vlist_t *vlistptr)
{
  vlist_check_ptr(__func__, vlistptr);

  int vlistID = vlistptr->self;

  vlistDelAtts(vlistID, CDI_GLOBAL);

  int nvars = vlistptr->nvars;

  for (int varID = 0; varID < nvars; varID++ )
    {
      if ( vlistptr->vars[varID].levinfo )  free(vlistptr->vars[varID].levinfo);
      if ( vlistptr->vars[varID].name )     free(vlistptr->vars[varID].name);
      if ( vlistptr->vars[varID].longname ) free(vlistptr->vars[varID].longname);
      if ( vlistptr->vars[varID].stdname )  free(vlistptr->vars[varID].stdname);
      if ( vlistptr->vars[varID].units )    free(vlistptr->vars[varID].units);

      if ( vlistptr->vars[varID].ensdata )  free(vlistptr->vars[varID].ensdata);

#if  defined  (HAVE_LIBGRIB_API)
      int i;
      for (i=0; i<vlistptr->vars[varID].opt_grib_int_nentries; i++) {
	if ( vlistptr->vars[varID].opt_grib_int_keyword[i] )
	  free(vlistptr->vars[varID].opt_grib_int_keyword[i]);
      }
      for (i=0; i<vlistptr->vars[varID].opt_grib_dbl_nentries; i++) {
	if ( vlistptr->vars[varID].opt_grib_dbl_keyword[i] )
	  free(vlistptr->vars[varID].opt_grib_dbl_keyword[i]);
      }
#endif

      vlistDelAtts(vlistID, varID);
    }

  if ( vlistptr->vars ) free(vlistptr->vars);

  vlist_delete_entry(vlistptr);
}


/*
@Function  vlistDestroy
@Title     Destroy a variable list

@Prototype void vlistDestroy(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}

@EndFunction
*/
void vlistDestroy(int vlistID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_delete(vlistptr);
}

/*
@Function  vlistCopy
@Title     Copy a variable list

@Prototype void vlistCopy(int vlistID2, int vlistID1)
@Parameter
    @Item  vlistID2  Target variable list ID
    @Item  vlistID1  Source variable list ID

@Description
The function @func{vlistCopy} copies all entries from vlistID1 to vlistID2.

@EndFunction
*/
void vlistCopy(int vlistID2, int vlistID1)
{
  vlist_t *vlistptr1, *vlistptr2;

  vlistptr1 = vlist_to_pointer(vlistID1);
  vlistptr2 = vlist_to_pointer(vlistID2);

  vlist_check_ptr(__func__, vlistptr1);
  vlist_check_ptr(__func__, vlistptr2);

  vlist_copy(vlistptr2, vlistptr1);

  vlistCopyVarAtts(vlistID1, CDI_GLOBAL, vlistID2, CDI_GLOBAL);

  if ( vlistptr1->vars )
    {
      int nvars = vlistptr1->nvars;
      int nlevs, varID;

      //vlistptr2->varsAllocated = nvars;
      vlistptr2->vars = (var_t *) malloc(vlistptr2->varsAllocated*sizeof(var_t));
      memcpy(vlistptr2->vars, vlistptr1->vars, vlistptr2->varsAllocated*sizeof(var_t));

      for ( varID = 0; varID < nvars; varID++ )
        {
          if ( vlistptr1->vars[varID].name )
            vlistptr2->vars[varID].name = strdupx(vlistptr1->vars[varID].name);

          if ( vlistptr1->vars[varID].longname )
            vlistptr2->vars[varID].longname = strdupx(vlistptr1->vars[varID].longname);

          if ( vlistptr1->vars[varID].stdname )
            vlistptr2->vars[varID].stdname = strdupx(vlistptr1->vars[varID].stdname);

          if ( vlistptr1->vars[varID].units )
            vlistptr2->vars[varID].units = strdupx(vlistptr1->vars[varID].units);

          if ( vlistptr1->vars[varID].ensdata )
            {
              vlistptr2->vars[varID].ensdata = (ensinfo_t *) malloc(sizeof(ensinfo_t));
              memcpy(vlistptr2->vars[varID].ensdata,
                     vlistptr1->vars[varID].ensdata, sizeof(ensinfo_t));
            }
#if  defined  (HAVE_LIBGRIB_API)
          /* ---------------------------------- */
          /* Local change: 2013-01-28, FP (DWD) */
          /* ---------------------------------- */

	  int i;
	  vlistptr2->vars[varID].opt_grib_int_nentries = vlistptr1->vars[varID].opt_grib_int_nentries;
	  for (i=0; i<vlistptr1->vars[varID].opt_grib_int_nentries; i++) {
	    if ( vlistptr1->vars[varID].opt_grib_int_keyword[i] ) {
	      vlistptr2->vars[varID].opt_grib_int_keyword[i] = strdupx(vlistptr1->vars[varID].opt_grib_int_keyword[i]);
	      vlistptr2->vars[varID].opt_grib_int_val[i]     = vlistptr1->vars[varID].opt_grib_int_val[i];
	    }
	  }
	  vlistptr2->vars[varID].opt_grib_dbl_nentries = vlistptr1->vars[varID].opt_grib_dbl_nentries;
	  for (i=0; i<vlistptr1->vars[varID].opt_grib_dbl_nentries; i++) {
	    if ( vlistptr1->vars[varID].opt_grib_dbl_keyword[i] ) {
	      vlistptr2->vars[varID].opt_grib_dbl_keyword[i] = strdupx(vlistptr1->vars[varID].opt_grib_dbl_keyword[i]);
	      vlistptr2->vars[varID].opt_grib_dbl_val[i]     = vlistptr1->vars[varID].opt_grib_dbl_val[i];
	    }
	  }
#endif

	  vlistptr2->vars[varID].atts.nelems = 0;
	  vlistCopyVarAtts(vlistID1, varID, vlistID2, varID);

          nlevs = vlistptr1->vars[varID].nlevs;
          vlistptr2->vars[varID].levinfo = (levinfo_t *) malloc(nlevs*sizeof(levinfo_t));
          memcpy(vlistptr2->vars[varID].levinfo,
                 vlistptr1->vars[varID].levinfo, nlevs*sizeof(levinfo_t));
	}
    }
}

/*
@Function  vlistDuplicate
@Title     Duplicate a variable list

@Prototype int vlistDuplicate(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}

@Description
The function @func{vlistDuplicate} duplicates the variable list from vlistID1.

@Result
@func{vlistDuplicate} returns an identifier to the duplicated variable list.

@EndFunction
*/
int vlistDuplicate(int vlistID)
{
  int vlistIDnew;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  vlistIDnew = vlistCreate();

  vlistCopy(vlistIDnew, vlistID);

  return (vlistIDnew);
}


void vlistClearFlag(int vlistID)
{
  int varID, levID;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  for ( varID = 0; varID < vlistptr->nvars; varID++ )
    {
      vlistptr->vars[varID].flag = FALSE;
      for ( levID = 0; levID < vlistptr->vars[varID].nlevs; levID++ )
        {
          vlistptr->vars[varID].levinfo[levID].flag = FALSE;
        }
    }
}

static
int vlist_generate_zaxis(int vlistID, int zaxistype, int nlevels, double *levels,
                         double *lbounds, double *ubounds, int vctsize, const double *vct)
{
  int zaxisdefined;
  int nzaxis;
  int zaxisID = CDI_UNDEFID;
  int index;
  int zaxisglobdefined = 0;
  int has_bounds = FALSE;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  zaxisdefined = 0;
  nzaxis = vlistptr->nzaxis;

  if ( lbounds && ubounds ) has_bounds = TRUE;

  for ( index = 0; index < nzaxis; ++index )
    {
      zaxisID = vlistptr->zaxisIDs[index];

      if ( zaxisCompare(zaxisID, zaxistype, nlevels, has_bounds, levels, NULL, NULL, 0) == 0 )
        {
          zaxisdefined = 1;
          break;
        }
    }

  if ( ! zaxisdefined )
    {
      nzaxis = zaxisSize();
      if ( nzaxis > 0 )
        {
          int *zaxisIndexList;
          zaxisIndexList = (int *) malloc ( nzaxis * sizeof ( int ));
          zaxisGetIndexList ( nzaxis, zaxisIndexList );
          for ( index = 0; index < nzaxis; ++index )
            {
              zaxisID = zaxisIndexList[index];
              if ( zaxisCompare(zaxisID, zaxistype, nlevels, has_bounds, levels, NULL, NULL, 0) == 0 )
                {
                  zaxisglobdefined = 1;
                  break;
                }
            }
          if ( zaxisIndexList ) free ( zaxisIndexList );
        }
    }

  if ( ! zaxisdefined )
    {
      if ( ! zaxisglobdefined )
	{
	  zaxisID = zaxisCreate(zaxistype, nlevels);
	  zaxisDefLevels(zaxisID, levels);
	  if ( has_bounds )
	    {
	      zaxisDefLbounds(zaxisID, lbounds);
	      zaxisDefUbounds(zaxisID, ubounds);
	    }

	  if ( zaxistype == ZAXIS_HYBRID )
	    {
	      if ( vctsize > 0 )
		zaxisDefVct(zaxisID, vctsize, vct);
	      else
		Warning("VCT missing");
	    }
	}

      nzaxis = vlistptr->nzaxis;
      vlistptr->zaxisIDs[nzaxis] = zaxisID;
      vlistptr->nzaxis++;
    }

  return (zaxisID);
}

/*
@Function  vlistCopyFlag
@Title     Copy some entries of a variable list

@Prototype void vlistCopyFlag(int vlistID2, int vlistID1)
@Parameter
    @Item  vlistID2  Target variable list ID
    @Item  vlistID1  Source variable list ID

@Description
The function @func{vlistCopyFlag} copies all entries with a flag from vlistID1 to vlistID2.

@EndFunction
*/
void vlistCopyFlag(int vlistID2, int vlistID1)
{
  vlist_t *vlistptr1, *vlistptr2;

  vlistptr1 = vlist_to_pointer(vlistID1);
  vlistptr2 = vlist_to_pointer(vlistID2);

  vlist_check_ptr(__func__, vlistptr1);
  vlist_check_ptr(__func__, vlistptr2);

  vlist_copy(vlistptr2, vlistptr1);

  vlistCopyVarAtts(vlistID1, CDI_GLOBAL, vlistID2, CDI_GLOBAL);

  if ( vlistptr1->vars )
    {
      int nvars = vlistptr1->nvars;
      int nvars2 = 0, levID2;
      int nlevs, nlevs2, levID, varID, varID2;
      int gridID, zaxisID;
      int index;

      vlistptr2->ngrids = 0;
      vlistptr2->nzaxis = 0;

      for ( varID = 0; varID < nvars; varID++ )
        if ( vlistptr1->vars[varID].flag ) nvars2++;

      vlistptr2->nvars = nvars2;
      vlistptr2->varsAllocated = nvars2;
      if ( nvars2 > 0 )
        vlistptr2->vars  = (var_t *) malloc(nvars2*sizeof(var_t));
      else
        vlistptr2->vars  = NULL;

      varID2 = 0;
      for ( varID = 0; varID < nvars; varID++ )
	if ( vlistptr1->vars[varID].flag )
	  {
	    vlistptr2->vars[varID2].flag = FALSE;
	    zaxisID = vlistptr1->vars[varID].zaxisID;
	    gridID  = vlistptr1->vars[varID].gridID;

	    memcpy(&vlistptr2->vars[varID2], &vlistptr1->vars[varID], sizeof(var_t));

	    vlistptr1->vars[varID].fvarID = varID2;
	    vlistptr2->vars[varID2].fvarID = varID;

	    vlistptr2->vars[varID2].mvarID = varID2;

	    if ( vlistptr1->vars[varID].name )
	      vlistptr2->vars[varID2].name = strdupx(vlistptr1->vars[varID].name);

	    if ( vlistptr1->vars[varID].longname )
	      vlistptr2->vars[varID2].longname = strdupx(vlistptr1->vars[varID].longname);

	    if ( vlistptr1->vars[varID].stdname )
	      vlistptr2->vars[varID2].stdname = strdupx(vlistptr1->vars[varID].stdname);

	    if ( vlistptr1->vars[varID].units )
	      vlistptr2->vars[varID2].units = strdupx(vlistptr1->vars[varID].units);

            if ( vlistptr1->vars[varID].ensdata )
              {
                vlistptr2->vars[varID2].ensdata = (ensinfo_t *) malloc(sizeof(ensinfo_t));
                memcpy(vlistptr2->vars[varID2].ensdata,
                       vlistptr1->vars[varID].ensdata, sizeof(ensinfo_t));
              }

#if  defined  (HAVE_LIBGRIB_API)
	    /* ---------------------------------- */
	    /* Local change: 2013-01-28, FP (DWD) */
	    /* ---------------------------------- */

	    int i;
	    vlistptr2->vars[varID2].opt_grib_int_nentries = vlistptr1->vars[varID].opt_grib_int_nentries;
	    for (i=0; i<vlistptr1->vars[varID].opt_grib_int_nentries; i++) {
	      if ( vlistptr1->vars[varID].opt_grib_int_keyword[i] ) {
		vlistptr2->vars[varID2].opt_grib_int_keyword[i] = strdupx(vlistptr1->vars[varID].opt_grib_int_keyword[i]);
		vlistptr2->vars[varID2].opt_grib_int_val[i]     = vlistptr1->vars[varID].opt_grib_int_val[i];
	      }
	    }
	    vlistptr2->vars[varID2].opt_grib_dbl_nentries = vlistptr1->vars[varID].opt_grib_dbl_nentries;
	    for (i=0; i<vlistptr1->vars[varID].opt_grib_dbl_nentries; i++) {
	      if ( vlistptr1->vars[varID].opt_grib_dbl_keyword[i] ) {
		vlistptr2->vars[varID2].opt_grib_dbl_keyword[i] = strdupx(vlistptr1->vars[varID].opt_grib_dbl_keyword[i]);
		vlistptr2->vars[varID2].opt_grib_dbl_val[i]     = vlistptr1->vars[varID].opt_grib_dbl_val[i];
	      }
	    }
#endif

	    vlistptr2->vars[varID2].atts.nelems = 0;
	    vlistCopyVarAtts(vlistID1, varID, vlistID2, varID2);

	    nlevs  = vlistptr1->vars[varID].nlevs;
	    nlevs2 = 0;
	    for ( levID = 0; levID < nlevs; levID++ )
	      if ( vlistptr1->vars[varID].levinfo[levID].flag ) nlevs2++;

	    vlistptr2->vars[varID2].levinfo = (levinfo_t *) malloc(nlevs2*sizeof(levinfo_t));

	    if ( nlevs != nlevs2 )
	      {
		int zaxisType;
		int zaxisID2;
		int nvct = 0;
		double *levels;
		double *lbounds = NULL, *ubounds = NULL;
		const double *vct = NULL;
                char ctemp[CDI_MAX_NAME];

		zaxisID = vlistptr1->vars[varID].zaxisID;
		levels = (double *) malloc(nlevs2*sizeof(double));
		levID2 = 0;
		for ( levID = 0; levID < nlevs; ++levID )
		  if ( vlistptr1->vars[varID].levinfo[levID].flag )
		    {
		      vlistptr1->vars[varID].levinfo[levID].flevelID = levID2;
		      vlistptr1->vars[varID].levinfo[levID].mlevelID = levID2;
		      levels[levID2++] = zaxisInqLevel(zaxisID, levID);
		    }

		zaxisType = zaxisInqType(zaxisID);

		if ( zaxisType == ZAXIS_HYBRID )
		  {
		    nvct = zaxisInqVctSize(zaxisID);
		    vct  = zaxisInqVctPtr(zaxisID);
		  }

                if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
                  {
                    double *lbounds1, *ubounds1;
                    lbounds1 = (double *) malloc(nlevs*sizeof(double));
                    ubounds1 = (double *) malloc(nlevs*sizeof(double));

                    zaxisInqLbounds(zaxisID, lbounds1);
                    zaxisInqUbounds(zaxisID, ubounds1);

                    lbounds = (double *) malloc(nlevs2*sizeof(double));
                    ubounds = (double *) malloc(nlevs2*sizeof(double));

                    levID2 = 0;
                    for ( levID = 0; levID < nlevs; ++levID )
                      if ( vlistptr1->vars[varID].levinfo[levID].flag )
                        {
                          lbounds[levID2] = lbounds1[levID];
                          ubounds[levID2] = ubounds1[levID];
                          levID2++;
                        }

                    free(lbounds1);
                    free(ubounds1);
                  }

		zaxisID2 = vlist_generate_zaxis(vlistID2, zaxisType, nlevs2, levels, lbounds, ubounds, nvct, vct);
		free(levels);
                if ( lbounds ) free(lbounds);
                if ( ubounds ) free(ubounds);

                zaxisInqName(zaxisID, ctemp);
                zaxisDefName(zaxisID2, ctemp);
                zaxisInqLongname(zaxisID, ctemp);
                zaxisDefLongname(zaxisID2, ctemp);
                zaxisInqUnits(zaxisID, ctemp);
                zaxisDefUnits(zaxisID2, ctemp);

		zaxisID = zaxisID2;
		vlistptr2->vars[varID2].zaxisID = zaxisID2;
		vlistptr2->vars[varID2].nlevs   = nlevs2;
	      }

	    for ( levID = 0; levID < nlevs2; levID++ )
	      {
		vlistptr2->vars[varID2].levinfo[levID].flag  = FALSE;
		vlistptr2->vars[varID2].levinfo[levID].index = -1;
	      }

	    levID2 = 0;
	    for ( levID = 0; levID < nlevs; levID++ )
	      if ( vlistptr1->vars[varID].levinfo[levID].flag )
		{
		  vlistptr2->vars[varID2].levinfo[levID2].flevelID = levID;
		  vlistptr2->vars[varID2].levinfo[levID2].mlevelID = levID;
		  levID2++;
		}

	    for ( index = 0; index <vlistptr2->ngrids; index++ )
	      if (vlistptr2->gridIDs[index] == gridID ) break;

	    if ( index == vlistptr2->ngrids )
	      {
		vlistptr2->gridIDs[vlistptr2->ngrids++] = gridID;
		if (vlistptr2->ngrids >= MAX_GRIDS_PS )
		  Error("Internal Problem! More than %d grids.", MAX_GRIDS_PS);
	      }

	    for ( index = 0; index < vlistptr2->nzaxis; index++ )
	      if ( vlistptr2->zaxisIDs[index] == zaxisID ) break;

	    if ( index == vlistptr2->nzaxis )
	      {
		vlistptr2->zaxisIDs[vlistptr2->nzaxis++] = zaxisID;
		if (vlistptr2->nzaxis >= MAX_ZAXES_PS )
		  Error("Internal Problem! More than %d zaxis.", MAX_ZAXES_PS);
	      }

	    varID2++;
	  }
    }
}

/*
@Function  vlistCat
@Title     Concatenate two variable lists

@Prototype void vlistCat(int vlistID2, int vlistID1)
@Parameter
    @Item  vlistID2  Target variable list ID
    @Item  vlistID1  Source variable list ID

@Description
Concatenate the variable list vlistID1 at the end of vlistID2.

@EndFunction
*/
void vlistCat(int vlistID2, int vlistID1)
{
  int nvars, nvars1, nvars2;
  int varID, varID2, nlevs;
  int index, gridID, zaxisID;
  vlist_t *vlistptr1, *vlistptr2;

  vlistptr1 = vlist_to_pointer(vlistID1);
  vlistptr2 = vlist_to_pointer(vlistID2);

  vlist_check_ptr(__func__, vlistptr1);
  vlist_check_ptr(__func__, vlistptr2);

  nvars1 = vlistptr1->nvars;
  nvars2 = vlistptr2->nvars;
  nvars = nvars1 + nvars2;
  vlistptr2->nvars = nvars;

  if ( nvars > vlistptr2->varsAllocated )
    {
      vlistptr2->varsAllocated = nvars;
      vlistptr2->vars = (var_t *) realloc(vlistptr2->vars, nvars*sizeof(var_t));
    }
  memcpy(vlistptr2->vars+nvars2, vlistptr1->vars, nvars1*sizeof(var_t));

  for ( varID = 0; varID < nvars1; varID++ )
    {
      varID2 = varID + nvars2;
      vlistptr1->vars[varID].fvarID = varID2;
      vlistptr2->vars[varID2].fvarID = varID;

      vlistptr1->vars[varID].mvarID = varID2;
      vlistptr2->vars[varID2].mvarID = varID;

      if ( vlistptr1->vars[varID].param < 0 )
	{
	  int pnum, pcat, pdis;
	  cdiDecodeParam(vlistptr1->vars[varID].param, &pnum, &pcat, &pdis);
	  pnum = -(varID2+1);
	  vlistptr2->vars[varID2].param = cdiEncodeParam(pnum, pcat, pdis);
	}

      if ( vlistptr1->vars[varID].name )
        vlistptr2->vars[varID2].name = strdupx(vlistptr1->vars[varID].name);

      if ( vlistptr1->vars[varID].longname )
        vlistptr2->vars[varID2].longname = strdupx(vlistptr1->vars[varID].longname);

      if ( vlistptr1->vars[varID].stdname )
        vlistptr2->vars[varID2].stdname = strdupx(vlistptr1->vars[varID].stdname);

      if ( vlistptr1->vars[varID].units )
        vlistptr2->vars[varID2].units = strdupx(vlistptr1->vars[varID].units);

      nlevs = vlistptr1->vars[varID].nlevs;
      vlistptr2->vars[varID2].levinfo = (levinfo_t *) malloc(nlevs*sizeof(levinfo_t));
      memcpy(vlistptr2->vars[varID2].levinfo, vlistptr1->vars[varID].levinfo, nlevs*sizeof(levinfo_t));

      if ( vlistptr1->vars[varID].ensdata )
        {
          vlistptr2->vars[varID2].ensdata = (ensinfo_t *) malloc(sizeof(ensinfo_t));
          memcpy(vlistptr2->vars[varID2].ensdata, vlistptr1->vars[varID].ensdata, sizeof(ensinfo_t));
        }

#if  defined  (HAVE_LIBGRIB_API)
      /* ---------------------------------- */
      /* Local change: 2013-01-28, FP (DWD) */
      /* ---------------------------------- */

      int i;
      vlistptr2->vars[varID2].opt_grib_int_nentries = vlistptr1->vars[varID].opt_grib_int_nentries;
      for (i=0; i<vlistptr1->vars[varID].opt_grib_int_nentries; i++) {
	if ( vlistptr1->vars[varID].opt_grib_int_keyword[i] ) {
	  vlistptr2->vars[varID2].opt_grib_int_keyword[i] = strdupx(vlistptr1->vars[varID].opt_grib_int_keyword[i]);
	  vlistptr2->vars[varID2].opt_grib_int_val[i]     = vlistptr1->vars[varID].opt_grib_int_val[i];
	}
      }
      vlistptr2->vars[varID2].opt_grib_dbl_nentries = vlistptr1->vars[varID].opt_grib_dbl_nentries;
      for (i=0; i<vlistptr1->vars[varID].opt_grib_dbl_nentries; i++) {
	if ( vlistptr1->vars[varID].opt_grib_dbl_keyword[i] ) {
	  vlistptr2->vars[varID2].opt_grib_dbl_keyword[i] = strdupx(vlistptr1->vars[varID].opt_grib_dbl_keyword[i]);
	  vlistptr2->vars[varID2].opt_grib_dbl_val[i]     = vlistptr1->vars[varID].opt_grib_dbl_val[i];
	}
      }
#endif

      vlistptr2->vars[varID2].atts.nelems = 0;
      vlistCopyVarAtts(vlistID1, varID, vlistID2, varID2);

      gridID = vlistptr1->vars[varID].gridID;
      for ( index = 0; index < vlistptr2->ngrids; index++ )
        if ( gridID == vlistptr2->gridIDs[index] ) break;

      if ( index == vlistptr2->ngrids )
	{
	  vlistptr2->gridIDs[vlistptr2->ngrids++] = gridID;
	  if ( vlistptr2->ngrids >= MAX_GRIDS_PS )
	    Error("Internal Problem! More than %d grids.", MAX_GRIDS_PS);
	}

      zaxisID = vlistptr1->vars[varID].zaxisID;
      for ( index = 0; index < vlistptr2->nzaxis; index++ )
        if ( zaxisID == vlistptr2->zaxisIDs[index] ) break;

      if ( index == vlistptr2->nzaxis )
	{
	  vlistptr2->zaxisIDs[vlistptr2->nzaxis++] = zaxisID;
	  if ( vlistptr2->nzaxis >= MAX_ZAXES_PS )
	    Error("Internal Problem! More than %d zaxis.", MAX_ZAXES_PS);
	}
    }
}

/*
@Function  vlistMerge
@Title     Merge two variable lists

@Prototype void vlistMerge(int vlistID2, int vlistID1)
@Parameter
    @Item  vlistID2  Target variable list ID
    @Item  vlistID1  Source variable list ID

@Description
Merge the variable list vlistID1 to the variable list vlistID2.

@EndFunction
*/
void vlistMerge(int vlistID2, int vlistID1)
{
  int nvars1, nvars2;
  int varID = 0, varID2, levID, nlevs, nlevs1, nlevs2;
  int index, zaxisID;
  int zaxisID1, zaxisID2;
  int *lvar;
  double *levels;
  vlist_t *vlistptr1, *vlistptr2;

  vlistptr1 = vlist_to_pointer(vlistID1);
  vlistptr2 = vlist_to_pointer(vlistID2);

  vlist_check_ptr(__func__, vlistptr1);
  vlist_check_ptr(__func__, vlistptr2);

  nvars1 = vlistptr1->nvars;
  nvars2 = vlistptr2->nvars;

  if ( nvars1 == nvars2 )
    {
      for ( varID = 0; varID < nvars2; varID++ )
	{
	  if ( vlistptr1->vars[varID].name && vlistptr2->vars[varID].name )
	    {
	      if ( strcmp(vlistptr1->vars[varID].name,
			  vlistptr2->vars[varID].name) != 0 ) break;
	    }
	  else
	    {
	      if ( vlistptr1->vars[varID].param != vlistptr2->vars[varID].param )
		break;
	    }
	}
    }

  if ( varID == nvars2 ) /* same variables in vlistID1 and vlistID2 */
    {
      for ( varID = 0; varID < nvars2; varID++ )
        {
          vlistptr1->vars[varID].fvarID = varID;
          vlistptr2->vars[varID].fvarID = varID;

          vlistptr1->vars[varID].mvarID = varID;
          vlistptr2->vars[varID].mvarID = varID;

          nlevs1 = vlistptr1->vars[varID].nlevs;
          nlevs2 = vlistptr2->vars[varID].nlevs;

          nlevs = nlevs1 + nlevs2;

          vlistptr2->vars[varID].nlevs = nlevs;
          /*
          fprintf(stderr, "var %d %d %d %d %d\n", varID, nlevs1, nlevs2, nlevs, sizeof(levinfo_t));
          */
          vlistptr2->vars[varID].levinfo =
            (levinfo_t *) realloc(vlistptr2->vars[varID].levinfo, nlevs*sizeof(levinfo_t));

	  memcpy(vlistptr2->vars[varID].levinfo+nlevs2,
		 vlistptr1->vars[varID].levinfo, nlevs1*sizeof(levinfo_t));

	  for ( levID = 0; levID < nlevs1; levID++ )
	    {
	      vlistptr1->vars[varID].levinfo[levID].mlevelID = nlevs2 + levID;
	    }
	}

      lvar = (int *) malloc(nvars2*sizeof(int));
      for ( varID = 0; varID < nvars2; varID++ ) lvar[varID] = FALSE;

      for ( varID = 0; varID < nvars2; varID++ )
        {
          if ( lvar[varID] == TRUE ) continue;

          zaxisID1 = vlistptr1->vars[varID].zaxisID;
          zaxisID2 = vlistptr2->vars[varID].zaxisID;
          /*
          nlevs1 = vlistptr1->vars[varID].nlevs;
          nlevs2 = vlistptr2->vars[varID].nlevs;
          */
          nlevs1 = zaxisInqSize(zaxisID1);
          nlevs2 = zaxisInqSize(zaxisID2);
          /*
          fprintf(stderr, "zaxis %d %d %d %d\n", zaxisID1, zaxisID2, nlevs1, nlevs2);
          */
          nlevs = nlevs1 + nlevs2;

          zaxisID = zaxisDuplicate(zaxisID2);

          zaxisResize(zaxisID, nlevs);

          levels = (double *) malloc(nlevs1*sizeof(double));

          zaxisInqLevels(zaxisID1, levels);
          /*
          for ( levID = 0; levID < nlevs1; levID++ )
            fprintf(stderr, "%d %d %d %d %d %g\n", varID, levID, nlevs1, nlevs2, vlistptr2->vars[varID].nlevs, levels[levID]);
          */
          for ( levID = 0; levID < nlevs1; levID++ )
            zaxisDefLevel(zaxisID, nlevs2+levID, levels[levID]);

          free(levels);

          for ( index = 0; index < vlistptr2->nzaxis; index++ )
            if ( vlistptr2->zaxisIDs[index] == zaxisID2 )
              vlistptr2->zaxisIDs[index] = zaxisID;

          for ( varID2 = 0; varID2 < nvars2; varID2++ )
            if ( lvar[varID2] == FALSE && vlistptr2->vars[varID2].zaxisID == zaxisID2 )
              {
                vlistptr2->vars[varID2].zaxisID = zaxisID;
                lvar[varID2] = TRUE;
              }
        }

      free(lvar);
    }
  else
    {
      vlistCat(vlistID2, vlistID1);
    }
}

/*
@Function  vlistNvars
@Title     Number of variables in a variable list

@Prototype int vlistNvars(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}

@Description
The function @func{vlistNvars} returns the number of variables in the variable list vlistID.

@Result
@func{vlistNvars} returns the number of variables in a variable list.

@EndFunction
*/
int vlistNvars(int vlistID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  return (vlistptr->nvars);
}


int vlistNrecs(int vlistID)
{
  int varID, nrecs = 0;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  for ( varID = 0; varID < vlistptr->nvars; varID++ )
    nrecs +=  vlistptr->vars[varID].nlevs;

  return (nrecs);
}


int vlistNumber(int vlistID)
{
  int varID, number, number2, datatype;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  datatype = vlistptr->vars[0].datatype;
  if (  datatype== DATATYPE_CPX32 || datatype == DATATYPE_CPX64 )
    number = CDI_COMP;
  else
    number = CDI_REAL;

  for ( varID = 1; varID < vlistptr->nvars; varID++ )
    {
      datatype = vlistptr->vars[varID].datatype;
      if ( datatype == DATATYPE_CPX32 || datatype == DATATYPE_CPX64 )
        number2 = CDI_COMP;
      else
        number2 = CDI_REAL;

      if ( number2 != number )
        {
          number = CDI_BOTH;
          break;
        }
    }

  return (number);
}

/*
@Function  vlistNgrids
@Title     Number of grids in a variable list

@Prototype int vlistNgrids(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}

@Description
The function @func{vlistNgrids} returns the number of grids in the variable list vlistID.

@Result
@func{vlistNgrids} returns the number of grids in a variable list.

@EndFunction
*/
int vlistNgrids(int vlistID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  return (vlistptr->ngrids);
}

/*
@Function  vlistNzaxis
@Title     Number of zaxis in a variable list

@Prototype int vlistNzaxis(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}

@Description
The function @func{vlistNzaxis} returns the number of zaxis in the variable list vlistID.

@Result
@func{vlistNzaxis} returns the number of zaxis in a variable list.

@EndFunction
*/
int vlistNzaxis(int vlistID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  return (vlistptr->nzaxis);
}


void vlistDefNtsteps(int vlistID, int nts)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  if ( reshGetStatus ( vlistID, &vlist_ops ) == CLOSED )
    {
      xwarning("%s", "Operation not executed." );
      return;
    }

  vlistptr->ntsteps = nts;
}


int vlistNtsteps(int vlistID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  return (vlistptr->ntsteps);
}

static void
vlistPrintKernel(vlist_t *vlistptr, FILE * fp )
{
  int nvars, flag, index;
  int varID, fvarID, mvarID, flevID, mlevID, levID;
  int param, gridID, zaxisID, tsteptype, nlevs;
  int dtype;

  int iorank;

  char paramstr[32];
  char *name, *longname, *units;
  double level;

  vlist_check_ptr(__func__, vlistptr);

  fprintf ( fp, "#\n# vlistID %d\n#\n", vlistptr->self);

  nvars = vlistptr->nvars;

  fprintf ( fp, "nvars   %d\n", nvars);
  fprintf ( fp, "ngrids  %d\n", vlistptr->ngrids);
  fprintf ( fp, "nzaxis  %d\n", vlistptr->nzaxis);
  fprintf ( fp, "taxisID %d\n", vlistptr->taxisID);
  fprintf ( fp, "instID  %d\n", vlistptr->instID);
  fprintf ( fp, "modelID %d\n", vlistptr->modelID);
  fprintf ( fp, "tableID %d\n", vlistptr->tableID);

  if ( nvars > 0 )
    {
      fprintf(fp, " varID param    gridID zaxisID tsteptype nlevel flag "
              " name     longname iorank\n");
      for ( varID = 0; varID < nvars; varID++ )
        {
          param    = vlistptr->vars[varID].param;
          gridID   = vlistptr->vars[varID].gridID;
          zaxisID  = vlistptr->vars[varID].zaxisID;
	  tsteptype= vlistptr->vars[varID].tsteptype;
          nlevs    = vlistptr->vars[varID].nlevs;
          name     = vlistptr->vars[varID].name;
          longname = vlistptr->vars[varID].longname;
          units    = vlistptr->vars[varID].units;
          flag     = vlistptr->vars[varID].flag;
          iorank   = vlistptr->vars[varID].iorank;

          cdiParamToString(param, paramstr, sizeof(paramstr));
          fprintf(fp, "%6d %-8s %6d %6d %6d %6d %5d %-8s"
                  " %s %6d",
                  varID, paramstr, gridID, zaxisID, tsteptype, nlevs, flag,
                  name ? name : "", longname ? longname : "",
                  iorank);

          if ( units ) fprintf ( fp, "   [%s]", units);
          fprintf ( fp, "\n");
        }

      fprintf(fp, "\n");
      fprintf(fp, " varID  levID fvarID flevID mvarID mlevID  index  dtype  flag  level\n");
      for ( varID = 0; varID < nvars; varID++ )
        {
          nlevs    = vlistptr->vars[varID].nlevs;
          zaxisID  = vlistptr->vars[varID].zaxisID;
          fvarID   = vlistptr->vars[varID].fvarID;
          mvarID   = vlistptr->vars[varID].mvarID;
          dtype    = vlistptr->vars[varID].datatype;
          for ( levID = 0; levID < nlevs; levID++ )
            {
              flevID = vlistptr->vars[varID].levinfo[levID].flevelID;
              mlevID = vlistptr->vars[varID].levinfo[levID].mlevelID;
              index  = vlistptr->vars[varID].levinfo[levID].index;
              flag   = vlistptr->vars[varID].levinfo[levID].flag;
              level  = zaxisInqLevel(zaxisID, levID);
              fprintf(fp, "%6d %6d %6d %6d %6d %6d %6d %6d %5d  %.9g\n",
                      varID, levID, fvarID, flevID, mvarID, mlevID, index,
                      dtype, flag, level);
            }
        }

      fprintf(fp, "\n");
      fprintf(fp, " varID  size iorank\n");
      for ( varID = 0; varID < nvars; varID++ )
        fprintf(fp, "%3d %8d %6d\n", varID,
                vlistptr->vars[varID].nlevs
                * gridInqSize(vlistptr->vars[varID].gridID),
                vlistptr->vars[varID].iorank);
    }
}


void vlistPrint(int vlistID)
{
  vlist_t *vlistptr;

  if ( vlistID == CDI_UNDEFID ) return;

  vlistptr = vlist_to_pointer(vlistID);
  vlist_check_ptr(__func__, vlistptr);
  vlistPrintKernel(vlistptr, stdout);
}

/*
@Function  vlistDefTaxis
@Title     Define the time axis

@Prototype void vlistDefTaxis(int vlistID, int taxisID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate}

@Description
The function @func{vlistDefTaxis} defines the time axis of a variable list.

@EndFunction
*/
void vlistDefTaxis(int vlistID, int taxisID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  if ( reshGetStatus ( vlistID, &vlist_ops ) == CLOSED )
    {
      xwarning("%s", "Operation not executed." );
      return;
    }

  vlistptr->taxisID = taxisID;
}

/*
@Function  vlistInqTaxis
@Title     Get the time axis

@Prototype int vlistInqTaxis(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}

@Description
The function @func{vlistInqTaxis} returns the time axis of a variable list.

@Result
@func{vlistInqTaxis} returns an identifier to the time axis.

@EndFunction
*/
int vlistInqTaxis(int vlistID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  return (vlistptr->taxisID);
}


void  vlistDefTable(int vlistID, int tableID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  if ( reshGetStatus ( vlistID, &vlist_ops ) == CLOSED )
    {
      xwarning("%s", "Operation not executed." );
      return;
    }

  vlistptr->tableID = tableID;
}


int vlistInqTable(int vlistID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  return (vlistptr->tableID);
}


void vlistDefInstitut(int vlistID, int instID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  if ( reshGetStatus ( vlistID, &vlist_ops ) == CLOSED )
    {
      xwarning("%s", "Operation not executed." );

      xdebug("%s", "");
      return;
    }

  vlistptr->instID = instID;
}


int vlistInqInstitut(int vlistID)
{
  int varID, instID;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  instID = vlistptr->instID;

  if ( instID == CDI_UNDEFID )
    {
      instID  = vlistInqVarInstitut(vlistID, 0);

      for ( varID = 1; varID < vlistptr->nvars; varID++ )
        if ( instID != vlistInqVarInstitut(vlistID, varID) )
          {
            instID = CDI_UNDEFID;
            break;
      }
      vlistDefInstitut(vlistID, instID);
    }

  return (instID);
}


void vlistDefModel(int vlistID, int modelID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  if ( reshGetStatus ( vlistID, &vlist_ops ) == CLOSED )
    {
      xwarning("%s", "Operation not executed." );
      return;
    }

  vlistptr->modelID = modelID;
}


int vlistInqModel(int vlistID)
{
  int varID, modelID;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  modelID = vlistptr->modelID;

  if ( modelID == CDI_UNDEFID )
    {
      modelID = vlistInqVarModel(vlistID, 0);

      for ( varID = 1; varID < vlistptr->nvars; varID++ )
        if ( modelID != vlistInqVarModel(vlistID, varID) )
          {
            modelID = CDI_UNDEFID;
            break;
          }

      vlistDefModel(vlistID, modelID);
    }

  return (modelID);
}


int vlistGridsizeMax(int vlistID)
{
  int gridsize, gridsizemax = 0;
  int gridID, index;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  for ( index = 0 ; index < vlistptr->ngrids ; index++ )
    {
      gridID = vlistptr->gridIDs[index];
      gridsize = gridInqSize(gridID);
      if ( gridsize > gridsizemax ) gridsizemax = gridsize;
    }

  return (gridsizemax);
}


int vlistGrid(int vlistID, int index)
{
  int gridID = CDI_UNDEFID;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  if ( index < vlistptr->ngrids && index >= 0 )
    gridID = vlistptr->gridIDs[index];

  return (gridID);
}


int vlistGridIndex(int vlistID, int gridID)
{
  int index;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  for ( index = 0 ; index < vlistptr->ngrids ; index++ )
    if ( gridID == vlistptr->gridIDs[index] ) break;

  if ( index == vlistptr->ngrids ) index = -1;

  return (index);
}


void vlistChangeGridIndex(int vlistID, int index, int gridID)
{
  int gridIDold;
  int varID, nvars;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  if ( reshGetStatus ( vlistID, &vlist_ops ) == CLOSED )
    {
      xwarning("%s", "Operation not executed." );
      return;
    }

  gridIDold = vlistptr->gridIDs[index];
  vlistptr->gridIDs[index] = gridID;

  nvars = vlistptr->nvars;
  for ( varID = 0; varID < nvars; varID++ )
    if ( vlistptr->vars[varID].gridID == gridIDold )
      vlistptr->vars[varID].gridID = gridID;
}


void vlistChangeGrid(int vlistID, int gridID1, int gridID2)
{
  int varID, nvars;
  int index, ngrids;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  if ( reshGetStatus ( vlistID, &vlist_ops ) == CLOSED )
    {
      xwarning("%s", "Operation not executed." );
      return;
    }

  ngrids = vlistptr->ngrids;
  for ( index = 0; index < ngrids; index++ )
    {
      if ( vlistptr->gridIDs[index] == gridID1 )
        {
          vlistptr->gridIDs[index] = gridID2;
          break;
        }
    }

  nvars = vlistptr->nvars;
  for ( varID = 0; varID < nvars; varID++ )
    if ( vlistptr->vars[varID].gridID == gridID1 )
      vlistptr->vars[varID].gridID = gridID2;
}


int vlistZaxis(int vlistID, int index)
{
  int zaxisID = CDI_UNDEFID;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  if ( index < vlistptr->nzaxis && index >= 0 )
    zaxisID = vlistptr->zaxisIDs[index];

  return (zaxisID);
}

int vlistZaxisIndex(int vlistID, int zaxisID)
{
  int index;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  for ( index = 0 ; index < vlistptr->nzaxis ; index++ )
    if ( zaxisID == vlistptr->zaxisIDs[index] ) break;

  if ( index == vlistptr->nzaxis ) index = -1;

  return (index);
}


void vlistChangeZaxisIndex(int vlistID, int index, int zaxisID)
{
  int zaxisIDold;
  int varID, nvars;
  int nlevs, levID;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  if ( reshGetStatus ( vlistID, &vlist_ops ) == CLOSED )
    {
      xwarning("%s", "Operation not executed." );
      return;
    }

  zaxisIDold = vlistptr->zaxisIDs[index];
  vlistptr->zaxisIDs[index] = zaxisID;

  nvars = vlistptr->nvars;
  for ( varID = 0; varID < nvars; varID++ )
    if ( vlistptr->vars[varID].zaxisID == zaxisIDold )
      {
        vlistptr->vars[varID].zaxisID = zaxisID;

        nlevs = zaxisInqSize(zaxisID);
        if ( nlevs != vlistptr->vars[varID].nlevs )
          {
            vlistptr->vars[varID].nlevs   = nlevs;
            vlistptr->vars[varID].levinfo = (levinfo_t *) realloc(vlistptr->vars[varID].levinfo,
                                                                     nlevs*sizeof(levinfo_t));

            for ( levID = 0; levID < nlevs; levID++ )
              {
                vlistptr->vars[varID].levinfo[levID].flevelID = levID;
                vlistptr->vars[varID].levinfo[levID].mlevelID = levID;
                vlistptr->vars[varID].levinfo[levID].index    = -1;
                vlistptr->vars[varID].levinfo[levID].flag     = FALSE;
              }
          }
      }
}


void vlistChangeZaxis(int vlistID, int zaxisID1, int zaxisID2)
{
  int varID, nvars;
  int index, nzaxis;
  int nlevs, levID;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  if ( reshGetStatus ( vlistID, &vlist_ops ) == CLOSED )
    {
      xwarning("%s", "Operation not executed." );
      return;
    }

  nzaxis = vlistptr->nzaxis;
  for ( index = 0; index < nzaxis; index++ )
    {
      if ( vlistptr->zaxisIDs[index] == zaxisID1 )
        {
          vlistptr->zaxisIDs[index] = zaxisID2;
          break;
        }
    }

  nvars = vlistptr->nvars;
  for ( varID = 0; varID < nvars; varID++ )
    if ( vlistptr->vars[varID].zaxisID == zaxisID1 )
      {
        vlistptr->vars[varID].zaxisID = zaxisID2;

        nlevs = zaxisInqSize(zaxisID2);
        if ( nlevs != vlistptr->vars[varID].nlevs )
          {
            vlistptr->vars[varID].nlevs   = nlevs;
            vlistptr->vars[varID].levinfo = (levinfo_t *) realloc(vlistptr->vars[varID].levinfo,
                                                                     nlevs*sizeof(levinfo_t));

            for ( levID = 0; levID < nlevs; levID++ )
              {
                vlistptr->vars[varID].levinfo[levID].flevelID = levID;
                vlistptr->vars[varID].levinfo[levID].mlevelID = levID;
                vlistptr->vars[varID].levinfo[levID].index    = -1;
                vlistptr->vars[varID].levinfo[levID].flag     = FALSE;
              }
          }
      }
}


int vlistHasTime(int vlistID)
{
  int varID;
  int hastime = FALSE;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlist_check_ptr(__func__, vlistptr);

  for ( varID = 0; varID <  vlistptr->nvars; varID++ )
    if ( vlistptr->vars[varID].tsteptype != TSTEP_CONSTANT )
      {
        hastime = TRUE;
        break;
      }

  return (hastime);
}

enum {
  vlist_nints=7,
};

static int
vlistTxCode ( void )
{
  return VLIST;
}


static
int  vlistGetSizeP ( void * vlistptr, void *context)
{
  int txsize, varID;
  vlist_t *p = vlistptr;
  txsize = serializeGetSize(vlist_nints, DATATYPE_INT, context);
  txsize += vlistAttsGetSize(p, CDI_GLOBAL, context);
  for ( varID = 0; varID <  p->nvars; varID++ )
    txsize += vlistVarGetSize(p, varID, context);
  return txsize;
}


static
void vlistPackP ( void * vlistptr, void * buf, int size, int *position,
                  void *context )
{
  int varID, tempbuf[vlist_nints];
  vlist_t *p = vlistptr;
  tempbuf[0] = p->self;
  tempbuf[1] = p->nvars;
  tempbuf[2] = p->ntsteps;
  tempbuf[3] = p->taxisID;
  tempbuf[4] = p->tableID;
  tempbuf[5] = p->instID;
  tempbuf[6] = p->modelID;
  serializePack(tempbuf, vlist_nints, DATATYPE_INT, buf, size, position, context);
  vlistAttsPack(p, CDI_GLOBAL, buf, size, position, context);
  for ( varID = 0; varID < p->nvars; varID++ )
    {
      vlistVarPack(p, varID, buf, size, position, context);
    }
}

void vlistUnpack(char * buf, int size, int *position, int nspTarget, void *context)
{
  int newvlist;
  int varID, tempbuf[vlist_nints];
  serializeUnpack(buf, size, position, tempbuf, vlist_nints, DATATYPE_INT, context);
  newvlist = vlistCreate();
  /* xassert(newvlist == tempbuf[0]); */
  vlistDefTaxis ( newvlist, namespaceAdaptKey ( tempbuf[3], nspTarget ));
  vlistDefTable(newvlist, tempbuf[4]);
  vlistDefInstitut ( newvlist, namespaceAdaptKey ( tempbuf[5], nspTarget ));
  vlistDefModel ( newvlist, namespaceAdaptKey ( tempbuf[6], nspTarget ));
  vlistAttsUnpack(newvlist, CDI_GLOBAL, buf, size, position, context);
  for ( varID = 0; varID < tempbuf[1]; varID++ )
    vlistVarUnpack(newvlist, buf, size, position, nspTarget, context);
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
