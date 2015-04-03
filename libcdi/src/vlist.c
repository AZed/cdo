#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "error.h"
#include "vlist.h"
#include "zaxis.h"
#include "varscan.h"
#include "namespace.h"
#include "resource_handle.h"
#include "vlist_var.h"
#include "vlist_att.h"

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


static int
vlist_compare(vlist_t *a, vlist_t *b)
{
  int diff;
  diff = (a->nvars != b->nvars) | (a->ngrids != b->ngrids)
    | (a->nzaxis != b->nzaxis) | (a->instID != b->instID)
    | (a->modelID != b->modelID) | (a->tableID != b->tableID)
    | (a->ntsteps != b->ntsteps) | (a->atts.nelems != b->atts.nelems);
  int nvars = a->nvars;
  for (int varID = 0; varID < nvars; ++varID)
    diff |= vlistVarCompare(a, varID, b, varID);
  size_t natts = a->atts.nelems;
  for (size_t attID = 0; attID < natts; ++attID)
    diff |= vlist_att_compare(a, CDI_GLOBAL, b, CDI_GLOBAL, (int)attID);
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

#if !defined(__cplusplus)
const
#endif
resOps vlistOps = {
  (valCompareFunc)vlist_compare,
  (valDestroyFunc)vlist_delete,
  (valPrintFunc)vlistPrintKernel
  , vlistGetSizeP,
  vlistPackP,
  vlistTxCode
};


vlist_t *vlist_to_pointer(int vlistID)
{
  VLIST_INIT();
  return (vlist_t*) reshGetVal(vlistID, &vlistOps );
}

static
void vlist_init_entry(vlist_t *vlistptr)
{
  vlistptr->locked         = 0;
  vlistptr->self           = CDI_UNDEFID;
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
vlist_t *vlist_new_entry(cdiResH resH)
{
  vlist_t *vlistptr = (vlist_t*) xmalloc(sizeof(vlist_t));
  vlist_init_entry(vlistptr);
  if (resH == CDI_UNDEFID)
    vlistptr->self = reshPut(vlistptr, &vlistOps);
  else
    {
      vlistptr->self = resH;
      reshReplace(resH, vlistptr, &vlistOps);
    }
  return (vlistptr);
}

static
void vlist_delete_entry(vlist_t *vlistptr)
{
  int idx;

  idx = vlistptr->self;

  reshRemove(idx, &vlistOps );

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
  int vlistID2 = vlistptr2->self;
  memcpy(vlistptr2, vlistptr1, sizeof(vlist_t));
  vlistptr2->atts.nelems = 0;
  vlistptr2->self = vlistID2;
}

void vlist_lock(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if ( !vlistptr->locked )
    {
      vlistptr->locked = 1;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


void vlist_unlock(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if ( vlistptr->locked )
    {
      vlistptr->locked = 0;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
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
  cdiInitialize();

  VLIST_INIT();

  vlist_t *vlistptr = vlist_new_entry(CDI_UNDEFID);
  return (vlistptr->self);
}

static void
vlist_delete(vlist_t *vlistptr)
{
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
      for (int i=0; i<vlistptr->vars[varID].opt_grib_int_nentries; i++) {
	if ( vlistptr->vars[varID].opt_grib_int_keyword[i] )
	  free(vlistptr->vars[varID].opt_grib_int_keyword[i]);
      }
      for (int i=0; i<vlistptr->vars[varID].opt_grib_dbl_nentries; i++) {
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
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.

@EndFunction
*/
void vlistDestroy(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if ( vlistptr->locked )
    Warning("Destroying of a locked object (vlistID=%d) failed!", vlistID);
  else
    vlist_delete(vlistptr);
}

/*
@Function  vlistCopy
@Title     Copy a variable list

@Prototype void vlistCopy(int vlistID2, int vlistID1)
@Parameter
    @Item  vlistID2  Target variable list ID.
    @Item  vlistID1  Source variable list ID.

@Description
The function @func{vlistCopy} copies all entries from vlistID1 to vlistID2.

@EndFunction
*/
void vlistCopy(int vlistID2, int vlistID1)
{
  vlist_t *vlistptr1, *vlistptr2;

  vlistptr1 = vlist_to_pointer(vlistID1);
  vlistptr2 = vlist_to_pointer(vlistID2);

  var_t *vlist2vars = vlistptr2->vars;
  vlist_copy(vlistptr2, vlistptr1);

  vlistCopyVarAtts(vlistID1, CDI_GLOBAL, vlistID2, CDI_GLOBAL);

  if ( vlistptr1->vars )
    {
      int nvars = vlistptr1->nvars;

      //vlistptr2->varsAllocated = nvars;
      vlistptr2->vars
        = xrealloc(vlist2vars,
                   (size_t)vlistptr2->varsAllocated * sizeof (var_t));
      memcpy(vlistptr2->vars, vlistptr1->vars,
             (size_t)vlistptr2->varsAllocated * sizeof (var_t));

      for ( int varID = 0; varID < nvars; varID++ )
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

	  vlistptr2->vars[varID].opt_grib_int_nentries = vlistptr1->vars[varID].opt_grib_int_nentries;
	  for (int i=0; i<vlistptr1->vars[varID].opt_grib_int_nentries; i++) {
	    if ( vlistptr1->vars[varID].opt_grib_int_keyword[i] ) {
	      vlistptr2->vars[varID].opt_grib_int_keyword[i] = strdupx(vlistptr1->vars[varID].opt_grib_int_keyword[i]);
	      vlistptr2->vars[varID].opt_grib_int_val[i]     = vlistptr1->vars[varID].opt_grib_int_val[i];
	      vlistptr2->vars[varID].opt_grib_int_update[i]  = TRUE;
	    }
	  }
	  vlistptr2->vars[varID].opt_grib_dbl_nentries = vlistptr1->vars[varID].opt_grib_dbl_nentries;
	  for (int i=0; i<vlistptr1->vars[varID].opt_grib_dbl_nentries; i++) {
	    if ( vlistptr1->vars[varID].opt_grib_dbl_keyword[i] ) {
	      vlistptr2->vars[varID].opt_grib_dbl_keyword[i] = strdupx(vlistptr1->vars[varID].opt_grib_dbl_keyword[i]);
	      vlistptr2->vars[varID].opt_grib_dbl_val[i]     = vlistptr1->vars[varID].opt_grib_dbl_val[i];
	      vlistptr2->vars[varID].opt_grib_dbl_update[i]  = TRUE;
	    }
	  }
#endif

	  vlistptr2->vars[varID].atts.nelems = 0;
	  vlistCopyVarAtts(vlistID1, varID, vlistID2, varID);

          if ( vlistptr1->vars[varID].levinfo )
            {
              size_t nlevs
                = (size_t)zaxisInqSize(vlistptr1->vars[varID].zaxisID);
              vlistptr2->vars[varID].levinfo
                = xmalloc(nlevs * sizeof (levinfo_t));
              memcpy(vlistptr2->vars[varID].levinfo,
                     vlistptr1->vars[varID].levinfo,
                     nlevs * sizeof (levinfo_t));
            }
	}
    }
}

/*
@Function  vlistDuplicate
@Title     Duplicate a variable list

@Prototype int vlistDuplicate(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.

@Description
The function @func{vlistDuplicate} duplicates the variable list from vlistID1.

@Result
@func{vlistDuplicate} returns an identifier to the duplicated variable list.

@EndFunction
*/
int vlistDuplicate(int vlistID)
{
  int vlistIDnew = vlistCreate();
  vlistCopy(vlistIDnew, vlistID);
  return (vlistIDnew);
}


void vlistClearFlag(int vlistID)
{
  int varID, levID;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for ( varID = 0; varID < vlistptr->nvars; varID++ )
    {
      vlistptr->vars[varID].flag = FALSE;
      if ( vlistptr->vars[varID].levinfo )
        {
          int nlevs = zaxisInqSize(vlistptr->vars[varID].zaxisID);
          for ( levID = 0; levID < nlevs; levID++ )
            {
              vlistptr->vars[varID].levinfo[levID].flag = FALSE;
            }
        }
    }
}

static
int vlist_generate_zaxis(int vlistID, int zaxistype, int nlevels, const double *levels,
                         const double *lbounds, const double *ubounds, int vctsize, const double *vct)
{
  int zaxisID = CDI_UNDEFID;
  int zaxisglobdefined = 0;
  int has_bounds = FALSE;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  int zaxisdefined = 0;
  int nzaxis = vlistptr->nzaxis;

  if ( lbounds && ubounds ) has_bounds = TRUE;

  for ( int index = 0; index < nzaxis; ++index )
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
          int *zaxisIndexList = (int *)xmalloc((size_t)nzaxis * sizeof (int));
          reshLock();
          zaxisGetIndexList ( nzaxis, zaxisIndexList );
          for ( int index = 0; index < nzaxis; ++index )
            {
              zaxisID = zaxisIndexList[index];
              if ( zaxisCompare(zaxisID, zaxistype, nlevels, has_bounds, levels, NULL, NULL, 0) == 0 )
                {
                  zaxisglobdefined = 1;
                  break;
                }
            }
          reshUnlock();
          free(zaxisIndexList);
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
    @Item  vlistID2  Target variable list ID.
    @Item  vlistID1  Source variable list ID.

@Description
The function @func{vlistCopyFlag} copies all entries with a flag from vlistID1 to vlistID2.

@EndFunction
*/
void vlistCopyFlag(int vlistID2, int vlistID1)
{
  vlist_t *vlistptr1 = vlist_to_pointer(vlistID1),
    *vlistptr2 = vlist_to_pointer(vlistID2);
  vlist_copy(vlistptr2, vlistptr1);

  vlistCopyVarAtts(vlistID1, CDI_GLOBAL, vlistID2, CDI_GLOBAL);

  if ( vlistptr1->vars )
    {
      int nvars = vlistptr1->nvars;
      int nvars2 = 0;
      int varID2;

      vlistptr2->ngrids = 0;
      vlistptr2->nzaxis = 0;

      for ( int varID = 0; varID < nvars; varID++ )
        nvars2 += (vlistptr1->vars[varID].flag != 0);

      vlistptr2->nvars = nvars2;
      vlistptr2->varsAllocated = nvars2;
      if ( nvars2 > 0 )
        vlistptr2->vars  = (var_t *)xmalloc((size_t)nvars2*sizeof(var_t));
      else
        vlistptr2->vars  = NULL;

      varID2 = 0;
      for ( int varID = 0; varID < nvars; varID++ )
	if ( vlistptr1->vars[varID].flag )
	  {
	    vlistptr2->vars[varID2].flag = FALSE;
	    int zaxisID = vlistptr1->vars[varID].zaxisID;
	    int gridID  = vlistptr1->vars[varID].gridID;

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
                vlistptr2->vars[varID2].ensdata = (ensinfo_t *)xmalloc(sizeof(ensinfo_t));
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
                vlistptr2->vars[varID2].opt_grib_int_update[i]  = TRUE;
	      }
	    }
	    vlistptr2->vars[varID2].opt_grib_dbl_nentries = vlistptr1->vars[varID].opt_grib_dbl_nentries;
	    for (i=0; i<vlistptr1->vars[varID].opt_grib_dbl_nentries; i++) {
	      if ( vlistptr1->vars[varID].opt_grib_dbl_keyword[i] ) {
		vlistptr2->vars[varID2].opt_grib_dbl_keyword[i] = strdupx(vlistptr1->vars[varID].opt_grib_dbl_keyword[i]);
		vlistptr2->vars[varID2].opt_grib_dbl_val[i]     = vlistptr1->vars[varID].opt_grib_dbl_val[i];
                vlistptr2->vars[varID2].opt_grib_dbl_update[i]  = TRUE;
	      }
	    }
#endif

	    vlistptr2->vars[varID2].atts.nelems = 0;
	    vlistCopyVarAtts(vlistID1, varID, vlistID2, varID2);

	    int nlevs  = zaxisInqSize(vlistptr1->vars[varID].zaxisID);
	    int nlevs2 = 0;
            if ( vlistptr1->vars[varID].levinfo )
              for ( int levID = 0; levID < nlevs; levID++ )
                nlevs2 += (vlistptr1->vars[varID].levinfo[levID].flag != 0);

	    vlistptr2->vars[varID2].levinfo = (levinfo_t *)xmalloc((size_t)nlevs2 * sizeof (levinfo_t));

	    if ( nlevs != nlevs2 )
	      {
		int nvct = 0;
		double *lbounds = NULL, *ubounds = NULL;
		const double *vct = NULL;
                char ctemp[CDI_MAX_NAME];

		zaxisID = vlistptr1->vars[varID].zaxisID;
		double *levels = (double *)xmalloc((size_t)nlevs2 * sizeof (double));
                int levID2 = 0;
                if (!vlistptr1->vars[varID].levinfo)
                  cdiVlistCreateVarLevInfo(vlistptr1, varID);
                for ( int levID = 0; levID < nlevs; ++levID )
                  if ( vlistptr1->vars[varID].levinfo[levID].flag )
                    {
                      vlistptr1->vars[varID].levinfo[levID].flevelID = levID2;
                      vlistptr1->vars[varID].levinfo[levID].mlevelID = levID2;
                      levels[levID2++] = zaxisInqLevel(zaxisID, levID);
                    }
		int zaxisType = zaxisInqType(zaxisID);

		if ( zaxisType == ZAXIS_HYBRID )
		  {
		    nvct = zaxisInqVctSize(zaxisID);
		    vct  = zaxisInqVctPtr(zaxisID);
		  }

                if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
                  {
                    lbounds = (double *)xmalloc(2 * (size_t)nlevs2 * sizeof (double));
                    ubounds = lbounds + nlevs2;

                    double *lbounds1 = (double *)xmalloc(2 * (size_t)nlevs * sizeof (double)),
                      *ubounds1 = lbounds1 + nlevs;

                    zaxisInqLbounds(zaxisID, lbounds1);
                    zaxisInqUbounds(zaxisID, ubounds1);

                    int levID2 = 0;
                    for ( int levID = 0; levID < nlevs; ++levID )
                      if ( vlistptr1->vars[varID].levinfo[levID].flag )
                        {
                          lbounds[levID2] = lbounds1[levID];
                          ubounds[levID2] = ubounds1[levID];
                          levID2++;
                        }

                    free(lbounds1);
                  }

		int zaxisID2 = vlist_generate_zaxis(vlistID2, zaxisType, nlevs2, levels, lbounds, ubounds, nvct, vct);
		free(levels);
                free(lbounds);

                zaxisInqName(zaxisID, ctemp);
                zaxisDefName(zaxisID2, ctemp);
                zaxisInqLongname(zaxisID, ctemp);
                zaxisDefLongname(zaxisID2, ctemp);
                zaxisInqUnits(zaxisID, ctemp);
                zaxisDefUnits(zaxisID2, ctemp);

		zaxisID = zaxisID2;
		vlistptr2->vars[varID2].zaxisID = zaxisID2;
	      }

	    for ( int levID = 0; levID < nlevs2; levID++ )
	      {
		vlistptr2->vars[varID2].levinfo[levID].flag  = FALSE;
		vlistptr2->vars[varID2].levinfo[levID].index = -1;
	      }

	    int levID2 = 0;
	    for ( int levID = 0; levID < nlevs; levID++ )
	      if ( vlistptr1->vars[varID].levinfo[levID].flag )
		{
		  vlistptr2->vars[varID2].levinfo[levID2].flevelID = levID;
		  vlistptr2->vars[varID2].levinfo[levID2].mlevelID = levID;
		  levID2++;
		}

            vlistAdd2GridIDs(vlistptr2, gridID);
            vlistAdd2ZaxisIDs(vlistptr2, zaxisID);

	    varID2++;
	  }
    }
}

/*
@Function  vlistCat
@Title     Concatenate two variable lists

@Prototype void vlistCat(int vlistID2, int vlistID1)
@Parameter
    @Item  vlistID2  Target variable list ID.
    @Item  vlistID1  Source variable list ID.

@Description
Concatenate the variable list vlistID1 at the end of vlistID2.

@EndFunction
*/
void vlistCat(int vlistID2, int vlistID1)
{
  vlist_t *vlistptr1 = vlist_to_pointer(vlistID1),
    *vlistptr2 = vlist_to_pointer(vlistID2);

  int nvars1 = vlistptr1->nvars;
  int nvars2 = vlistptr2->nvars;
  int nvars = nvars1 + nvars2;
  vlistptr2->nvars = nvars;

  if ( nvars > vlistptr2->varsAllocated )
    {
      vlistptr2->varsAllocated = nvars;
      vlistptr2->vars = xrealloc(vlistptr2->vars,
                                 (size_t)nvars * sizeof (var_t));
    }
  memcpy(vlistptr2->vars+nvars2, vlistptr1->vars,
         (size_t)nvars1 * sizeof (var_t));

  for (int varID = 0; varID < nvars1; varID++ )
    {
      int varID2 = varID + nvars2;
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

      int nlevs = zaxisInqSize(vlistptr1->vars[varID].zaxisID);
      if (vlistptr1->vars[varID].levinfo)
        {
          vlistptr2->vars[varID2].levinfo
            = (levinfo_t *)xmalloc((size_t)nlevs * sizeof (levinfo_t));
          memcpy(vlistptr2->vars[varID2].levinfo,
                 vlistptr1->vars[varID].levinfo,
                 (size_t)nlevs * sizeof (levinfo_t));
        }

      if ( vlistptr1->vars[varID].ensdata )
        {
          vlistptr2->vars[varID2].ensdata = (ensinfo_t *) malloc(sizeof(ensinfo_t));
          memcpy(vlistptr2->vars[varID2].ensdata, vlistptr1->vars[varID].ensdata, sizeof(ensinfo_t));
        }

#if  defined  (HAVE_LIBGRIB_API)
      /* ---------------------------------- */
      /* Local change: 2013-01-28, FP (DWD) */
      /* ---------------------------------- */

      vlistptr2->vars[varID2].opt_grib_int_nentries = vlistptr1->vars[varID].opt_grib_int_nentries;
      int n = vlistptr1->vars[varID].opt_grib_int_nentries;
      for (int i = 0; i < n; ++i) {
	if ( vlistptr1->vars[varID].opt_grib_int_keyword[i] ) {
	  vlistptr2->vars[varID2].opt_grib_int_keyword[i] = strdupx(vlistptr1->vars[varID].opt_grib_int_keyword[i]);
	  vlistptr2->vars[varID2].opt_grib_int_val[i]     = vlistptr1->vars[varID].opt_grib_int_val[i];
          vlistptr2->vars[varID2].opt_grib_int_update[i]  = TRUE;
	}
      }
      vlistptr2->vars[varID2].opt_grib_dbl_nentries = vlistptr1->vars[varID].opt_grib_dbl_nentries;
      n = vlistptr1->vars[varID].opt_grib_dbl_nentries;
      for (int i = 0; i < n; i++) {
	if ( vlistptr1->vars[varID].opt_grib_dbl_keyword[i] ) {
	  vlistptr2->vars[varID2].opt_grib_dbl_keyword[i] = strdupx(vlistptr1->vars[varID].opt_grib_dbl_keyword[i]);
	  vlistptr2->vars[varID2].opt_grib_dbl_val[i]     = vlistptr1->vars[varID].opt_grib_dbl_val[i];
          vlistptr2->vars[varID2].opt_grib_dbl_update[i]  = TRUE;
	}
      }
#endif

      vlistptr2->vars[varID2].atts.nelems = 0;
      vlistCopyVarAtts(vlistID1, varID, vlistID2, varID2);

      vlistAdd2GridIDs(vlistptr2, vlistptr1->vars[varID].gridID);
      vlistAdd2ZaxisIDs(vlistptr2, vlistptr1->vars[varID].zaxisID);
    }
}

/*
@Function  vlistMerge
@Title     Merge two variable lists

@Prototype void vlistMerge(int vlistID2, int vlistID1)
@Parameter
    @Item  vlistID2  Target variable list ID.
    @Item  vlistID1  Source variable list ID.

@Description
Merge the variable list vlistID1 to the variable list vlistID2.

@EndFunction
*/
void vlistMerge(int vlistID2, int vlistID1)
{
  int varID = 0;
  vlist_t *vlistptr1 = vlist_to_pointer(vlistID1),
    *vlistptr2 = vlist_to_pointer(vlistID2);

  int nvars1 = vlistptr1->nvars;
  int nvars2 = vlistptr2->nvars;

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

          int nlevs1 = zaxisInqSize(vlistptr1->vars[varID].zaxisID);
          int nlevs2 = zaxisInqSize(vlistptr2->vars[varID].zaxisID);

          int nlevs = nlevs1 + nlevs2;

          /*
          fprintf(stderr, "var %d %d %d %d %d\n", varID, nlevs1, nlevs2, nlevs, sizeof(levinfo_t));
          */
          if (vlistptr1->vars[varID].levinfo)
            {
              vlistptr2->vars[varID].levinfo =
                (levinfo_t*)xrealloc(vlistptr2->vars[varID].levinfo,
                                     (size_t)nlevs * sizeof(levinfo_t));

              memcpy(vlistptr2->vars[varID].levinfo+nlevs2,
                     vlistptr1->vars[varID].levinfo,
                     (size_t)nlevs1 * sizeof (levinfo_t));
            }
          else
            cdiVlistCreateVarLevInfo(vlistptr1, varID);
	  for ( int levID = 0; levID < nlevs1; levID++ )
	    {
	      vlistptr1->vars[varID].levinfo[levID].mlevelID = nlevs2 + levID;
	    }
	}

      int *lvar = (int *)xcalloc((size_t)nvars2, sizeof(int));

      for ( varID = 0; varID < nvars2; varID++ )
        {
          if ( lvar[varID] == TRUE ) continue;

          int zaxisID1 = vlistptr1->vars[varID].zaxisID;
          int zaxisID2 = vlistptr2->vars[varID].zaxisID;
          /*
          nlevs1 = zaxisInqSize(vlistptr1->vars[varID].zaxisID);
          nlevs2 = zaxisInqSize(vlistptr2->vars[varID].zaxisID);
          */
          int nlevs1 = zaxisInqSize(zaxisID1);
          int nlevs2 = zaxisInqSize(zaxisID2);
          /*
          fprintf(stderr, "zaxis %d %d %d %d\n", zaxisID1, zaxisID2, nlevs1, nlevs2);
          */
          int nlevs = nlevs1 + nlevs2;

          int zaxisID = zaxisDuplicate(zaxisID2);

          zaxisResize(zaxisID, nlevs);

          double *levels = (double *)xmalloc((size_t)nlevs1 * sizeof(double));

          zaxisInqLevels(zaxisID1, levels);
          /*
          for ( levID = 0; levID < nlevs1; levID++ )
            fprintf(stderr, "%d %d %d %d %d %g\n", varID, levID, nlevs1, nlevs2, vlistptr2->vars[varID].nlevs, levels[levID]);
          */
          for ( int levID = 0; levID < nlevs1; levID++ )
            zaxisDefLevel(zaxisID, nlevs2+levID, levels[levID]);

          free(levels);

          for ( int index = 0; index < vlistptr2->nzaxis; index++ )
            if ( vlistptr2->zaxisIDs[index] == zaxisID2 )
              vlistptr2->zaxisIDs[index] = zaxisID;

          for ( int varID2 = 0; varID2 < nvars2; varID2++ )
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
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.

@Description
The function @func{vlistNvars} returns the number of variables in the variable list vlistID.

@Result
@func{vlistNvars} returns the number of variables in a variable list.

@EndFunction
*/
int vlistNvars(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return (vlistptr->nvars);
}


int vlistNrecs(int vlistID)
{
  int nrecs = 0;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for ( int varID = 0; varID < vlistptr->nvars; varID++ )
    nrecs +=  zaxisInqSize(vlistptr->vars[varID].zaxisID);

  return (nrecs);
}


int vlistNumber(int vlistID)
{
  int number, number2, datatype;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  datatype = vlistptr->vars[0].datatype;
  if (  datatype== DATATYPE_CPX32 || datatype == DATATYPE_CPX64 )
    number = CDI_COMP;
  else
    number = CDI_REAL;

  for ( int varID = 1; varID < vlistptr->nvars; varID++ )
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
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.

@Description
The function @func{vlistNgrids} returns the number of grids in the variable list vlistID.

@Result
@func{vlistNgrids} returns the number of grids in a variable list.

@EndFunction
*/
int vlistNgrids(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->ngrids);
}

/*
@Function  vlistNzaxis
@Title     Number of zaxis in a variable list

@Prototype int vlistNzaxis(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.

@Description
The function @func{vlistNzaxis} returns the number of zaxis in the variable list vlistID.

@Result
@func{vlistNzaxis} returns the number of zaxis in a variable list.

@EndFunction
*/
int vlistNzaxis(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->nzaxis);
}


void vlistDefNtsteps(int vlistID, int nts)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->ntsteps != nts)
    {
      vlistptr->ntsteps = nts;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

// This function is used in CDO!
int vlistNtsteps(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  return (int)vlistptr->ntsteps;
}

static void
vlistPrintKernel(vlist_t *vlistptr, FILE * fp )
{
  char paramstr[32];

  fprintf ( fp, "#\n# vlistID %d\n#\n", vlistptr->self);

  int nvars = vlistptr->nvars;

  fprintf(fp, "nvars   %d\n"
          "ngrids  %d\n"
          "nzaxis  %d\n"
          "taxisID %d\n"
          "instID  %d\n"
          "modelID %d\n"
          "tableID %d\n",
          nvars, vlistptr->ngrids, vlistptr->nzaxis, vlistptr->taxisID,
          vlistptr->instID, vlistptr->modelID, vlistptr->tableID);

  if ( nvars > 0 )
    {
      fprintf(fp, " varID param    gridID zaxisID tsteptype flag "
              " name     longname iorank\n");
      for ( int varID = 0; varID < nvars; varID++ )
        {
          int param = vlistptr->vars[varID].param;
          int gridID = vlistptr->vars[varID].gridID;
          int zaxisID = vlistptr->vars[varID].zaxisID;
	  int tsteptype = vlistptr->vars[varID].tsteptype;
          const char *name = vlistptr->vars[varID].name;
          const char *longname = vlistptr->vars[varID].longname;
          const char *units = vlistptr->vars[varID].units;
          int flag = vlistptr->vars[varID].flag;
          int iorank = vlistptr->vars[varID].iorank;

          cdiParamToString(param, paramstr, sizeof(paramstr));
          fprintf(fp, "%6d %-8s %6d %6d %6d %5d %-8s"
                  " %s %6d",
                  varID, paramstr, gridID, zaxisID, tsteptype, flag,
                  name ? name : "", longname ? longname : "",
                  iorank);

          if ( units ) fprintf(fp, "   [%s]", units);
          fputs("\n", fp);
        }

      fputs("\n"
            " varID  levID fvarID flevID mvarID mlevID  index  dtype  flag  level\n", fp);
      for ( int varID = 0; varID < nvars; varID++ )
        {
          int zaxisID = vlistptr->vars[varID].zaxisID;
          int nlevs = zaxisInqSize(zaxisID);
          int fvarID = vlistptr->vars[varID].fvarID;
          int mvarID = vlistptr->vars[varID].mvarID;
          int dtype    = vlistptr->vars[varID].datatype;
          for ( int levID = 0; levID < nlevs; levID++ )
            {
              levinfo_t li;
              if (vlistptr->vars[varID].levinfo)
                li = vlistptr->vars[varID].levinfo[levID];
              else
                li = DEFAULT_LEVINFO(levID);
              int flevID = li.flevelID;
              int mlevID = li.mlevelID;
              int index  = li.index;
              int flag   = li.flag;
              double level  = zaxisInqLevel(zaxisID, levID);

              fprintf(fp, "%6d %6d %6d %6d %6d %6d %6d %6d %5d  %.9g\n",
                      varID, levID, fvarID, flevID, mvarID, mlevID, index,
                      dtype, flag, level);
            }
        }

      fputs("\n"
            " varID  size iorank\n", fp);
      for ( int varID = 0; varID < nvars; varID++ )
        fprintf(fp, "%3d %8d %6d\n", varID,
                zaxisInqSize(vlistptr->vars[varID].zaxisID)
                * gridInqSize(vlistptr->vars[varID].gridID),
                vlistptr->vars[varID].iorank);
    }
}


void vlistPrint(int vlistID)
{
  if ( vlistID == CDI_UNDEFID ) return;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  vlistPrintKernel(vlistptr, stdout);
}

/*
@Function  vlistDefTaxis
@Title     Define the time axis

@Prototype void vlistDefTaxis(int vlistID, int taxisID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate}.

@Description
The function @func{vlistDefTaxis} defines the time axis of a variable list.

@EndFunction
*/
void vlistDefTaxis(int vlistID, int taxisID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->taxisID != taxisID)
    {
      vlistptr->taxisID = taxisID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  vlistInqTaxis
@Title     Get the time axis

@Prototype int vlistInqTaxis(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.

@Description
The function @func{vlistInqTaxis} returns the time axis of a variable list.

@Result
@func{vlistInqTaxis} returns an identifier to the time axis.

@EndFunction
*/
int vlistInqTaxis(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->taxisID);
}


void vlistDefTable(int vlistID, int tableID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->tableID != tableID)
    {
      vlistptr->tableID = tableID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqTable(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->tableID);
}


void vlistDefInstitut(int vlistID, int instID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->instID != instID)
    {
      vlistptr->instID = instID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqInstitut(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int instID = vlistptr->instID;

  if ( instID == CDI_UNDEFID )
    {
      instID  = vlistInqVarInstitut(vlistID, 0);

      for ( int varID = 1; varID < vlistptr->nvars; varID++ )
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
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->modelID != modelID)
    {
      vlistptr->modelID = modelID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqModel(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int modelID = vlistptr->modelID;

  if ( modelID == CDI_UNDEFID )
    {
      modelID = vlistInqVarModel(vlistID, 0);

      for ( int varID = 1; varID < vlistptr->nvars; varID++ )
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
  int gridsizemax = 0;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for ( int index = 0 ; index < vlistptr->ngrids ; index++ )
    {
      int gridID = vlistptr->gridIDs[index];
      int gridsize = gridInqSize(gridID);
      if ( gridsize > gridsizemax ) gridsizemax = gridsize;
    }

  return (gridsizemax);
}


int vlistGrid(int vlistID, int index)
{
  int gridID = CDI_UNDEFID;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if ( index < vlistptr->ngrids && index >= 0 )
    gridID = vlistptr->gridIDs[index];

  return (gridID);
}


int vlistGridIndex(int vlistID, int gridID)
{
  int index;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for ( index = 0 ; index < vlistptr->ngrids ; index++ )
    if ( gridID == vlistptr->gridIDs[index] ) break;

  if ( index == vlistptr->ngrids ) index = -1;

  return (index);
}


void vlistChangeGridIndex(int vlistID, int index, int gridID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int gridIDold = vlistptr->gridIDs[index];
  if (gridIDold != gridID)
    {
      vlistptr->gridIDs[index] = gridID;

      int nvars = vlistptr->nvars;
      for ( int varID = 0; varID < nvars; varID++ )
        if ( vlistptr->vars[varID].gridID == gridIDold )
          vlistptr->vars[varID].gridID = gridID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


void vlistChangeGrid(int vlistID, int gridID1, int gridID2)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (gridID1 != gridID2)
    {
      int ngrids = vlistptr->ngrids;
      for ( int index = 0; index < ngrids; index++ )
        {
          if ( vlistptr->gridIDs[index] == gridID1 )
            {
              vlistptr->gridIDs[index] = gridID2;
              break;
            }
        }
      int nvars = vlistptr->nvars;
      for ( int varID = 0; varID < nvars; varID++ )
        if ( vlistptr->vars[varID].gridID == gridID1 )
          vlistptr->vars[varID].gridID = gridID2;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistZaxis(int vlistID, int index)
{
  int zaxisID = CDI_UNDEFID;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if ( index < vlistptr->nzaxis && index >= 0 )
    zaxisID = vlistptr->zaxisIDs[index];

  return (zaxisID);
}

int vlistZaxisIndex(int vlistID, int zaxisID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int index;
  for ( index = 0 ; index < vlistptr->nzaxis ; index++ )
    if ( zaxisID == vlistptr->zaxisIDs[index] ) break;

  if ( index == vlistptr->nzaxis ) index = -1;

  return (index);
}


void vlistChangeZaxisIndex(int vlistID, int index, int zaxisID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int zaxisIDold = vlistptr->zaxisIDs[index];
  if (zaxisIDold != zaxisID)
    {
      vlistptr->zaxisIDs[index] = zaxisID;

      int nlevs = zaxisInqSize(zaxisID),
        nlevsOld = zaxisInqSize(zaxisIDold);
      int nvars = vlistptr->nvars;
      for ( int varID = 0; varID < nvars; varID++ )
        if ( vlistptr->vars[varID].zaxisID == zaxisIDold )
          {
            vlistptr->vars[varID].zaxisID = zaxisID;
            if ( vlistptr->vars[varID].levinfo && nlevs != nlevsOld )
              {
                vlistptr->vars[varID].levinfo = (levinfo_t *)xrealloc(vlistptr->vars[varID].levinfo, (size_t)nlevs * sizeof (levinfo_t));

                for ( int levID = 0; levID < nlevs; levID++ )
                  vlistptr->vars[varID].levinfo[levID] = DEFAULT_LEVINFO(levID);
              }
          }
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


void vlistChangeZaxis(int vlistID, int zaxisID1, int zaxisID2)
{
  int nlevs1 = zaxisInqSize(zaxisID1), nlevs2 = zaxisInqSize(zaxisID2);
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int nzaxis = vlistptr->nzaxis;
  for ( int index = 0; index < nzaxis; index++ )
    {
      if ( vlistptr->zaxisIDs[index] == zaxisID1 )
        {
          vlistptr->zaxisIDs[index] = zaxisID2;
          break;
        }
    }

  int nvars = vlistptr->nvars;
  for ( int varID = 0; varID < nvars; varID++ )
    if ( vlistptr->vars[varID].zaxisID == zaxisID1 )
      {
        vlistptr->vars[varID].zaxisID = zaxisID2;

        if ( vlistptr->vars[varID].levinfo && nlevs2 != nlevs1 )
          {
            vlistptr->vars[varID].levinfo
              = (levinfo_t *)xrealloc(vlistptr->vars[varID].levinfo,
                                      (size_t)nlevs2 * sizeof(levinfo_t));

            for ( int levID = 0; levID < nlevs2; levID++ )
              vlistptr->vars[varID].levinfo[levID] = DEFAULT_LEVINFO(levID);
          }
      }
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}


int vlistHasTime(int vlistID)
{
  int hastime = FALSE;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for ( int varID = 0; varID <  vlistptr->nvars; varID++ )
    if ( vlistptr->vars[varID].tsteptype != TSTEP_CONSTANT )
      {
        hastime = TRUE;
        break;
      }

  return (hastime);
}

enum {
  vlist_nints=6,
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
  vlist_t *p = (vlist_t*) vlistptr;
  txsize = serializeGetSize(vlist_nints, DATATYPE_INT, context);
  txsize += serializeGetSize(1, DATATYPE_LONG, context);
  txsize += vlistAttsGetSize(p, CDI_GLOBAL, context);
  for ( varID = 0; varID <  p->nvars; varID++ )
    txsize += vlistVarGetPackSize(p, varID, context);
  return txsize;
}


static
void vlistPackP ( void * vlistptr, void * buf, int size, int *position,
                  void *context )
{
  int varID, tempbuf[vlist_nints];
  vlist_t *p = (vlist_t*) vlistptr;
  tempbuf[0] = p->self;
  tempbuf[1] = p->nvars;
  tempbuf[2] = p->taxisID;
  tempbuf[3] = p->tableID;
  tempbuf[4] = p->instID;
  tempbuf[5] = p->modelID;
  serializePack(tempbuf, vlist_nints, DATATYPE_INT, buf, size, position, context);
  serializePack(&p->ntsteps, 1, DATATYPE_LONG, buf, size, position, context);

  vlistAttsPack(p, CDI_GLOBAL, buf, size, position, context);
  for ( varID = 0; varID < p->nvars; varID++ )
    {
      vlistVarPack(p, varID, buf, size, position, context);
    }
}

void vlistUnpack(char * buf, int size, int *position, int originNamespace,
                 void *context, int force_id)
{
  int tempbuf[vlist_nints];
  serializeUnpack(buf, size, position, tempbuf, vlist_nints, DATATYPE_INT, context);
  int nvars = tempbuf[1];
  int targetID = namespaceAdaptKey(tempbuf[0], originNamespace);
  vlist_t *p = vlist_new_entry(force_id?targetID:CDI_UNDEFID);
  xassert(!force_id || p->self == targetID);
  if (!force_id)
    targetID = p->self;
  p->taxisID = namespaceAdaptKey(tempbuf[2], originNamespace);
  p->tableID = tempbuf[3];
  p->instID = namespaceAdaptKey(tempbuf[4], originNamespace);
  p->modelID = namespaceAdaptKey(tempbuf[5], originNamespace);
  serializeUnpack(buf, size, position, &p->ntsteps, 1, DATATYPE_LONG, context);
  vlistAttsUnpack(targetID, CDI_GLOBAL, buf, size, position, context);
  for (int varID = 0; varID < nvars; varID++ )
    vlistVarUnpack(targetID, buf, size, position, originNamespace, context);
}


void vlist_check_contents(int vlistID)
{
  int index, nzaxis, zaxisID;

  nzaxis = vlistNzaxis(vlistID);

  for ( index = 0; index < nzaxis; index++ )
    {
      zaxisID = vlistZaxis(vlistID, index);
      if ( zaxisInqType(zaxisID) == ZAXIS_GENERIC )
	cdiCheckZaxis(zaxisID);
    }
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
