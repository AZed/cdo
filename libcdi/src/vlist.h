#ifndef _VLIST_H
#define _VLIST_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stddef.h>  /* size_t */

#ifndef _CDI_LIMITS_H
#  include "cdi_limits.h"
#endif

#define VALIDMISS 1.e+303

/*
 * CDI attribute
 */
typedef struct {
  size_t    xsz;	  /* amount of space at xvalue                      */
  size_t    namesz;       /* size of name                                   */
  char     *name;         /* attribute name                                 */
  int       indtype;	  /* internal data type of xvalue (INT, FLT or TXT) */
  int       exdtype;      /* external data type                             */
                          /* indtype    exdtype                             */
                          /* TXT        TXT                                 */
                          /* INT        INT16, INT32                        */
                          /* FLT        FLT32, FLT64                        */
  size_t    nelems;    	  /* number of elements                             */
  void     *xvalue;       /* the actual data                                */
} cdi_att_t;


typedef struct {
  size_t     nalloc;		/* number allocated >= nelems */
  size_t     nelems;		/* length of the array */
  cdi_att_t  value[MAX_ATTRIBUTES];
} cdi_atts_t;


typedef struct
{
  int      flag;
  int      index;
  int      mlevelID;
  int      flevelID;
}
levinfo_t;

#define DEFAULT_LEVINFO(levID) \
  (levinfo_t){ 0, -1, levID, levID}
/*
#define DEFAULT_LEVINFO(levID) \
  (levinfo_t){ .flag = 0, .index = -1, .flevelID = levID, .mlevelID = levID}
*/
typedef struct
{
  int ens_index;
  int ens_count;
  int forecast_init_type;
}
ensinfo_t;

/* ---------------------------------- */
/* Local change: 2013-01-28, FP (DWD) */
/* ---------------------------------- */

/* Length of optional keyword/value pair list */
#define MAX_OPT_GRIB_ENTRIES 50


typedef struct
{
  int         flag;
  int         isUsed;
  int         mvarID;
  int         fvarID;
  int         param;
  int         gridID;
  int         zaxisID;
  int         tsteptype; /* TSTEP_* */
  int         datatype;  /* DATATYPE_PACKX for GRIB data, else DATATYPE_FLT32 or DATATYPE_FLT64 */
  int         instID;
  int         modelID;
  int         tableID;
  int         timave;
  int         timaccu;
  int         typeOfGeneratingProcess;
  int         chunktype;
  int         xyz;
  int         missvalused; /* TRUE if missval is defined */
  int         lvalidrange;
  char       *name;
  char       *longname;
  char       *stdname;
  char       *units;
  char       *extra;
  double      missval;
  double      scalefactor;
  double      addoffset;
  double      validrange[2];
  levinfo_t  *levinfo;
  int         comptype;     // compression type
  int         complevel;    // compression level
  ensinfo_t  *ensdata;      /* Ensemble information */
  cdi_atts_t  atts;
  int         iorank;
#if  defined  (HAVE_LIBGRIB_API)
  /* ---------------------------------- */
  /* Local change: 2013-01-28, FP (DWD) */
  /* ---------------------------------- */

  /* (Optional) list of keyword/double value pairs */
  int    opt_grib_dbl_nentries;
  char*  opt_grib_dbl_keyword[MAX_OPT_GRIB_ENTRIES];
  double opt_grib_dbl_val[MAX_OPT_GRIB_ENTRIES];
  /* (Optional) list of keyword/integer value pairs */
  int    opt_grib_int_nentries;
  char*  opt_grib_int_keyword[MAX_OPT_GRIB_ENTRIES];
  int    opt_grib_int_val[MAX_OPT_GRIB_ENTRIES];
#endif
}
var_t;


typedef struct
{
  int         self;
  int         nvars;        /* number of variables                */
  int         ngrids;
  int         nzaxis;
  int         ntsteps;
  int         taxisID;
  int         tableID;
  int         instID;
  int         modelID;
  int         varsAllocated;
  int         gridIDs[MAX_GRIDS_PS];
  int         zaxisIDs[MAX_ZAXES_PS];
  var_t      *vars;
  cdi_atts_t  atts;
}
vlist_t;


vlist_t *vlist_to_pointer(int vlistID);
const char *vlistInqVarNamePtr(int vlistID, int varID);
const char *vlistInqVarLongnamePtr(int vlistID, int varID);
const char *vlistInqVarStdnamePtr(int vlistID, int varID);
const char *vlistInqVarUnitsPtr(int vlistID, int varID);
void     vlistDestroyVarName(int vlistID, int varID);
void     vlistDestroyVarLongname(int vlistID, int varID);
void     vlistDestroyVarUnits(int vlistID, int varID);
void     vlistDefVarTsteptype(int vlistID, int varID, int tsteptype);
int      vlistInqVarMissvalUsed(int vlistID, int varID);
int      vlistHasTime(int vlistID);

int      vlistDelAtts(int vlistID, int varID);
int      vlistCopyVarAtts(int vlistID1, int varID_1, int vlistID2, int varID_2);

void     vlistUnpack(char * buffer, int bufferSize, int * pos,
                     int originNamespace, void *context, int force_id);

/*      vlistDefVarValidrange: Define the valid range of a Variable */
void    vlistDefVarValidrange(int vlistID, int varID, const double *validrange);

/*      vlistInqVarValidrange: Get the valid range of a Variable */
int     vlistInqVarValidrange(int vlistID, int varID, double *validrange);

int vlist_att_compare(vlist_t *a, int varIDA, vlist_t *b, int varIDB,
                      int attnum);

#if  defined  (HAVE_LIBGRIB_API)
extern int   cdiNAdditionalGRIBKeys;
extern char* cdiAdditionalGRIBKeys[];
#endif

#endif  /* _VLIST_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
