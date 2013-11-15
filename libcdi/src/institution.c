#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "resource_handle.h"
#include "pio_util.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "namespace.h"
#include "serialize.h"

#undef  UNDEFID
#define UNDEFID  -1

int ECMWF  = UNDEFID;
int MPIMET = UNDEFID;
int DWD    = UNDEFID;
int MCH    = UNDEFID;

typedef struct
{
  int    self;
  int    used;
  int    center;
  int    subcenter;
  char  *name;
  char  *longname;
}
institute_t;


static int    instituteCompareP ( void * instituteptr1, void * instituteptr2 );
static void   instituteDestroyP ( void * instituteptr );
static void   institutePrintP   ( void * instituteptr, FILE * fp );
static int    instituteGetSizeP ( void * instituteptr, void *context );
static void   institutePackP    ( void * instituteptr, void *buf, int size, int *position, void *context );
static int    instituteTxCode   ( void );

resOps instituteOps = { instituteCompareP, instituteDestroyP, institutePrintP
                        ,instituteGetSizeP, institutePackP, instituteTxCode
 };

static int * instituteInitializedNsp;

static
void instituteDefaultValue ( institute_t * instituteptr )
{
  instituteptr->self       = UNDEFID;
  instituteptr->used       = 0;
  instituteptr->center     = UNDEFID;
  instituteptr->subcenter  = UNDEFID;
  instituteptr->name       = NULL;
  instituteptr->longname   = NULL;
}

static
institute_t * instituteNewEntry ( void )
{
  institute_t *instituteptr;

  instituteptr = ( institute_t * ) xmalloc ( sizeof ( institute_t ));

  instituteDefaultValue ( instituteptr );
  instituteptr->self = reshPut (( void * ) instituteptr, &instituteOps );
  instituteptr->used = 1;

  return  instituteptr;
}

static
void instituteDefaultEntries ( void )
{
  cdiResH resH[12];
  int i;

  resH[0]  = ECMWF   = institutDef( 98,   0, "ECMWF",     "European Centre for Medium-Range Weather Forecasts");
  resH[1]  = MPIMET  = institutDef( 98, 232, "MPIMET",    "Max-Planck-Institute for Meteorology");
  resH[2]  =           institutDef( 98, 255, "MPIMET",    "Max-Planck-Institute for Meteorology");
  resH[3]  =           institutDef( 98, 232, "MPIMET",    "Max-Planck Institute for Meteorology");
  resH[4]  =           institutDef( 78, 255, "DWD",       "Deutscher Wetterdienst");
  resH[5]  = MCH     = institutDef(215, 255, "MCH",       "MeteoSwiss");
  resH[6]  =           institutDef(  7,   0, "NCEP",      "National Centers for Environmental Prediction");
  resH[7]  =           institutDef(  7,   1, "NCEP",      "National Centers for Environmental Prediction");
  resH[8]  =           institutDef( 60,   0, "NCAR",      "National Center for Atmospheric Research");
  resH[9]  =           institutDef( 74,   0, "METOFFICE", "U.K. Met Office");
  resH[10] =           institutDef( 97,   0, "ESA",       "European Space Agency ");
  resH[11] =           institutDef( 99,   0, "KNMI",      "Royal Netherlands Meteorological Institute");
  /*     (void) institutDef(  0,   0, "IPSL", "IPSL (Institut Pierre Simon Laplace, Paris, France)"); */

  for ( i = 0; i < 12 ; i++ )
    reshSetStatus(resH[i], &instituteOps, SUSPENDED);
}

static
void instituteFinalize ( void )
{
  free ( instituteInitializedNsp );
}

static
void instituteInit (void)
{
  static int instituteInitialized = 0;
  int nsp, nspc;

  nspc = namespaceGetNumber ();

  if ( !instituteInitialized )
    {
      instituteInitialized = 1;
      instituteInitializedNsp = xcalloc ( 1, nspc * sizeof ( int ));
      atexit ( instituteFinalize );
    }

  nsp = namespaceGetActive ();
  if ( instituteInitializedNsp[nsp] ) return;
  instituteInitializedNsp[nsp] = 1;

  instituteDefaultEntries();
}


int instituteCount ( void )
{
  return reshCountType ( &instituteOps );
}


int instituteCompareKernel ( institute_t *  ip1, institute_t * ip2 )
{
  int differ = 0;
  size_t len;

  if ( ip1->name )
    {
      if ( ip1->center    > 0 && ip2->center    != ip1->center )    differ = 1;
      if ( ip1->subcenter > 0 && ip2->subcenter != ip1->subcenter ) differ = 1;

      if ( !differ )
        {
          if ( ip2->name )
            {
              len = strlen(ip2->name);
              if ( memcmp(ip2->name, ip1->name, len)) differ = 1;
            }
        }
    }
  else if ( ip1->longname )
    {
      if ( ip2->longname )
        {
          len = strlen(ip2->longname);
          if ( memcmp(ip2->longname, ip1->longname, len)) differ = 1;
        }
    }
  else
    {
      if ( !( ip2->center    == ip1->center &&
              ip2->subcenter == ip1->subcenter )) differ = 1;
    }

  return differ;
}


static int instituteCompareP ( void *  instituteptr1, void * instituteptr2 )
{
  institute_t * i1, * i2;

  i1 = ( institute_t * ) instituteptr1;
  i2 = ( institute_t * ) instituteptr2;

  xassert(i1);
  xassert(i2);

  return instituteCompareKernel ( i1, i2 );
}

struct instLoc
{
  institute_t *ip;
  int id;
};

static enum cdiApplyRet
findInstitute(int id, void *res, void *data)
{
  institute_t * ip1 = ((struct instLoc *)data)->ip;
  institute_t * ip2 = res;
  if (ip2->used && !instituteCompareKernel(ip1, ip2))
    {
      ((struct instLoc *)data)->id = id;
      return CDI_APPLY_STOP;
    }
  else
    return CDI_APPLY_GO_ON;
}


int institutInq(int center, int subcenter, const char *name, const char *longname)
{
  instituteInit ();

  institute_t * ip_ref = xmalloc(sizeof (*ip_ref));
  ip_ref->self       = UNDEFID;
  ip_ref->used       = 0;
  ip_ref->center     = center;
  ip_ref->subcenter  = subcenter;
  ip_ref->name       = name && name[0] ? (char *)name : NULL;
  ip_ref->longname   = longname && longname[0] ? (char *)longname : NULL;

  struct instLoc state = { .ip = ip_ref, .id = UNDEFID };
  cdiResHFilterApply(&instituteOps, findInstitute, &state);

  free(ip_ref);

  return state.id;
}


int institutDef(int center, int subcenter, const char *name, const char *longname)
{
  institute_t * instituteptr;

  instituteInit ();

  instituteptr = instituteNewEntry();

  instituteptr->center    = center;
  instituteptr->subcenter = subcenter;
  if ( name && *name )         instituteptr->name     = strdupx(name);
  if ( longname && *longname ) instituteptr->longname = strdupx(longname);

  return instituteptr->self;
}


int institutInqCenter(int instID)
{
  institute_t * instituteptr = NULL;

  instituteInit ();

  if ( instID != UNDEFID )
    instituteptr = ( institute_t * ) reshGetVal ( instID, &instituteOps );

  return  instituteptr ? instituteptr->center : UNDEFID;
}


int institutInqSubcenter(int instID)
{
  institute_t * instituteptr = NULL;

  instituteInit ();

  if ( instID != UNDEFID )
    instituteptr = ( institute_t * ) reshGetVal ( instID, &instituteOps );

  return instituteptr ? instituteptr->subcenter: UNDEFID;
}


char *institutInqNamePtr(int instID)
{
  institute_t * instituteptr = NULL;

  instituteInit ();

  if ( instID != UNDEFID )
    instituteptr = ( institute_t * ) reshGetVal ( instID, &instituteOps );

  return instituteptr ? instituteptr->name : NULL;
}


char *institutInqLongnamePtr(int instID)
{
  institute_t * instituteptr = NULL;

  instituteInit ();

  if ( instID != UNDEFID )
    instituteptr = ( institute_t * ) reshGetVal ( instID, &instituteOps );

  return instituteptr ? instituteptr->longname : NULL;
}

static enum cdiApplyRet
activeInstitutes(int id, void *res, void *data)
{
  if (res && ((institute_t *)res)->used)
    ++(*(int *)data);
  return CDI_APPLY_GO_ON;
}

int institutInqNumber(void)
{
  int instNum = 0;

  instituteInit ();

  cdiResHFilterApply(&instituteOps, activeInstitutes, &instNum);
  return instNum;
}


void instituteDestroyP ( void * instituteptr )
{
  int id;
  institute_t * i1 = ( institute_t * ) instituteptr;

  xassert ( i1 );

  id = i1->self;

  if ( instituteptr ) free ( instituteptr );

  reshRemove ( id, &instituteOps );
}


void institutePrintP   ( void * instituteptr, FILE * fp )
{
  institute_t * ip = ( institute_t * ) instituteptr;

  if ( !ip ) return;

  fprintf ( fp, "#\n");
  fprintf ( fp, "# instituteID %d\n", ip->self);
  fprintf ( fp, "#\n");
  fprintf ( fp, "self          = %d\n", ip->self );
  fprintf ( fp, "used          = %d\n", ip->used );
  fprintf ( fp, "center        = %d\n", ip->center );
  fprintf ( fp, "subcenter     = %d\n", ip->subcenter );
  fprintf ( fp, "name          = %s\n", ip->name ? ip->name : "NN" );
  fprintf ( fp, "longname      = %s\n", ip->longname ? ip->longname : "NN" );
}


static int
instituteTxCode ( void )
{
  return INSTITUTE;
}

enum {
  institute_nints = 5,
};

static int instituteGetSizeP ( void * instituteptr, void *context)
{
  institute_t *p = instituteptr;
  int txsize = serializeGetSize(institute_nints, DATATYPE_INT, context)
    + serializeGetSize(strlen(p->name) + 1, DATATYPE_TXT, context)
    + serializeGetSize(strlen(p->longname) + 1, DATATYPE_TXT, context);
  return txsize;
}

static void institutePackP(void * instituteptr, void *buf, int size, int *position, void *context)
{
  institute_t *p = instituteptr;
  int tempbuf[institute_nints];
  tempbuf[0] = p->self;
  tempbuf[1] = p->center;
  tempbuf[2] = p->subcenter;
  tempbuf[3] = (int)strlen(p->name) + 1;
  tempbuf[4] = (int)strlen(p->longname) + 1;
  serializePack(tempbuf, institute_nints, DATATYPE_INT, buf, size, position, context);
  serializePack(p->name, tempbuf[3], DATATYPE_TXT, buf, size, position, context);
  serializePack(p->longname, tempbuf[4], DATATYPE_TXT, buf, size, position, context);
}

int instituteUnpack(void *buf, int size, int *position, int nspTarget, void *context)
{
  int tempbuf[institute_nints];
  int instituteID;
  char *name, *longname;
  serializeUnpack(buf, size, position, tempbuf, institute_nints, DATATYPE_INT, context);
  name = xmalloc(tempbuf[3]);
  longname = xmalloc(tempbuf[4]);
  serializeUnpack(buf, size, position, name, tempbuf[3], DATATYPE_TXT, context);
  serializeUnpack(buf, size, position, longname, tempbuf[4], DATATYPE_TXT, context);
  instituteID = institutDef(tempbuf[1], tempbuf[2], name, longname);
  // FIXME: this should work, once all types are transferred
  //xassert(instituteID == tempbuf[0]);
  return instituteID;
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
