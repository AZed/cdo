#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "dmemory.h"
#include "cdi.h"
#include "stream_int.h"

#undef  UNDEFID
#define UNDEFID  -1

int ECMWF  = UNDEFID;
int MPIMET = UNDEFID;
int DWD    = UNDEFID;
int MCH    = UNDEFID;

typedef struct
{
  int    used;  
  int    center;  
  int    subcenter;
  char  *name;
  char  *longname;
}
Institut;

Institut *instituts;

static int institutsSize = 0;
static int institutsNum  = 0;
static int InstitutsInit = 0;


void institutsDefault(void);
void institutsInit(void);


void institutsInitEntry(int instID)
{
  if ( instID < 0 || instID >= institutsSize )
    Error("instID %d undefined!", instID);

  instituts[instID].used       = 0;
  instituts[instID].center     = UNDEFID;
  instituts[instID].subcenter  = UNDEFID;
  instituts[instID].name       = NULL;
  instituts[instID].longname   = NULL;
}


int institutsNewEntry(void)
{
  int instID = 0;

  /*
    Look for a free slot in instituts.
    (Create the table the first time through).
  */
  if ( !institutsSize )
    {
      int i;

      institutsSize = 32;
      instituts = (Institut *) malloc(institutsSize*sizeof(Institut));
      if( instituts == NULL )
	{
          Message("institutsSize = %d", institutsSize);
	  SysError("Allocation of Institut failed");
	}

      for( i = 0; i < institutsSize; i++ )
	institutsInitEntry(i);
    }
  else
    {
      while( instID < institutsSize )
	{
	  if ( instituts[instID].used == 0 ) break;
	  instID++;
	}
    }
  /*
    If the table overflows, double its size.
  */
  if ( instID == institutsSize )
    {
      int i;

      institutsSize = 2*institutsSize;
      instituts = (Institut *) realloc(instituts, institutsSize*sizeof(Institut));
      if( instituts == NULL )
	{
          Message("institutsSize = %d", institutsSize);
	  SysError("Reallocation of Institut failed");
	}

      for( i = instID; i < institutsSize; i++ )
	institutsInitEntry(i);
    }

  instituts[instID].used = 1;
  institutsNum++;

  return (instID);
}


void institutsInit(void)
{
  InstitutsInit = 1;

  institutsDefault();
}


int institutInq(int center, int subcenter, const char *name, const char *longname)
{
  int instID;
  size_t len;
  int found;

  if ( ! InstitutsInit ) institutsInit();

  for( instID = 0; instID < institutsSize; instID++ )
    {
      if ( instituts[instID].used )
	{
	  if ( name )
	    {
	      found = 1;
	      if ( center    > 0 && instituts[instID].center    != center )    found = 0;
	      if ( subcenter > 0 && instituts[instID].subcenter != subcenter ) found = 0;

	      if ( found )
		{
		  if ( instituts[instID].name )
		    {
		      len = strlen(instituts[instID].name);
		      if ( memcmp(instituts[instID].name, name, len) == 0 ) break;
		    }
		}
	    }
	  else if ( longname )
	    {
	      if ( instituts[instID].longname )
		{
		  len = strlen(instituts[instID].longname);
		  if ( memcmp(instituts[instID].longname, longname, len) == 0 ) break;
		}
	    }
	  else
	    {
	      if ( instituts[instID].center    == center &&
		   instituts[instID].subcenter == subcenter ) break;
	    }
	}
    }

  if ( instID == institutsSize ) instID = UNDEFID;

  return (instID);
}


int institutDef(int center, int subcenter, const char *name, const char *longname)
{
  int instID;

  if ( ! InstitutsInit ) institutsInit();
  /*
  instID = institutInq(center, subcenter, name, longname);

  if ( instID == UNDEFID )
  */
    {
      instID = institutsNewEntry();

      instituts[instID].center    = center;
      instituts[instID].subcenter = subcenter;

      if ( name )     instituts[instID].name     = strdupx(name);
      if ( longname ) instituts[instID].longname = strdupx(longname);
    }

  return (instID);
}


void institutionCheckID(const char *caller, int instID)
{
  if ( instID < 0 || instID >= institutsSize )
    Errorc("instID %d undefined!", instID);

  if ( ! instituts[instID].used )
    Errorc("instID %d undefined!", instID);
}


int institutInqCenter(int instID)
{
  int center = UNDEFID;

  if ( ! InstitutsInit ) institutsInit();

  if ( instID != UNDEFID )
    {
      institutionCheckID(__func__, instID);

      center = instituts[instID].center;
    }

  return (center);
}


int institutInqSubcenter(int instID)
{
  int subcenter = UNDEFID;

  if ( ! InstitutsInit ) institutsInit();

  if ( instID != UNDEFID )
    {
      institutionCheckID(__func__, instID);

      subcenter = instituts[instID].subcenter;
    }

  return (subcenter);
}


char *institutInqNamePtr(int instID)
{
  char *name = NULL;

  if ( ! InstitutsInit ) institutsInit();

  if ( instID != UNDEFID )
    {
      institutionCheckID(__func__, instID);

      if ( instituts[instID].name )
	name = instituts[instID].name;
    }

  return (name);
}


char *institutInqLongnamePtr(int instID)
{
  char *name = NULL;

  if ( ! InstitutsInit ) institutsInit();

  if ( instID != UNDEFID )
    {
      institutionCheckID(__func__, instID);

      if ( instituts[instID].longname )
	name = instituts[instID].longname;
    }

  return (name);
}


int institutInqNumber(void)
{
  if ( ! InstitutsInit ) institutsInit();

  return (institutsNum);
}


void institutsDefault(void)
{
  ECMWF   = institutDef( 98,   0, "ECMWF",     "European Centre for Medium-Range Weather Forecasts");
  MPIMET  = institutDef( 98, 232, "MPIMET",    "Max-Planck-Institute for Meteorology");
     (void) institutDef( 98, 255, "MPIMET",    "Max-Planck-Institute for Meteorology");
     (void) institutDef( 98, 232, "MPIMET",    "Max-Planck Institute for Meteorology");
  DWD     = institutDef( 78, 255, "DWD",       "Deutscher Wetterdienst");
  MCH     = institutDef(215, 255, "MCH",       "MeteoSwiss");
     (void) institutDef(  7,   0, "NCEP",      "National Centers for Environmental Prediction");
     (void) institutDef(  7,   1, "NCEP",      "National Centers for Environmental Prediction");
     (void) institutDef( 60,   0, "NCAR",      "National Center for Atmospheric Research");
     (void) institutDef( 74,   0, "METOFFICE", "U.K. Met Office");
     (void) institutDef( 97,   0, "ESA",       "European Space Agency ");
     (void) institutDef( 99,   0, "KNMI",      "Royal Netherlands Meteorological Institute");
     /*     (void) institutDef(  0,   0, "IPSL",      "IPSL (Institut Pierre Simon Laplace, Paris, France)"); */
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
