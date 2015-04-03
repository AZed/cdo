#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <ctype.h>

#include "dmemory.h"
#include "cdi.h"
#include "stream_int.h"

#undef  UNDEFID
#define UNDEFID -1

/*int TableDefine = 0; */ /* Define new table also if the entry already exist */
                          /* This is needed for createtable */

#include "tablepar.h"
#include "table.h"

#define MAX_TABLE  256
#define MAX_PARS   1024

typedef struct
{
  int    used;  
  PAR   *pars;
  int    npars;
  int    modelID;
  int    number;
  char  *name;
} 
PARTAB;

static PARTAB parTable[MAX_TABLE];
static int  parTableSize = MAX_TABLE;
static int  parTableNum  = 0;
static int  ParTableInit = 0;

static char *tablePath = NULL;

void tableDefModelID(int tableID, int modelID);
void tableDefNum(int tableID, int tablenum);


void tableDefEntry(int tableID, int id, const char *name,
		   const char *longname, const char *units)
{
  int item;

  item = parTable[tableID].npars++;
  parTable[tableID].pars[item].id       = id;
  parTable[tableID].pars[item].name     = NULL;
  parTable[tableID].pars[item].longname = NULL;
  parTable[tableID].pars[item].units    = NULL;

  if ( name )
    if ( strlen(name) > 0 )
      parTable[tableID].pars[item].name     = strdupx(name);
  if ( longname )
    if ( strlen(longname) > 0 )
      parTable[tableID].pars[item].longname = strdupx(longname);
  if ( units )
    if ( strlen(units) > 0 )
      parTable[tableID].pars[item].units    = strdupx(units);
}

void tableLink(int tableID, PAR *pars, int npars)
{
  int item;

  for ( item = 0; item < npars; item++ )
    {
      parTable[tableID].pars[item].id       = pars[item].id;
      parTable[tableID].pars[item].name     = pars[item].name;
      parTable[tableID].pars[item].longname = pars[item].longname;
      parTable[tableID].pars[item].units    = pars[item].units;
    }

  parTable[tableID].npars = npars;
}

void parTableInitEntry(int tableID)
{
  parTable[tableID].used    = 0;
  parTable[tableID].pars    = NULL;
  parTable[tableID].npars   = 0;
  parTable[tableID].modelID = UNDEFID;
  parTable[tableID].number  = UNDEFID;
  parTable[tableID].name    = NULL;
}

void tableGetPath(void)
{
  char *path;

  path = getenv("TABLEPATH");

  if ( path ) tablePath = strdupx(path);
  /*
  printf("tablePath = %s\n", tablePath);
  */
}

void parTableInit(void)
{
  ParTableInit = 1;

  if ( cdiPartabIntern )
    tableDefault();

  tableGetPath();
}

int tableNewEntry()
{
  int tableID = 0;
  static int init = 0;

  if ( ! init )
    {
      for ( tableID = 0; tableID < parTableSize; tableID++ )
	parTableInitEntry(tableID);
      init = 1;
    }

  /*
    Look for a free slot in parTable.
  */
  for ( tableID = 0; tableID < parTableSize; tableID++ )
    {
      if ( ! parTable[tableID].used ) break;
    }

  if ( tableID == parTableSize )
    Error("no more entries!");

  parTable[tableID].used = 1;
  parTableNum++;

  return (tableID);
}

int decodeForm1(char *pline, char *name, char *longname, char *units)
{
  /* Format 1 : code name add mult longname [units] */
  double add, mult;
  int level;
  char *pstart, *pend;
  long len;

  level = strtol(pline, &pline, 10);
  while ( isspace((int) *pline) ) pline++;

  pstart = pline;
  while ( ! (isspace((int) *pline) || *pline == 0) ) pline++;
  len = pline - pstart;
  if ( len > 0 )
    {
      memcpy(name, pstart, len);
      name[len] = 0;
    }
  else
    return (0);

  len = strlen(pline);
  if ( len == 0 ) return (0);

  add  = strtod(pline, &pline);
  mult = strtod(pline, &pline);

  while ( isspace((int) *pline) ) pline++;

  len = strlen(pline);
  if ( len > 0)
    {
      pstart = pline;
      pend = strrchr(pline, '[');
      if ( pend )
	pend--;
      else
	pend = pstart + len;
      while ( isspace((int) *pend) ) pend--;
      len = pend - pstart + 1;
      if ( len > 0 )
	{
	  memcpy(longname, pstart, len);
	  longname[len] = 0;
	}
      pstart = strrchr(pline, '[');
      if ( pstart )
	{
	  pstart++;
	  while ( isspace((int) *pstart) ) pstart++;
	  pend = strchr(pstart, ']');
	  if ( ! pend ) return (0);
	  pend--;
	  while ( isspace((int) *pend) ) pend--;
	  len = pend - pstart + 1;
	  if ( len > 0 )
	    {
	      memcpy(units, pstart, len);
	      units[len] = 0;
	    }	  
	}
    }
 
  return (0);
}

int decodeForm2(char *pline, char *name, char *longname, char *units)
{
  /* Format 2 : code | name | longname | units */
  char *pend;
  long len;

  pline = strchr(pline, '|');
  pline++;

  while ( isspace((int) *pline) ) pline++;
  pend = strchr(pline, '|');
  if ( ! pend )
    {
      pend = pline;
      while ( ! isspace((int) *pend) ) pend++;
      len = pend - pline;
      if ( len > 0 )
	{
	  memcpy(name, pline, len);
	  name[len] = 0;
	}
      return (0);
    }
  else
    {
      pend--;
      while ( isspace((int) *pend) ) pend--;
      len = pend - pline + 1;
      if ( len > 0 )
	{
	  memcpy(name, pline, len);
	  name[len] = 0;
	}
    }

  pline = strchr(pline, '|');
  pline++;
  while ( isspace((int) *pline) ) pline++;
  pend = strchr(pline, '|');
  if ( !pend ) pend = strchr(pline, 0);
  pend--;
  while ( isspace((int) *pend) ) pend--;
  len = pend - pline + 1;
  if ( len > 0 )
    {
      memcpy(longname, pline, len);
      longname[len] = 0;
    }

  pline = strchr(pline, '|');
  if ( pline )
    {
      pline++;
      while ( isspace((int) *pline) ) pline++;
      pend = strchr(pline, '|');
      if ( !pend ) pend = strchr(pline, 0);
      pend--;
      while ( isspace((int) *pend) ) pend--;
      len = pend - pline + 1;
      if ( len < 0 ) len = 0;
      memcpy(units, pline, len);
      units[len] = 0;
    }

  return (0);
}

int tableRead(const char *tablefile)
{
  char line[1024], *pline;
  int lnr = 0;
  long len;
  int id;
  char name[256], longname[256], units[256];
  int tableID = UNDEFID;
  int err;
  char *tablename;
  FILE *tablefp;

  tablefp = fopen(tablefile, "r");
  if ( tablefp == NULL ) return (tableID);

  tablename = strrchr(tablefile, '/');
  if ( tablename == 0 ) tablename = (char *) tablefile;
  else                  tablename++;

  tableID = tableDef(-1, 0, tablename);

  while ( fgets(line, 1023, tablefp) )
    {
      len = strlen(line);
      if ( line[len-1] == '\n' ) line[len-1] = '\0';
      lnr++;
      id       = CDI_UNDEFID;
      name[0]     = 0;
      longname[0] = 0;
      units[0]    = 0;
      if ( line[0] == '#' ) continue;
      pline = line;

      len = strlen(pline);
      if ( len < 4 ) continue;
      while ( isspace((int) *pline) ) pline++;
      id = atoi(pline);
      /*
      if ( id > 255 ) id -= 256;
      */
      if ( id == 0 ) continue;

      while ( isdigit((int) *pline) ) pline++; 

      if ( strchr(pline, '|') )
	err = decodeForm2(pline, name, longname, units);
      else
	err = decodeForm1(pline, name, longname, units);

      if ( err ) continue;

      if ( strlen(name) == 0 ) sprintf(name, "var%d", id);

      tableDefEntry(tableID, id, name, longname, units);
    }

  return (tableID);
}

int tableFromEnv(int modelID, int tablenum)
{
  int tableID = UNDEFID;
  char tablename[256] = {'\0'};
  int tablenamefound = 0;

  if ( modelInqNamePtr(modelID) )
    {
      strcpy(tablename, modelInqNamePtr(modelID));
      if ( tablenum )
	{
	  int len = strlen(tablename);
	  sprintf(tablename+len, "_%03d", tablenum);
	}
      tablenamefound = 1;
    }
  else
    {
      int instID = modelInqInstitut(modelID);
      if ( instID != UNDEFID )
	{
	  if ( institutInqNamePtr(instID) )
	    {
	      strcpy(tablename, institutInqNamePtr(instID));
	      if ( tablenum )
		{
		  int len = strlen(tablename);
		  sprintf(tablename+len, "_%03d", tablenum);
		}
	      tablenamefound = 1;
	    }
	}
    }

  if ( tablenamefound )
    {
      int lenp = 0, lenf;
      char *tablefile = NULL;
      if ( tablePath )
	lenp = strlen(tablePath);
      lenf = strlen(tablename);
      /* if (tablePath) printf("tablePath = %s\n", tablePath); */
      /* if (tablename) printf("tableName = %s\n", tablename); */
      tablefile = (char *) malloc(lenp+lenf+3);
      if ( tablePath )
	{
	  strcpy(tablefile, tablePath);
	  strcat(tablefile, "/");
	}
      else
	tablefile[0] = '\0';
      strcat(tablefile, tablename);
      /* if (tablefile) printf("tableFile = %s\n", tablefile); */

      tableID = tableRead(tablefile);
      if ( tableID != UNDEFID )
	{
	  tableDefModelID(tableID, modelID);
	  tableDefNum(tableID, tablenum);
	}
      /* printf("tableID = %d %s\n", tableID, tablefile); */

      free(tablefile);
    }

  return (tableID);
}

int tableInq(int modelID, int tablenum, const char *tablename)
{
  int tableID = UNDEFID;
  int modelID2 = UNDEFID, i, len;
  char tablefile[256] = {'\0'};

  if ( ! ParTableInit ) parTableInit();

  if ( tablename )
    {
      size_t len;
      strcpy(tablefile, tablename);
      /*
      printf("tableInq: tablefile = >%s<\n", tablefile);
      */
      /* search for internal table */
      for ( tableID = 0; tableID < MAX_TABLE; tableID++ )
	{
	  if ( parTable[tableID].used && parTable[tableID].name )
	    {
	      /* len = strlen(parTable[tableID].name); */
	      len = strlen(tablename);
	      if ( memcmp(parTable[tableID].name, tablename, len) == 0 ) break;
	    }
	}
      if ( tableID == MAX_TABLE ) tableID = UNDEFID;
      if ( CDI_Debug )
	Message("tableID = %d tablename = %s", tableID, tablename);
    }
  else
    {
      for ( tableID = 0; tableID < MAX_TABLE; tableID++ )
	{
	  if ( parTable[tableID].used )
	    {	  
	      if ( parTable[tableID].modelID == modelID &&
		   parTable[tableID].number  == tablenum ) break;
	    }
	}
  
      if ( tableID == MAX_TABLE ) tableID = UNDEFID;

      if ( tableID == UNDEFID )
	{
	  if ( modelID != UNDEFID )
	    {
	      if ( modelInqNamePtr(modelID) )
		{
		  strcpy(tablefile, modelInqNamePtr(modelID));
		  len = strlen(tablefile);
		  for ( i = 0; i < len; i++)
		    if ( tablefile[i] == '.' ) tablefile[i] = '\0';
		  modelID2 = modelInq(-1, 0, tablefile);
		}
	    }
	  if ( modelID2 != UNDEFID )
	    for ( tableID = 0; tableID < MAX_TABLE; tableID++ )
	      {
		if ( parTable[tableID].used )
		  {
		    if ( parTable[tableID].modelID == modelID2 &&
			 parTable[tableID].number  == tablenum ) break;
		  }
	      }
	}

      if ( tableID == MAX_TABLE ) tableID = UNDEFID;

      if ( tableID == UNDEFID && modelID != UNDEFID )
	tableID = tableFromEnv(modelID, tablenum);

      if ( CDI_Debug )
	if ( tablename )
	  Message("tableID = %d tablename = %s", tableID, tablename);
    }

  return (tableID);
}

int tableDef(int modelID, int tablenum, const char *tablename)
{
  int tableID = UNDEFID;

  if ( ! ParTableInit ) parTableInit();
  /*
  if ( ! (modelID == UNDEFID && tablenum == 0) )
    tableID = tableInq(modelID, tablenum, tablename);
    */
  if ( tableID == UNDEFID )
    {
      tableID = tableNewEntry();

      parTable[tableID].modelID = modelID;
      parTable[tableID].number  = tablenum;
      if ( tablename ) 
	parTable[tableID].name = strdupx(tablename);

      parTable[tableID].pars = (PAR *) malloc(MAX_PARS * sizeof(PAR));
    }

  return (tableID);
}

void tableDefModelID(int tableID, int modelID)
{
  parTable[tableID].modelID = modelID;
}

void tableDefNum(int tableID, int tablenum)
{
  parTable[tableID].number  = tablenum;
}

int tableInqNum(int tableID)
{
  int number = 0;

  if ( tableID >= 0 && tableID < MAX_TABLE )
    number = parTable[tableID].number;

  return (number);
}

int tableInqModel(int tableID)
{
  int modelID = -1;

  if ( tableID >= 0 && tableID < MAX_TABLE )
    modelID = parTable[tableID].modelID;

  return (modelID);
}

void partabCheckID(int item)
{
  if ( item < 0 || item >= parTableSize )
    Error("item %d undefined!", item);

  if ( ! parTable[item].name )
    Error("item %d name undefined!", item);
}

char *tableInqNamePtr(int tableID)
{
  char *tablename = NULL;

  if ( CDI_Debug )
    Message("tableID = %d", tableID);

  if ( ! ParTableInit ) parTableInit();

  if ( tableID >= 0 && tableID < parTableSize )
    if ( parTable[tableID].name )
      tablename = parTable[tableID].name;

  return (tablename);
}

void tableWrite(const char *ptfile, int tableID)
{
  int item, npars;
  int lenname, lenlname, lenunits;
  int maxname = 4, maxlname = 10, maxunits = 2;
  FILE *ptfp;
  int tablenum, modelID, instID = CDI_UNDEFID;
  int center = 0, subcenter = 0;
  char *name, *longname, *units;
  char *instnameptr = NULL, *modelnameptr = NULL;

  if ( CDI_Debug )
    Message("write parameter table %d to %s", tableID, ptfile);

  if ( tableID == UNDEFID )
    {
      Warning("parameter table ID undefined");
      return;
    }

  partabCheckID(tableID);

  ptfp = fopen(ptfile, "w");

  npars = parTable[tableID].npars;

  for ( item = 0; item < npars; item++)
    {
      if ( parTable[tableID].pars[item].name )
	{
	  lenname  = strlen(parTable[tableID].pars[item].name);
	  if ( lenname  > maxname )  maxname  = lenname;
	}

      if ( parTable[tableID].pars[item].longname )
	{
	  lenlname = strlen(parTable[tableID].pars[item].longname);
	  if ( lenlname > maxlname ) maxlname = lenlname;
	}

      if ( parTable[tableID].pars[item].units )
	{
	  lenunits = strlen(parTable[tableID].pars[item].units);
	  if ( lenunits > maxunits ) maxunits = lenunits;
	}
    }

  tablenum = tableInqNum(tableID);
  modelID = parTable[tableID].modelID;
  if ( modelID != CDI_UNDEFID )
    {
      modelnameptr = modelInqNamePtr(modelID);
      instID = modelInqInstitut(modelID);
    }
  if ( instID != CDI_UNDEFID )
    {
      center = institutInqCenter(instID);
      subcenter = institutInqSubcenter(instID);
      instnameptr = institutInqNamePtr(instID);
    }

  fprintf(ptfp, "# Parameter table\n");
  fprintf(ptfp, "#\n");
  if ( tablenum )
    fprintf(ptfp, "# TABLE_ID=%d\n", tablenum);
  fprintf(ptfp, "# TABLE_NAME=%s\n", parTable[tableID].name);
  if ( modelnameptr )
    fprintf(ptfp, "# TABLE_MODEL=%s\n", modelnameptr);
  if ( instnameptr )
    fprintf(ptfp, "# TABLE_INSTITUT=%s\n", instnameptr);
  if ( center )
    fprintf(ptfp, "# TABLE_CENTER=%d\n", center);
  if ( subcenter )
    fprintf(ptfp, "# TABLE_SUBCENTER=%d\n", subcenter);
  fprintf(ptfp, "#\n");
  fprintf(ptfp, "#\n");
  fprintf(ptfp, "# id       = parameter ID\n");
  fprintf(ptfp, "# name     = variable name\n");
  fprintf(ptfp, "# title    = long name (description)\n");
  fprintf(ptfp, "# units    = variable units\n");
  fprintf(ptfp, "#\n");
  fprintf(ptfp, "# The format of each record is:\n");
  fprintf(ptfp, "#\n");
  fprintf(ptfp, "# id | %-*s | %-*s | %-*s\n",
	  maxname,  "name",
	  maxlname, "title",
	  maxunits, "units");
	  
  for ( item = 0; item < npars; item++)
    {
      name = parTable[tableID].pars[item].name;
      longname = parTable[tableID].pars[item].longname;
      units = parTable[tableID].pars[item].units;
      if ( name == NULL ) name = " ";
      if ( longname == NULL ) longname = " ";
      if ( units == NULL ) units = " ";
      fprintf(ptfp, "%4d | %-*s | %-*s | %-*s\n",
	      parTable[tableID].pars[item].id,
	      maxname, name,
	      maxlname, longname,
	      maxunits, units);
    }

  fclose(ptfp);
}


void tableWriteC(const char *filename, int tableID)
{
  char chelp[] = "";
  int item, npars;
  int lenname, lenlname, lenunits;
  int maxname = 0, maxlname = 0, maxunits = 0;
  char tablename[256];
  int len, i;
  FILE *ptfp;

  if ( CDI_Debug )
    Message("write parameter table %d to %s", tableID, filename);

  if ( tableID == UNDEFID )
    {
      Warning("parameter table ID undefined");
      return;
    }

  partabCheckID(tableID);

  ptfp = fopen(filename, "w");

  npars = parTable[tableID].npars;

  for ( item = 0; item < npars; item++)
    {
      if ( parTable[tableID].pars[item].name )
	{
	  lenname  = strlen(parTable[tableID].pars[item].name);
	  if ( lenname  > maxname )  maxname  = lenname;
	}

      if ( parTable[tableID].pars[item].longname )
	{
	  lenlname = strlen(parTable[tableID].pars[item].longname);
	  if ( lenlname > maxlname ) maxlname = lenlname;
	}

      if ( parTable[tableID].pars[item].units )
	{
	  lenunits = strlen(parTable[tableID].pars[item].units);
	  if ( lenunits > maxunits ) maxunits = lenunits;
	}
    }

  strcpy(tablename, parTable[tableID].name);
  len = strlen(tablename);

  for ( i = 0; i < len; i++ )
    if ( tablename[i] == '.' ) tablename[i] = '_';

  fprintf(ptfp, "static PAR %s[] = {\n", tablename);
	  
  for ( item = 0; item < npars; item++ )
    {
      len = strlen(parTable[tableID].pars[item].name);
      fprintf(ptfp, "  {%4d, \"%s\", %-*s",
	      parTable[tableID].pars[item].id,
	      parTable[tableID].pars[item].name, maxname-len, chelp);

      if ( parTable[tableID].pars[item].longname )
	len = strlen(parTable[tableID].pars[item].longname);
      else
	len = 0;

      if ( len == 0 )
	fprintf(ptfp, " NULL, %-*s", maxlname-3, chelp);
      else
	fprintf(ptfp, "\"%s\", %-*s",
		parTable[tableID].pars[item].longname, maxlname-len, chelp);

      if ( parTable[tableID].pars[item].units )
	len = strlen(parTable[tableID].pars[item].units);
      else
	len = 0;

      if ( len == 0 )
	fprintf(ptfp, " NULL %-*s},\n", maxunits-3, chelp);
      else
	fprintf(ptfp, "\"%s\" %-*s},\n",
		parTable[tableID].pars[item].units,
		maxunits-len, chelp);
    }

  fprintf(ptfp, "};\n\n");

  fclose(ptfp);
}


int tableInqParCode(int tableID, char *varname, int *code)
{
  int item, npars;
  int err = 0;

  npars = parTable[tableID].npars;

  if ( tableID == UNDEFID || varname == NULL )
    {
      err = 1;
    }
  else
    {
      for ( item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].name )
	    if ( strcmp(parTable[tableID].pars[item].name, varname) == 0 )
	      {
		*code = parTable[tableID].pars[item].id;
		break;
	      }
	}
      if ( item == npars ) err = 1;
    }

  return (err);
}


int tableInqParName(int tableID, int code, char *varname)
{
  int item, npars;
  int err = 0;

  npars = parTable[tableID].npars;

  if ( tableID == UNDEFID )
    {
      err = 1;
    }
  else
    {
      for ( item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].id == code )
	    {
	      if ( parTable[tableID].pars[item].name )
		strcpy(varname, parTable[tableID].pars[item].name);
	      break;
	    }
	}
      if ( item == npars ) err = 1;
    }

  return (err);
}


char *tableInqParNamePtr(int tableID, int code)
{
  char *name = NULL;
  int item, npars;

  if ( tableID != UNDEFID )
    {
      npars = parTable[tableID].npars;
      for ( item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].id == code )
	    {
	      name = parTable[tableID].pars[item].name;
	      break;
	    }
	}
    }

  return (name);
}


char *tableInqParLongnamePtr(int tableID, int code)
{
  char *longname = NULL;
  int item, npars;

  if ( tableID != UNDEFID )
    {
      npars = parTable[tableID].npars;
      for ( item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].id == code )
	    {
	      longname = parTable[tableID].pars[item].longname;
	      break;
	    }
	}
    }

  return (longname);
}


char *tableInqParUnitsPtr(int tableID, int code)
{
  char *units = NULL;
  int item, npars;

  if ( tableID != UNDEFID )
    {
      npars = parTable[tableID].npars;
      for ( item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].id == code )
	    {
	      units = parTable[tableID].pars[item].units;
	      break;
	    }
	}
    }

  return (units);
}


int tableInqParLongname(int tableID, int code, char *longname)
{
  int item, npars;
  int err = 0;

  npars = parTable[tableID].npars;

  if ( tableID == UNDEFID )
    {
      err = 1;
    }
  else
    {
      for ( item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].id == code )
	    {
	      if ( parTable[tableID].pars[item].longname )
		strcpy(longname, parTable[tableID].pars[item].longname);
	      break;
	    }
	}
      if ( item == npars ) err = 1;
    }

  return (err);
}


int tableInqParUnits(int tableID, int code, char *units)
{
  int item, npars;
  int err = 0;

  npars = parTable[tableID].npars;

  if ( tableID == UNDEFID )
    {
      err = 1;
    }
  else
    {
      for ( item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].id == code )
	    {
	      if ( parTable[tableID].pars[item].units )
		strcpy(units, parTable[tableID].pars[item].units);
	      break;
	    }
	}
      if ( item == npars ) err = 1;
    }

  return (err);
}


void tableInqPar(int tableID, int code, char *name, char *longname, char *units)
{
  int item, npars;

  npars = parTable[tableID].npars;

  for ( item = 0; item < npars; item++ )
    {
      if ( parTable[tableID].pars[item].id == code )
	{
	  if ( parTable[tableID].pars[item].name )
	    strcpy(name, parTable[tableID].pars[item].name);
	  if ( parTable[tableID].pars[item].longname )
	    strcpy(longname, parTable[tableID].pars[item].longname);
	  if ( parTable[tableID].pars[item].units )
	    strcpy(units, parTable[tableID].pars[item].units);
	  break;
	}
    }
}


int parInqID(int tableID, int code)
{
  int item, npars;

  npars = parTable[tableID].npars;

  for ( item = 0; item < npars; item++ )
    {
      if ( parTable[tableID].pars[item].id == code ) break;
    }

  if ( item == npars ) item = -1;

  return (item);
}

int tableInqNumber(void)
{
  if ( ! ParTableInit ) parTableInit();

  return (parTableNum);
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
