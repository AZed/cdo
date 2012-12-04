#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <limits.h>

#include "dmemory.h"
#include "cdi.h"
#include "stream_int.h"
#include "vlist.h"


static
void vlistvarInitEntry(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistptr->vars[varID].fvarID        = varID;
  vlistptr->vars[varID].mvarID        = varID;
  vlistptr->vars[varID].flag          = 0;
  vlistptr->vars[varID].param         = 0;
  vlistptr->vars[varID].timeID        = CDI_UNDEFID;
  vlistptr->vars[varID].datatype      = CDI_UNDEFID;
  vlistptr->vars[varID].tsteptype     = TSTEP_INSTANT;
  vlistptr->vars[varID].timave        = 0;
  vlistptr->vars[varID].timaccu       = 0;
  vlistptr->vars[varID].xyz           = 0;
  vlistptr->vars[varID].gridID        = CDI_UNDEFID;
  vlistptr->vars[varID].zaxisID       = CDI_UNDEFID;
  vlistptr->vars[varID].instID        = CDI_UNDEFID;
  vlistptr->vars[varID].modelID       = CDI_UNDEFID;
  vlistptr->vars[varID].tableID       = CDI_UNDEFID;
  vlistptr->vars[varID].missvalused   = FALSE;
  vlistptr->vars[varID].missval       = cdiDefaultMissval;
  vlistptr->vars[varID].addoffset     = 0.0;
  vlistptr->vars[varID].scalefactor   = 1.0;
  vlistptr->vars[varID].name          = NULL;
  vlistptr->vars[varID].longname      = NULL;
  vlistptr->vars[varID].stdname       = NULL;
  vlistptr->vars[varID].units         = NULL;
  vlistptr->vars[varID].nlevs         = 0;
  vlistptr->vars[varID].levinfo       = NULL;
  vlistptr->vars[varID].comptype      = COMPRESS_NONE;
  vlistptr->vars[varID].complevel     = 1;
  vlistptr->vars[varID].atts.nalloc   = MAX_ATTRIBUTES;
  vlistptr->vars[varID].atts.nelems   = 0;
  vlistptr->vars[varID].lvalidrange   = 0;
  vlistptr->vars[varID].validrange[0] = VALIDMISS;
  vlistptr->vars[varID].validrange[1] = VALIDMISS;
}

static
int vlistvarNewEntry(int vlistID)
{
  int varID = 0;
  int vlistvarSize;
  var_t *vlistvar;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistvarSize = vlistptr->varsAllocated;
  vlistvar     = vlistptr->vars;
  /*
    Look for a free slot in vlistvar.
    (Create the table the first time through).
  */
  if ( ! vlistvarSize )
    {
      int i;

      vlistvarSize = 2;
      vlistvar = (var_t *) malloc(vlistvarSize*sizeof(var_t));
      if ( vlistvar == NULL )
	{
          Message("vlistvarSize = %d", vlistvarSize);
	  SysError("Allocation of var_t failed");
	}

      for ( i = 0; i < vlistvarSize; i++ )
	vlistvar[i].isUsed = FALSE;
    }
  else
    {
      while ( varID < vlistvarSize )
	{
	  if ( ! vlistvar[varID].isUsed ) break;
	  varID++;
	}
    }
  /*
    If the table overflows, double its size.
  */
  if ( varID == vlistvarSize )
    {
      int i;

      vlistvarSize = 2*vlistvarSize;
      vlistvar = (var_t *) realloc(vlistvar, vlistvarSize*sizeof(var_t));
      if ( vlistvar == NULL )
	{
          Message("vlistvarSize = %d", vlistvarSize);
	  SysError("Reallocation of var_t failed");
	}
      varID = vlistvarSize/2;

      for ( i = varID; i < vlistvarSize; i++ )
	vlistvar[i].isUsed = FALSE;
    }

  vlistptr->varsAllocated = vlistvarSize;
  vlistptr->vars          = vlistvar;

  vlistvarInitEntry(vlistID, varID);

  vlistptr->vars[varID].isUsed = TRUE;

  return (varID);
}

static
void vlistCheckVarID(const char *caller, int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  if ( vlistptr == NULL )
    Errorc("vlist undefined!");

  if ( varID < 0 || varID >= vlistptr->nvars )
    Errorc("varID %d undefined!", varID);

  if ( ! vlistptr->vars[varID].isUsed )
    Errorc("varID %d undefined!", varID);
}

/*
@Function  vlistDefVar
@Title     Define a Variable

@Prototype int vlistDefVar(int vlistID, int gridID, int zaxisID, int timeID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate}.
    @Item  timeID   One of the set of predefined CDI time identifiers.
                    The valid CDI time identifiers are @func{TIME_CONSTANT} and @func{TIME_VARIABLE}.

@Description
The function @func{vlistDefVar} adds a new variable to vlistID.

@Result
@func{vlistDefVar} returns an identifier to the new variable.

@Example
Here is an example using @func{vlistCreate} to create a variable list
and add a variable with @func{vlistDefVar}.

@Source
#include "cdi.h"
   ...
int vlistID, varID;
   ...
vlistID = vlistCreate();
varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
   ...
streamDefVlist(streamID, vlistID);
   ...
vlistDestroy(vlistID);
   ...
@EndSource
@EndFunction
*/
int vlistDefVar(int vlistID, int gridID, int zaxisID, int timeID)
{
  int varID;
  int nlevs;
  int levID;
  int index;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  if ( CDI_Debug )
    Message("gridID = %d  zaxisID = %d  timeID = %d", gridID, zaxisID, timeID);

  varID = vlistvarNewEntry(vlistID);

  vlistptr->nvars++;

  vlistptr->vars[varID].gridID  = gridID;
  vlistptr->vars[varID].zaxisID = zaxisID;
  vlistptr->vars[varID].timeID  = timeID;

  if ( timeID != TIME_VARIABLE && timeID != TIME_CONSTANT )
    {
      Message("unexpected timeID %d. Set to TIME_VARIABLE", timeID);
      vlistptr->vars[varID].timeID = TIME_VARIABLE;	  
    }

  nlevs = zaxisInqSize(zaxisID);

  vlistptr->vars[varID].levinfo = (levinfo_t *) malloc(nlevs*sizeof(levinfo_t));

  for ( levID = 0; levID < nlevs; levID++ )
    {
      vlistptr->vars[varID].levinfo[levID].flag     = 0;
      vlistptr->vars[varID].levinfo[levID].index    = -1;
      vlistptr->vars[varID].levinfo[levID].flevelID = levID;
      vlistptr->vars[varID].levinfo[levID].mlevelID = levID;
    }

  vlistptr->vars[varID].nlevs = nlevs;

  for ( index = 0; index < vlistptr->ngrids; index++ )
    if ( gridID == vlistptr->gridIDs[index] ) break;

  if ( index == vlistptr->ngrids )
    {
      if ( vlistptr->ngrids + 1 >= MAX_GRIDS_PS )
	Error("Maximum of %d grids reached", MAX_GRIDS_PS);

      vlistptr->gridIDs[vlistptr->ngrids] = gridID;
      vlistptr->ngrids++;
    }

  for ( index = 0; index < vlistptr->nzaxis; index++ )
    if ( zaxisID == vlistptr->zaxisIDs[index] ) break;

  if ( index == vlistptr->nzaxis )
    {
      if ( vlistptr->nzaxis + 1 >= MAX_ZAXES_PS )
	Error("Maximum of %d zaxis reached", MAX_ZAXES_PS);

      vlistptr->zaxisIDs[vlistptr->nzaxis] = zaxisID;
      vlistptr->nzaxis++;
    }

  vlistptr->vars[varID].param = cdiEncodeParam(-(varID + 1), 255, 255);

  return (varID);
}

/*
@Function  vlistDefVarParam
@Title     Define the parameter number of a Variable

@Prototype void vlistDefVarParam(int vlistID, int varID, int param)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  param    Parameter number.

@Description
The function @func{vlistDefVarParam} defines the parameter number of a variable.

@EndFunction
*/
void vlistDefVarParam(int vlistID, int varID, int param)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  vlistptr->vars[varID].param = param;
}

/*
@Function  vlistDefVarCode
@Title     Define the code number of a Variable

@Prototype void vlistDefVarCode(int vlistID, int varID, int code)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  code     Code number.

@Description
The function @func{vlistDefVarCode} defines the code number of a variable.

@EndFunction
*/
void vlistDefVarCode(int vlistID, int varID, int code)
{
  vlist_t *vlistptr;
  int param, pnum, pcat, pdis;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  param = vlistptr->vars[varID].param;

  cdiDecodeParam(param, &pnum, &pcat, &pdis);
  
  vlistptr->vars[varID].param = cdiEncodeParam(code, pcat, pdis);
}


void vlistInqVar(int vlistID, int varID, int *gridID, int *zaxisID, int *timeID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  *gridID  = vlistptr->vars[varID].gridID;
  *zaxisID = vlistptr->vars[varID].zaxisID;
  *timeID  = vlistptr->vars[varID].timeID;

  return;
}

/*
@Function  vlistInqVarGrid
@Title     Get the Grid ID of a Variable

@Prototype int vlistInqVarGrid(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarGrid} returns the grid ID of a variable.

@Result
@func{vlistInqVarGrid} returns the grid ID of the variable.

@EndFunction
*/
int vlistInqVarGrid(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].gridID);
}

/*
@Function  vlistInqVarZaxis
@Title     Get the Zaxis ID of a Variable

@Prototype int vlistInqVarZaxis(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarZaxis} returns the zaxis ID of a variable.

@Result
@func{vlistInqVarZaxis} returns the zaxis ID of the variable.

@EndFunction
*/
int vlistInqVarZaxis(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].zaxisID);
}


int vlistInqVarTime(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].timeID);
}

/*
@Function  vlistInqVarParam
@Title     Get the parameter number of a Variable

@Prototype int vlistInqVarParam(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarParam} returns the parameter number of a variable.

@Result
@func{vlistInqVarParam} returns the parameter number of the variable.

@EndFunction
*/
int vlistInqVarParam(int vlistID, int varID)
{
  vlist_t *vlistptr;
  int param;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  param = vlistptr->vars[varID].param;

  return (param);
}

/*
@Function  vlistInqVarCode
@Title     Get the Code number of a Variable

@Prototype int vlistInqVarCode(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarCode} returns the code number of a variable.

@Result
@func{vlistInqVarCode} returns the code number of the variable.

@EndFunction
*/
int vlistInqVarCode(int vlistID, int varID)
{
  vlist_t *vlistptr;
  int param, code;
  int pdis, pcat, pnum;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  param = vlistptr->vars[varID].param;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);
  code = pnum;

  if ( code < 0 && vlistptr->vars[varID].tableID != -1 && vlistptr->vars[varID].name != NULL )
    {
      tableInqParCode(vlistptr->vars[varID].tableID, vlistptr->vars[varID].name, &code);
    }

  return (code);
}


const char *vlistInqVarNamePtr(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].name);
}


const char *vlistInqVarLongnamePtr(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].longname);
}


const char *vlistInqVarStdnamePtr(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].stdname);
}


const char *vlistInqVarUnitsPtr(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].units);
}

/*
@Function  vlistInqVarName
@Title     Get the name of a Variable

@Prototype void vlistInqVarName(int vlistID, int varID, char *name)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  name     Returned variable name. The caller must allocate space for the 
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarName} returns the name of a variable.

@Result
@func{vlistInqVarName} returns the name of the variable to the parameter name if available,
otherwise the result is an empty string.

@EndFunction
*/
void vlistInqVarName(int vlistID, int varID, char *name)
{
  int tableID;
  int param;
  int pdis, pcat, pnum;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( vlistptr->vars[varID].name == NULL )
    {
      param = vlistptr->vars[varID].param;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      if ( pdis == 255 )
	{
	  int code = pnum;
	  tableID = vlistptr->vars[varID].tableID;
	  if ( tableInqParName(tableID, code, name) != 0 )
	    sprintf(name, "var%d", code);
	}
      else
	{
	  sprintf(name, "param%d.%d.%d", pnum, pcat, pdis);
	}
    }  
  else
    strcpy(name, vlistptr->vars[varID].name);

  return;
}

/*
@Function  vlistInqVarLongname
@Title     Get the longname of a Variable

@Prototype void vlistInqVarLongname(int vlistID, int varID, char *longname)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  longname Long name of the variable. The caller must allocate space for the 
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarLongname} returns the longname of a variable if available,
otherwise the result is an empty string.

@Result
@func{vlistInqVaeLongname} returns the longname of the variable to the parameter longname.

@EndFunction
*/
void vlistInqVarLongname(int vlistID, int varID, char *longname)
{
  int tableID;
  int param;
  int pdis, pcat, pnum;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  longname[0] = '\0';

  if ( vlistptr->vars[varID].longname == NULL )
    {
      param = vlistptr->vars[varID].param;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      if ( pdis == 255 )
	{
	  int code = pnum;
	  tableID = vlistptr->vars[varID].tableID;
	  if ( tableInqParLongname(tableID, code, longname) != 0 )
	    longname[0] = '\0';
	}
    }  
  else
    strcpy(longname, vlistptr->vars[varID].longname);

  return;
}

/*
@Function  vlistInqVarStdname
@Title     Get the standard name of a Variable

@Prototype void vlistInqVarStdname(int vlistID, int varID, char *stdname)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  stdname  Standard name of the variable. The caller must allocate space for the 
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarStdname} returns the standard name of a variable if available,
otherwise the result is an empty string.

@Result
@func{vlistInqVarName} returns the standard name of the variable to the parameter stdname.

@EndFunction
*/
void vlistInqVarStdname(int vlistID, int varID, char *stdname)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( vlistptr->vars[varID].stdname == NULL )
    {
      stdname[0] = '\0';
    }  
  else
    strcpy(stdname, vlistptr->vars[varID].stdname);

  return;
}

/*
@Function  vlistInqVarUnits
@Title     Get the units of a Variable

@Prototype void vlistInqVarUnits(int vlistID, int varID, char *units)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  units    Units of the variable. The caller must allocate space for the 
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarUnits} returns the units of a variable if available,
otherwise the result is an empty string.

@Result
@func{vlistInqVarUnits} returns the units of the variable to the parameter units.

@EndFunction
*/
void vlistInqVarUnits(int vlistID, int varID, char *units)
{
  int tableID;
  int param;
  int pdis, pcat, pnum;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  units[0] = '\0';

  if ( vlistptr->vars[varID].units == NULL )
    {
      param = vlistptr->vars[varID].param;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      if ( pdis == 255 )
	{
	  int code = pnum;
	  tableID = vlistptr->vars[varID].tableID;
	  if ( tableInqParUnits(tableID, code, units) != 0 )
	    units[0] = '\0';
	}
    }  
  else
    strcpy(units, vlistptr->vars[varID].units);

  return;
}

/* not used 
int vlistInqVarID(int vlistID, int code)
{
  int varID;
  vlist_t *vlistptr;
  int param, pdis, pcat, pnum;

  vlistptr = vlist_to_pointer(vlistID);

  for ( varID = 0; varID < vlistptr->nvars; varID++ )
    {
      param = vlistptr->vars[varID].param;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      if ( pnum == code ) break;
    }

  if ( varID == vlistptr->nvars )
    {
      varID = CDI_UNDEFID;
    }

  return (varID);
}
*/

int vlistInqVarSize(int vlistID, int varID)
{
  int size;
  int zaxisID, gridID;
  int nlevs, gridsize;
  int timeID;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  vlistInqVar(vlistID, varID, &gridID, &zaxisID, &timeID);

  nlevs = zaxisInqSize(zaxisID);

  gridsize = gridInqSize(gridID);

  size = gridsize*nlevs;

  return (size);
}

/*
@Function  vlistInqVarDatatype
@Title     Get the data type of a Variable

@Prototype int vlistInqVarDatatype(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarDatatype} returns the data type of a variable.

@Result
@func{vlistInqVarDatatype} returns an identifier to the data type of the variable.
The valid CDI data types are @func{DATATYPE_PACK8}, @func{DATATYPE_PACK16}, @func{DATATYPE_PACK24},
@func{DATATYPE_FLT32}, @func{DATATYPE_FLT64}, @func{DATATYPE_INT8}, @func{DATATYPE_INT16} and 
@func{DATATYPE_INT32}.

@EndFunction
*/
int vlistInqVarDatatype(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].datatype);
}


int vlistInqVarNumber(int vlistID, int varID)
{
  vlist_t *vlistptr;
  int number = CDI_REAL;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( vlistptr->vars[varID].datatype == DATATYPE_CPX32 ||
       vlistptr->vars[varID].datatype == DATATYPE_CPX64 )
    number = CDI_COMP;

  return (number);
}

/*
@Function  vlistDefVarDatatype
@Title     Define the data type of a Variable

@Prototype void vlistDefVarDatatype(int vlistID, int varID, int datatype)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  datatype The data type identifier.
                    The valid CDI data types are @func{DATATYPE_PACK8}, @func{DATATYPE_PACK16},
                    @func{DATATYPE_PACK24}, @func{DATATYPE_FLT32}, @func{DATATYPE_FLT64},
                    @func{DATATYPE_INT8}, @func{DATATYPE_INT16} and @func{DATATYPE_INT32}.

@Description
The function @func{vlistDefVarDatatype} defines the data type of a variable.

@EndFunction
*/
void vlistDefVarDatatype(int vlistID, int varID, int datatype)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  vlistptr->vars[varID].datatype = datatype;
  
  if ( vlistptr->vars[varID].missvalused == FALSE )
    switch (datatype)
      {
      case DATATYPE_INT8:   vlistptr->vars[varID].missval = SCHAR_MIN; break;
      case DATATYPE_UINT8:  vlistptr->vars[varID].missval = UCHAR_MAX; break;
      case DATATYPE_INT16:  vlistptr->vars[varID].missval = SHRT_MIN;  break;
      case DATATYPE_UINT16: vlistptr->vars[varID].missval = USHRT_MAX; break;
      case DATATYPE_INT32:  vlistptr->vars[varID].missval = INT_MIN;   break;
      case DATATYPE_UINT32: vlistptr->vars[varID].missval = UINT_MAX;  break;
      }
}


void vlistDefVarInstitut(int vlistID, int varID, int instID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistptr->vars[varID].instID = instID;
}


int vlistInqVarInstitut(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->vars[varID].instID);
}


void vlistDefVarModel(int vlistID, int varID, int modelID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistptr->vars[varID].modelID = modelID;
}


int vlistInqVarModel(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->vars[varID].modelID);
}


void vlistDefVarTable(int vlistID, int varID, int tableID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistptr->vars[varID].tableID = tableID;

  {
    int param, pnum, pcat, pdis;
    int tablenum;
    tablenum = tableInqNum(tableID);

    param = vlistptr->vars[varID].param;

    cdiDecodeParam(param, &pnum, &pcat, &pdis);
  
    vlistptr->vars[varID].param = cdiEncodeParam(pnum, tablenum, pdis);
  }
}


int vlistInqVarTable(int vlistID, int varID)
{
  int tableID;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  tableID = vlistptr->vars[varID].tableID;

  return (tableID);
}

/*
@Function  vlistDefVarName
@Title     Define the name of a Variable

@Prototype void vlistDefVarName(int vlistID, int varID, const char *name)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  name     Name of the variable.

@Description
The function @func{vlistDefVarName} defines the name of a variable.

@EndFunction
*/
void vlistDefVarName(int vlistID, int varID, const char *name)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( name )
    {
      if ( vlistptr->vars[varID].name )
	{
	  free(vlistptr->vars[varID].name);
	  vlistptr->vars[varID].name = NULL;
	}

      vlistptr->vars[varID].name = strdupx(name);
    }
}

/*
@Function  vlistDefVarLongname
@Title     Define the long name of a Variable

@Prototype void vlistDefVarLongname(int vlistID, int varID, const char *longname)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  longname Long name of the variable.

@Description
The function @func{vlistDefVarLongname} defines the long name of a variable.

@EndFunction
*/
void vlistDefVarLongname(int vlistID, int varID, const char *longname)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( longname )
    {
      if ( vlistptr->vars[varID].longname )
	{
	  free(vlistptr->vars[varID].longname);
	  vlistptr->vars[varID].longname = 0;
	}

      vlistptr->vars[varID].longname = strdupx(longname);
    }
}

/*
@Function  vlistDefVarStdname
@Title     Define the standard name of a Variable

@Prototype void vlistDefVarStdname(int vlistID, int varID, const char *stdname)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  stdname  Standard name of the variable.

@Description
The function @func{vlistDefVarStdname} defines the standard name of a variable.

@EndFunction
*/
void vlistDefVarStdname(int vlistID, int varID, const char *stdname)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( stdname )
    {
      if ( vlistptr->vars[varID].stdname )
	{
	  free(vlistptr->vars[varID].stdname);
	  vlistptr->vars[varID].stdname = 0;
	}

      vlistptr->vars[varID].stdname = strdupx(stdname);
    }
}

/*
@Function  vlistDefVarUnits
@Title     Define the units of a Variable

@Prototype void vlistDefVarUnits(int vlistID, int varID, const char *units)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  units    Units of the variable.

@Description
The function @func{vlistDefVarUnits} defines the units of a variable.

@EndFunction
*/
void vlistDefVarUnits(int vlistID, int varID, const char *units)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( units )
    {
      if ( vlistptr->vars[varID].units )
	{
	  free(vlistptr->vars[varID].units);
	  vlistptr->vars[varID].units = 0;
	}

      vlistptr->vars[varID].units = strdupx(units);
    }
}

/*
@Function  vlistInqVarMissval
@Title     Get the missing value of a Variable

@Prototype double vlistInqVarMissval(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarMissval} returns the missing value of a variable.

@Result
@func{vlistInqVarMissval} returns the missing value of the variable.

@EndFunction
*/
double vlistInqVarMissval(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].missval);
}

/*
@Function  vlistDefVarMissval
@Title     Define the missing value of a Variable

@Prototype void vlistDefVarMissval(int vlistID, int varID, double missval)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  missval  Missing value.

@Description
The function @func{vlistDefVarMissval} defines the missing value of a variable.

@EndFunction
*/
void vlistDefVarMissval(int vlistID, int varID, double missval)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  vlistptr->vars[varID].missval = missval;
  vlistptr->vars[varID].missvalused = TRUE;
}


int vlistInqVarValidrange(int vlistID, int varID, double *validrange)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( validrange != NULL && vlistptr->vars[varID].lvalidrange )
    {
      validrange[0] = vlistptr->vars[varID].validrange[0];
      validrange[1] = vlistptr->vars[varID].validrange[1];
    }

  return (vlistptr->vars[varID].lvalidrange);
}


void vlistDefVarValidrange(int vlistID, int varID, const double *validrange)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  vlistptr->vars[varID].validrange[0] = validrange[0];
  vlistptr->vars[varID].validrange[1] = validrange[1];
  vlistptr->vars[varID].lvalidrange = TRUE;
}


double vlistInqVarScalefactor(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].scalefactor);
}


double vlistInqVarAddoffset(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].addoffset);
}

void vlistDefVarScalefactor(int vlistID, int varID, double scalefactor)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  vlistptr->vars[varID].scalefactor = scalefactor;
}


void vlistDefVarAddoffset(int vlistID, int varID, double addoffset)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  vlistptr->vars[varID].addoffset = addoffset;
}


void vlistDefVarTsteptype(int vlistID, int varID, int tsteptype)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistptr->vars[varID].tsteptype = tsteptype;
}


int vlistInqVarTsteptype(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->vars[varID].tsteptype);
}


void vlistDefVarTimave(int vlistID, int varID, int timave)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistptr->vars[varID].timave = timave;
}


int vlistInqVarTimave(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->vars[varID].timave);
}


void vlistDefVarTimaccu(int vlistID, int varID, int timaccu)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistptr->vars[varID].timaccu = timaccu;
}


int vlistInqVarTimaccu(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->vars[varID].timaccu);
}


void vlistDestroyVarName(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  if ( vlistptr->vars[varID].name )
    {
      free(vlistptr->vars[varID].name);
      vlistptr->vars[varID].name = NULL;
    }
}


void vlistDestroyVarLongname(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  if ( vlistptr->vars[varID].longname )
    {
      free(vlistptr->vars[varID].longname);
      vlistptr->vars[varID].longname = NULL;
    }
}


void vlistDestroyVarStdname(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  if ( vlistptr->vars[varID].stdname )
    {
      free(vlistptr->vars[varID].stdname);
      vlistptr->vars[varID].stdname = NULL;
    }
}


void vlistDestroyVarUnits(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  if ( vlistptr->vars[varID].units )
    {
      free(vlistptr->vars[varID].units);
      vlistptr->vars[varID].units = NULL;
    }
}


void vlistDefVarTime(int vlistID, int varID, int timeID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistptr->vars[varID].timeID = timeID;
}


int vlistInqVarMissvalUsed(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->vars[varID].missvalused);
}


void vlistDefFlag(int vlistID, int varID, int levID, int flag)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistptr->vars[varID].flag = flag;
  vlistptr->vars[varID].levinfo[levID].flag = flag;
}


int vlistInqFlag(int vlistID, int varID, int levID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->vars[varID].levinfo[levID].flag);
}


int vlistFindVar(int vlistID, int fvarID)
{
  int varID;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  for ( varID = 0; varID < vlistptr->nvars; varID++ )
    {
      if ( vlistptr->vars[varID].fvarID == fvarID ) break;
    }

  if ( varID == vlistptr->nvars )
    {
      varID = -1;
      Message("varID not found for fvarID %d in vlistID %d!", fvarID, vlistID);
    }

  return (varID);  
}


int vlistFindLevel(int vlistID, int fvarID, int flevelID)
{
  int varID;
  int levelID = -1;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  varID = vlistFindVar(vlistID, fvarID);

  if ( varID != -1 )
    {
      for ( levelID = 0; levelID < vlistptr->vars[varID].nlevs; levelID++ )
	{
	  if ( vlistptr->vars[varID].levinfo[levelID].flevelID == flevelID ) break;
	}

      if ( levelID == vlistptr->vars[varID].nlevs )
	{
	  levelID = -1;
	  Message("levelID not found for fvarID %d and levelID %d in vlistID %d!",
		  fvarID, flevelID, vlistID);
	}
    }

  return (levelID);  
}


int vlistMergedVar(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->vars[varID].mvarID);  
}


int vlistMergedLevel(int vlistID, int varID, int levelID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->vars[varID].levinfo[levelID].mlevelID);  
}


void vlistDefIndex(int vlistID, int varID, int levelID, int index)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistptr->vars[varID].levinfo[levelID].index = index;  
}


int vlistInqIndex(int vlistID, int varID, int levelID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->vars[varID].levinfo[levelID].index);  
}


void vlistChangeVarZaxis(int vlistID, int varID, int zaxisID)
{
  int nlevs1, nlevs2;
  int nvars, index;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  nlevs1 = zaxisInqSize(vlistptr->vars[varID].zaxisID);
  nlevs2 = zaxisInqSize(zaxisID);

  if ( nlevs1 != nlevs2 ) Error("Number of levels must not change!");

  nvars = vlistptr->nvars;
  for ( index = 0; index < nvars; index++ )
    if ( index != varID )
      if ( vlistptr->vars[index].zaxisID == vlistptr->vars[varID].zaxisID ) break;

  if ( index == nvars )
    {
      for ( index = 0; index < vlistptr->nzaxis; index++ )
	if ( vlistptr->zaxisIDs[index] == vlistptr->vars[varID].zaxisID )
	  vlistptr->zaxisIDs[index] = zaxisID;
    }
  else
    {
      for ( index = 0; index < vlistptr->nzaxis; index++ )
	if ( vlistptr->zaxisIDs[index] == zaxisID ) break;

      if ( index == vlistptr->nzaxis )
	{
	  if ( vlistptr->nzaxis + 1 >= MAX_ZAXES_PS )
	    Error("Maximum of %d zaxis reached", MAX_ZAXES_PS);

	  vlistptr->zaxisIDs[vlistptr->nzaxis] = zaxisID;
	  vlistptr->nzaxis++;
	}
    }
  
  vlistptr->vars[varID].zaxisID = zaxisID;
}


void vlistChangeVarGrid(int vlistID, int varID, int gridID)
{
  int nvars, index;
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  nvars = vlistptr->nvars;
  for ( index = 0; index < nvars; index++ )
    if ( index != varID )
      if ( vlistptr->vars[index].gridID == vlistptr->vars[varID].gridID ) break;

  if ( index == nvars )
    {
      for ( index = 0; index < vlistptr->ngrids; index++ )
	if ( vlistptr->gridIDs[index] == vlistptr->vars[varID].gridID )
	  vlistptr->gridIDs[index] = gridID;
    }
  else
    {
      for ( index = 0; index < vlistptr->ngrids; index++ )
	if ( vlistptr->gridIDs[index] == gridID ) break;

      if ( index == vlistptr->ngrids )
	{
	  if ( vlistptr->ngrids + 1 >= MAX_GRIDS_PS )
	    Error("Maximum of %d grids reached", MAX_GRIDS_PS);

	  vlistptr->gridIDs[vlistptr->ngrids] = gridID;
	  vlistptr->ngrids++;
	}
    }
  
  vlistptr->vars[varID].gridID = gridID;
}


void vlistDefVarCompType(int vlistID, int varID, int comptype)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  vlistptr->vars[varID].comptype = comptype;
}


int vlistInqVarCompType(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].comptype);
}


void vlistDefVarCompLevel(int vlistID, int varID, int complevel)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  vlistptr->vars[varID].complevel = complevel;
}


int vlistInqVarCompLevel(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].complevel);
}


void  vlistDefVarXYZ(int vlistID, int varID, int xyz)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  /* check xyz dimension order */
  {
    int dimorder[3];
    int dimx = 0, dimy = 0, dimz = 0;
    int posx = -1, posy = -1, posz = -1;
    dimorder[0] = xyz/100;
    dimorder[1] = (xyz-dimorder[0]*100)/10;
    dimorder[2] = (xyz-dimorder[0]*100-dimorder[1]*10);
    for ( int id = 0; id < 3; ++id )
      {
        if      ( dimorder[id] == 3 ) { dimz++; posz=id; }
        else if ( dimorder[id] == 2 ) { dimy++; posy=id; }
        else if ( dimorder[id] == 1 ) { dimx++; posx=id; }
      }
    if ( dimz > 1 || dimy > 1 || dimx > 1 ) xyz = 321; // ZYX
    else
      {
        int lchanged = 0;
        if ( dimz == 0 ) for ( int id = 0; id < 3; ++id ) if ( dimorder[id] == 0 ) {dimorder[id] = 3; lchanged++; break;}
        if ( dimy == 0 ) for ( int id = 0; id < 3; ++id ) if ( dimorder[id] == 0 ) {dimorder[id] = 2; lchanged++; break;}
        if ( dimx == 0 ) for ( int id = 0; id < 3; ++id ) if ( dimorder[id] == 0 ) {dimorder[id] = 1; lchanged++; break;}
        if ( lchanged ) xyz = dimorder[0]*100 + dimorder[1]*10 + dimorder[2];
      }
  }

  vlistptr->vars[varID].xyz = xyz;
}


int vlistInqVarXYZ(int vlistID, int varID)
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].xyz);
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
