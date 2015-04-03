#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "dmemory.h"

#include "cdi.h"
#include "stream_int.h"
#include "vlist.h"


static
cdi_atts_t *get_attsp(vlist_t *vlistptr, int varID)
{
  cdi_atts_t *attsp = NULL;

  if ( varID == CDI_GLOBAL )
    {
      attsp = &vlistptr->atts;
    }
  else
    {
      if ( varID >= 0 && varID < vlistptr->nvars )
	attsp = &(vlistptr->vars[varID].atts);
    }

  return (attsp);
}

static
cdi_att_t *find_att(cdi_atts_t *attsp, const char *name)
{
  cdi_att_t *attp;
  size_t attid;
  size_t slen;

  assert(attsp != NULL);

  if ( attsp->nelems == 0 ) return NULL;

  slen = strlen(name);

  for ( attid = 0; attid < attsp->nelems; attid++ )
    {
      attp = &(attsp->value[attid]);
      if ( attp->namesz == slen )
	if ( memcmp(attp->name, name, slen) == 0)
	  {
	    return (attp); /* Normal return */
	  }
    }

  return (NULL);
}

static
cdi_att_t *new_att(cdi_atts_t *attsp, const char *name)
{
  cdi_att_t *attp;
  size_t slen;

  assert(attsp != NULL);
  assert(name  != NULL);

  if ( attsp->nelems == attsp->nalloc ) return (NULL);

  attp = &(attsp->value[attsp->nelems]);
  attsp->nelems++;

  slen = strlen(name);

  attp->name = (char *) malloc(slen+1);
  memcpy(attp->name, name, slen+1);
  attp->namesz = slen;
  attp->xvalue = NULL;

  return (attp);
}

static
void fill_att(cdi_att_t *attp, int indtype, int exdtype, size_t nelems, size_t xsz, const void *xvalue)
{
  assert(attp != NULL);

  attp->xsz = xsz;
  attp->indtype = indtype;
  attp->exdtype = exdtype;
  attp->nelems  = nelems;

  if ( xsz > 0 )
    {
      attp->xvalue = (void *) realloc(attp->xvalue, xsz);
      memcpy(attp->xvalue, xvalue, xsz);
    }
}

/*
@Function  vlistInqNatts
@Title     Get number of variable attributes

@Prototype int vlistInqNatts(int vlistID, int varID, int *nattsp)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier, or CDI_GLOBAL for a global attribute.
    @Item  nattsp   Pointer to location for returned number of variable attributes.

@Description
The function @func{vlistInqNatts} gets the number of variable attributes assigned to this variable.

@EndFunction
*/
int vlistInqNatts(int vlistID, int varID, int *nattsp)
{
  int status = CDI_NOERR;
  vlist_t *vlistptr;
  cdi_atts_t *attsp;

  vlistptr = vlist_to_pointer(vlistID);
  
  attsp = get_attsp(vlistptr, varID);
  assert(attsp != NULL);

  *nattsp = attsp->nelems;

  return (status);
}

/*
@Function  vlistInqAtt
@Title     Get information about an attribute

@Prototype int vlistInqAtt(int vlistID, int varID, int attnum, char *name, int *typep, int *lenp)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier, or CDI_GLOBAL for a global attribute.
    @Item  attnum   Attribute number (from 0 to natts-1).
    @Item  name     Pointer to the location for the returned attribute name. The caller must allocate space for the 
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant CDI_MAX_NAME.
    @Item  typep    Pointer to location for returned attribute type.
    @Item  lenp     Pointer to location for returned attribute number.

@Description
The function @func{vlistInqNatts} gets information about an attribute.

@EndFunction
*/
int vlistInqAtt(int vlistID, int varID, int attnum, char *name, int *typep, int *lenp)
{
  int status = CDI_NOERR;
  vlist_t *vlistptr;
  cdi_att_t *attp = NULL;
  cdi_atts_t *attsp;

  assert(name != NULL);

  vlistptr = vlist_to_pointer(vlistID);

  attsp = get_attsp(vlistptr, varID);
  assert(attsp != NULL);

  if ( attnum >= 0 && attnum < (int)attsp->nelems )
    attp = &(attsp->value[attnum]);

  if ( attp != NULL ) /* name in use */
    {
      memcpy(name, attp->name, attp->namesz+1);
      *typep  = attp->exdtype;
      *lenp   = attp->nelems;
    }
  else
    {
      name[0] =  0;
      *typep  = -1;
      *lenp   =  0;
    }

  return (status);
}


int vlistDelAtts(int vlistID, int varID)
{
  int status = CDI_NOERR;
  vlist_t *vlistptr;
  cdi_att_t *attp = NULL;
  cdi_atts_t *attsp;
  int attid;

  vlistptr = vlist_to_pointer(vlistID);

  attsp = get_attsp(vlistptr, varID);
  assert(attsp != NULL);

  for ( attid = 0; attid < (int)attsp->nelems; attid++ )
    {
      attp = &(attsp->value[attid]);
      if ( attp->name   ) free(attp->name);
      if ( attp->xvalue ) free(attp->xvalue);
    }

  attsp->nelems = 0;

  return (status);
}


int vlistDelAtt(int vlistID, int varID, const char *name)
{
  int status = CDI_NOERR;

  fprintf(stderr, "vlistDelAtt not implemented!\n");

  return (status);
}

static
int vlist_def_att(int indtype, int exdtype, int vlistID, int varID, const char *name, size_t len, size_t xsz, const void *xp)
{
  int status = CDI_NOERR;
  vlist_t *vlistptr;
  cdi_att_t *attp;
  cdi_atts_t *attsp;

  if ( len != 0 && xp == NULL ) /* Null arg */
    {
      return (CDI_EINVAL);
    }

  vlistptr = vlist_to_pointer(vlistID);

  attsp = get_attsp(vlistptr, varID);
  assert(attsp != NULL);

  attp = find_att(attsp, name);
  if ( attp == NULL )
    attp = new_att(attsp, name);

  if ( attp != NULL )
    fill_att(attp, indtype, exdtype, len, xsz, xp);
  
  return (status);
}

static
int vlist_inq_att(int indtype, int vlistID, int varID, const char *name, size_t mxsz, void *xp)
{
  int status = CDI_NOERR;
  vlist_t *vlistptr;
  cdi_att_t *attp;
  cdi_atts_t *attsp;
  size_t xsz;

  if ( mxsz != 0 && xp == NULL ) /* Null arg */
    {
      return (CDI_EINVAL);
    }

  vlistptr = vlist_to_pointer(vlistID);

  attsp = get_attsp(vlistptr, varID);
  assert(attsp != NULL);

  attp = find_att(attsp, name);
  if ( attp != NULL ) /* name in use */
    {
      if ( attp->indtype == indtype )
	{
	  xsz = attp->xsz;
	  if ( mxsz < xsz ) xsz = mxsz;
	  if ( xsz > 0 )
	    memcpy(xp, attp->xvalue, xsz);
	}
      else
	{
	  Warning("Attribute %s has wrong data type!", name);
	}
    }
  else
    {
      Warning("Internal problem, attribute %s not found!", name);
    }

  return (status);
}


int vlistCopyVarAtts(int vlistID1, int varID_1, int vlistID2, int varID_2)
{
  int status = CDI_NOERR;
  vlist_t *vlistptr1;
  cdi_att_t *attp = NULL;
  cdi_atts_t *attsp1;
  int attid;

  vlistptr1 = vlist_to_pointer(vlistID1);

  attsp1 = get_attsp(vlistptr1, varID_1);
  assert(attsp1 != NULL);

  for ( attid = 0; attid < (int)attsp1->nelems; attid++ )
    {
      attp = &(attsp1->value[attid]);
      vlist_def_att(attp->indtype, attp->exdtype, vlistID2, varID_2, attp->name, attp->nelems, attp->xsz, attp->xvalue);
    }

  return (status);
}

/*
@Function  vlistDefAttInt
@Title     Define an integer attribute

@Prototype int vlistDefAttInt(int vlistID, int varID, const char *name, int type, int len, const int *ip)

@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier, or CDI_GLOBAL for a global attribute.
    @Item  name     Attribute name.
    @Item  type     External data type (DATATYPE_INT16 or DATATYPE_INT32).
    @Item  len      Number of values provided for the attribute.
    @Item  ip       Pointer to one or more integer values.

@Description
The function @func{vlistDefAttInt} defines an integer attribute.

@EndFunction
*/
int vlistDefAttInt(int vlistID, int varID, const char *name, int type, int len, const int *ip)
{
  int status;

  status = vlist_def_att(DATATYPE_INT, type, vlistID, varID, name, (size_t) len, len*sizeof(int), (const void *) ip);

  return (status);
}

/*
@Function  vlistDefAttFlt
@Title     Define a floating point attribute

@Prototype int vlistDefAttFlt(int vlistID, int varID, const char *name, int type, int len, const double *dp)

@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier, or CDI_GLOBAL for a global attribute.
    @Item  name     Attribute name.
    @Item  type     External data type (DATATYPE_FLT32 or DATATYPE_FLT64).
    @Item  len      Number of values provided for the attribute.
    @Item  dp       Pointer to one or more floating point values.

@Description
The function @func{vlistDefAttFlt} defines a floating point attribute.

@EndFunction
*/
int vlistDefAttFlt(int vlistID, int varID, const char *name, int type, int len, const double *dp)
{
  int status;

  status = vlist_def_att(DATATYPE_FLT, type, vlistID, varID, name, (size_t) len, len*sizeof(double), (const void *) dp);

  return (status);
}

/*
@Function  vlistDefAttTxt
@Title     Define a text attribute

@Prototype int vlistDefAttTxt(int vlistID, int varID, const char *name, int len, const char *tp)

@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier, or CDI_GLOBAL for a global attribute.
    @Item  name     Attribute name.
    @Item  len      Number of values provided for the attribute.
    @Item  tp       Pointer to one or more character values.

@Description
The function @func{vlistDefAttTxt} defines a text attribute.

@EndFunction
*/
int vlistDefAttTxt(int vlistID, int varID, const char *name, int len, const char *tp)
{
  int status;

  status = vlist_def_att(DATATYPE_TXT, DATATYPE_TXT, vlistID, varID, name, (size_t) len, len*sizeof(char), (const void *) tp);

  return (status);
}

/*
@Function  vlistInqAttInt
@Title     Get the value(s) of an integer attribute

@Prototype int vlistInqAttInt(int vlistID, int varID, const char *name, int mlen, int *ip)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier, or CDI_GLOBAL for a global attribute.
    @Item  name     Attribute name.
    @Item  mlen     Number of allocated values provided for the attribute.
    @Item  ip       Pointer location for returned integer attribute value(s).

@Description
The function @func{vlistInqAttInt} gets the values(s) of an integer attribute.

@EndFunction
*/
int vlistInqAttInt(int vlistID, int varID, const char *name, int mlen, int *ip)
{
  int status = CDI_NOERR;

  status = vlist_inq_att(DATATYPE_INT, vlistID, varID, name, mlen*sizeof(int), (void *) ip);

  return (CDI_NOERR);
}

/*
@Function  vlistInqAttFlt
@Title     Get the value(s) of a floating point attribute

@Prototype int vlistInqAttFlt(int vlistID, int varID, const char *name, int mlen, int *dp)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier, or CDI_GLOBAL for a global attribute.
    @Item  name     Attribute name.
    @Item  mlen     Number of allocated values provided for the attribute.
    @Item  dp       Pointer location for returned floating point attribute value(s).

@Description
The function @func{vlistInqAttFlt} gets the values(s) of a floating point attribute.

@EndFunction
*/
int vlistInqAttFlt(int vlistID, int varID, const char *name, int mlen, double *dp)
{
  int status = CDI_NOERR;

  status = vlist_inq_att(DATATYPE_FLT, vlistID, varID, name, mlen*sizeof(double), (void *) dp);

  return (status);
}

/*
@Function  vlistInqAttTxt
@Title     Get the value(s) of a text attribute

@Prototype int vlistInqAttTxt(int vlistID, int varID, const char *name, int mlen, int *tp)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier, or CDI_GLOBAL for a global attribute.
    @Item  name     Attribute name.
    @Item  mlen     Number of allocated values provided for the attribute.
    @Item  tp       Pointer location for returned text attribute value(s).

@Description
The function @func{vlistInqAttTxt} gets the values(s) of a text attribute.

@EndFunction
*/
int vlistInqAttTxt(int vlistID, int varID, const char *name, int mlen, char *tp)
{
  int status = CDI_NOERR;

  status = vlist_inq_att(DATATYPE_TXT, vlistID, varID, name, mlen*sizeof(char), (void *) tp);

  return (status);
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
