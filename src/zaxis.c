/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2010 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "error.h"

#define UNDEFID -1

#define MAX_LINE_LEN 65536


typedef struct {
  double *vals;
  double *lbounds;
  double *ubounds;
  double *vct;
  int     vctsize;
  int     type;
  int     size;
  char    name[CDI_MAX_NAME];
  char    longname[CDI_MAX_NAME];
  char    units[CDI_MAX_NAME];
}
zaxis_t;


void zaxisInit(zaxis_t *zaxis)
{
  zaxis->vals        = NULL;
  zaxis->lbounds     = NULL;
  zaxis->ubounds     = NULL;
  zaxis->vct         = NULL;
  zaxis->type        = UNDEFID;
  zaxis->vctsize     = 0;
  zaxis->size        = 0;
  zaxis->name[0]     = 0;
  zaxis->longname[0] = 0;
  zaxis->units[0]    = 0;
}


static int getoptname(char *optname, const char *optstring, int nopt)
{
  int i, nerr = 0;
  size_t namelen;
  const char *pname;
  const char *pend;

  pname = optstring;
  pend  = optstring;

  for ( i = 0; i < nopt; i++ )
    {
      pend = strchr(pname, ',');
      if ( pend == NULL )
	break;
      else
	pname = pend + 1;
    }

  if ( pend )
    {
      pend = strchr(pname, ',');
      if ( pend == NULL )
	namelen = strlen(pname);
      else
	namelen = pend - pname;

      memcpy(optname, pname, namelen);
      optname[namelen] = '\0';
    }
  else
    nerr = 1;

  return (nerr);
}


int zaxisDefine(zaxis_t zaxis)
{
  int zaxisID = UNDEFID;

  if ( zaxis.type == -1 ) Error("zaxistype undefined!");

  if ( zaxis.size == 0 ) Error("zaxis size undefined!");

  zaxisID = zaxisCreate(zaxis.type, zaxis.size);

  if ( zaxis.vals )
    {
      zaxisDefLevels(zaxisID, zaxis.vals);
      free(zaxis.vals);
    }
  if ( zaxis.lbounds )
    {
      zaxisDefLbounds(zaxisID, zaxis.lbounds);
      free(zaxis.lbounds);
    }
  if ( zaxis.ubounds )
    {
      zaxisDefUbounds(zaxisID, zaxis.ubounds);
      free(zaxis.ubounds);
    }

  if ( zaxis.name[0] )     zaxisDefName(zaxisID, zaxis.name);
  if ( zaxis.longname[0] ) zaxisDefLongname(zaxisID, zaxis.longname);
  if ( zaxis.units[0] )    zaxisDefUnits(zaxisID, zaxis.units);

  if ( zaxis.type == ZAXIS_HYBRID || zaxis.type == ZAXIS_HYBRID_HALF )
    {
      if ( zaxis.vctsize && zaxis.vct )
	zaxisDefVct(zaxisID, zaxis.vctsize, zaxis.vct);
      else
	Warning("vct undefined!");	    
    }

  return (zaxisID);
}


static char *skipSeparator(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' ) pline++;
  while ( isspace((int) *pline) ) pline++;

  return (pline);
}


int zaxisFromFile(FILE *gfp)
{
  char line[MAX_LINE_LEN], *pline;
  int zaxisID;
  zaxis_t zaxis;

  zaxisInit(&zaxis);

  while ( readline(gfp, line, MAX_LINE_LEN) )
    {
      if ( line[0] == '#' ) continue;
      if ( line[0] == '\0' ) continue;
      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( pline[0] == '\0' ) continue;
      if ( memcmp(pline, "zaxistype", 9) == 0 || 
	   memcmp(pline, "type", 4) == 0 )
	{
	  if ( *pline == 'z' )
	    pline = skipSeparator(pline + 9);
	  else
	    pline = skipSeparator(pline + 4);

	  if ( memcmp(pline, "pressure", 6) == 0 )
	    zaxis.type = ZAXIS_PRESSURE;
	  else if ( memcmp(pline, "hybrid_half", 11)  == 0 )
	    zaxis.type = ZAXIS_HYBRID_HALF;
	  else if ( memcmp(pline, "hybrid", 6)  == 0 )
	    zaxis.type = ZAXIS_HYBRID;
	  else if ( memcmp(pline, "height", 6) == 0 )
	    zaxis.type = ZAXIS_HEIGHT;
	  else if ( memcmp(pline, "depth below sea", 15) == 0 ||
		    memcmp(pline, "depth_below_sea", 15) == 0 )
	    zaxis.type = ZAXIS_DEPTH_BELOW_SEA;
	  else if ( memcmp(pline, "depth below land", 16) == 0 ||
		    memcmp(pline, "depth_below_land", 16) == 0 )
	    zaxis.type = ZAXIS_DEPTH_BELOW_LAND;
	  else if ( memcmp(pline, "isentropic", 10)  == 0 )
	    zaxis.type = ZAXIS_ISENTROPIC;
	  else if ( memcmp(pline, "surface", 7)  == 0 )
	    zaxis.type = ZAXIS_SURFACE;
	  else
	    Warning("Invalid zaxisname : %s", pline);
	}
      else if ( memcmp(pline, "size", 4)  == 0 )
	{
	  zaxis.size = atol(skipSeparator(pline + 4));
	}
      else if ( memcmp(pline, "vctsize", 7)  == 0 )
	{
	  zaxis.vctsize = atol(skipSeparator(pline + 7));
	}
      else if ( memcmp(pline, "name", 4)  == 0 )
	{
	  strcpy(zaxis.name, skipSeparator(pline + 4));
	}
      else if ( memcmp(pline, "longname", 8)  == 0 )
	{
	  strcpy(zaxis.longname, skipSeparator(pline + 8));
	}
      else if ( memcmp(pline, "units", 5)  == 0 )
	{
	  strcpy(zaxis.units, skipSeparator(pline + 5));
	}
      else if ( memcmp(pline, "levels", 6)  == 0 )
	{
	  int i;
	  double flev;

	  if ( zaxis.size > 0 )
	    {
	      pline = skipSeparator(pline + 6);
	  
	      zaxis.vals = (double *) malloc(zaxis.size*sizeof(double));
	      for ( i = 0; i < zaxis.size; i++ )
		{
		  pline = skipSeparator(pline);
		  if ( strlen(pline) == 0 )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning("Incomplete command: >levels<");
			  break;
			}
		      pline = line;
		      pline = skipSeparator(pline);
		    }
		  flev = 0;
		  sscanf(pline, "%lg", &flev);
		  zaxis.vals[i] = flev;
		  while ( isalnum((int) *pline) ||
			  isdigit((int) *pline) ||
			  ispunct((int) *pline) ) pline++;
		}
	    }
	  else
	    {
	      Warning("size undefined!");
	    }
	}
      else if ( memcmp(pline, "vct", 3)  == 0 )
	{
	  int i;
	  double flev;

	  if ( zaxis.vctsize > 0 )
	    {
	      pline = skipSeparator(pline + 3);
	  
	      zaxis.vct = (double *) malloc(zaxis.vctsize*sizeof(double));
	      for ( i = 0; i < zaxis.vctsize; i++ )
		{
		  pline = skipSeparator(pline);
		  if ( strlen(pline) == 0 )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning("Incomplete command: >vct<");
			  break;
			}
		      pline = line;
		      pline = skipSeparator(pline);
		    }
		  flev = 0;
		  sscanf(pline, "%lg", &flev);
		  zaxis.vct[i] = flev;
		  while ( isalnum((int) *pline) ||
			  isdigit((int) *pline) ||
			  ispunct((int) *pline) ) pline++;
		}
	    }
	  else
	    {
	      Warning("vctsize undefined!");
	    }
	}
      else if ( memcmp(pline, "lbounds", 7)  == 0 )
	{
	  int i;
	  double flev;

	  if ( zaxis.size > 0 )
	    {
	      pline = skipSeparator(pline + 7);
	  
	      zaxis.lbounds = (double *) malloc(zaxis.size*sizeof(double));
	      for ( i = 0; i < zaxis.size; i++ )
		{
		  pline = skipSeparator(pline);
		  if ( strlen(pline) == 0 )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning("Incomplete command: >lbounds<");
			  break;
			}
		      pline = line;
		      pline = skipSeparator(pline);
		    }
		  flev = 0;
		  sscanf(pline, "%lg", &flev);
		  zaxis.lbounds[i] = flev;
		  while ( isalnum((int) *pline) ||
			  isdigit((int) *pline) ||
			  ispunct((int) *pline) ) pline++;
		}
	    }
	  else
	    {
	      Warning("size undefined!");
	    }
	}
      else if ( memcmp(pline, "ubounds", 7)  == 0 )
	{
	  int i;
	  double flev;

	  if ( zaxis.size > 0 )
	    {
	      pline = skipSeparator(pline + 7);
	  
	      zaxis.ubounds = (double *) malloc(zaxis.size*sizeof(double));
	      for ( i = 0; i < zaxis.size; i++ )
		{
		  pline = skipSeparator(pline);
		  if ( strlen(pline) == 0 )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning("Incomplete command: >ubounds<");
			  break;
			}
		      pline = line;
		      pline = skipSeparator(pline);
		    }
		  flev = 0;
		  sscanf(pline, "%lg", &flev);
		  zaxis.ubounds[i] = flev;
		  while ( isalnum((int) *pline) ||
			  isdigit((int) *pline) ||
			  ispunct((int) *pline) ) pline++;
		}
	    }
	  else
	    {
	      Warning("size undefined!");
	    }
	}
      else
	Warning("Invalid zaxis command : >%s<", pline);
    }

  zaxisID = zaxisDefine(zaxis);

  return (zaxisID);
}


int zaxisFromName(const char *zaxisname)
{
  const char *pline;
  int zaxisID = UNDEFID;
  zaxis_t zaxis;

  zaxisInit(&zaxis);

  pline = zaxisname;
  if ( memcmp(pline, "surface", 7) == 0 ) /* surface */
    {
      zaxis.type = ZAXIS_SURFACE;
      zaxis.size = 1;
      zaxis.vals = (double *) malloc(zaxis.size*sizeof(double));
      zaxis.vals[0] = 0;
    }

  if ( zaxis.type != -1 ) zaxisID = zaxisDefine(zaxis);

  return (zaxisID);
}


int cdoDefineZaxis(const char *zaxisfile)
{
  FILE *zfp;
  int zaxisID = -1;

  zfp = fopen(zaxisfile, "r");
  if ( zfp == NULL )
    {
      zaxisID = zaxisFromName(zaxisfile);

      if ( zaxisID == -1 ) cdoAbort("Open failed on %s!", zaxisfile);
    }
  else
    {
      zaxisID = zaxisFromFile(zfp);
      fclose(zfp);
    }

  if ( zaxisID == -1 ) cdoAbort("Invalid zaxis description file %s!", zaxisfile);

  return (zaxisID);
}


void defineZaxis(const char *zaxisarg)
{
  char zaxisfile[4096];
  int nfile = 0;

  while ( getoptname(zaxisfile, zaxisarg, nfile++) == 0 )
    {      
      (void) cdoDefineZaxis(zaxisfile);
    }
}


int zaxis2ltype(int zaxisID)
{
  int ltype;
  int zaxistype;

  zaxistype = zaxisInqType(zaxisID);

  ltype = zaxisInqLtype(zaxisID);

  if ( ltype <= 0 ) ltype = ztype2ltype(zaxistype);

  return (ltype);
}


int ztype2ltype(int zaxistype)
{
  int ltype;

  ltype = -1;
  if      ( zaxistype == ZAXIS_SURFACE           )  ltype =   1;
  else if ( zaxistype == ZAXIS_PRESSURE          )  ltype = 100;
  else if ( zaxistype == ZAXIS_ALTITUDE          )  ltype = 103;
  else if ( zaxistype == ZAXIS_HEIGHT            )  ltype = 105;
  else if ( zaxistype == ZAXIS_SIGMA             )  ltype = 107;
  else if ( zaxistype == ZAXIS_HYBRID            )  ltype = 109;
  else if ( zaxistype == ZAXIS_HYBRID_HALF       )  ltype = 109;
  else if ( zaxistype == ZAXIS_DEPTH_BELOW_LAND  )  ltype = 111;
  else if ( zaxistype == ZAXIS_ISENTROPIC        )  ltype = 113;
  else if ( zaxistype == ZAXIS_DEPTH_BELOW_SEA   )  ltype = 160;
  else cdoWarning("zaxis type %d not supported", zaxistype);

  return (ltype);
}


int ltype2ztype(int ltype)
{
  int zaxistype = -1;

  if      ( ltype ==   1 ) zaxistype = ZAXIS_SURFACE;
  else if ( ltype ==   0 ) zaxistype = ZAXIS_SURFACE;
  else if ( ltype == 100 ) zaxistype = ZAXIS_PRESSURE;
  else if ( ltype == 103 ) zaxistype = ZAXIS_ALTITUDE;
  else if ( ltype == 105 ) zaxistype = ZAXIS_HEIGHT;
  else if ( ltype == 107 ) zaxistype = ZAXIS_SIGMA;
  else if ( ltype == 109 ) zaxistype = ZAXIS_HYBRID;
  else if ( ltype == 110 ) zaxistype = ZAXIS_HYBRID_HALF;
  else if ( ltype == 111 ) zaxistype = ZAXIS_DEPTH_BELOW_LAND;
  else if ( ltype == 113 ) zaxistype = ZAXIS_ISENTROPIC;
  else if ( ltype == 160 ) zaxistype = ZAXIS_DEPTH_BELOW_SEA;
  else                     zaxistype = ZAXIS_GENERIC;

  return (zaxistype);
}
