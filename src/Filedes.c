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

/*
   This module contains the following operators:

      Filedes    pardes          Parameter description
      Filedes    griddes         Grid description
      Filedes    vct             Vertical coordinate table
*/

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


static
void printAtts(int vlistID, int varID)
{
#define MAXATT 8192
  int natts, ia;
  char attname[1024];
  int atttype, attlen;
  char atttxt[MAXATT];
  int attint[MAXATT];
  double attflt[MAXATT];
  int i;

  vlistInqNatts(vlistID, varID, &natts);

  for ( ia = 0; ia < natts; ++ia )
    {
      vlistInqAtt(vlistID, varID, ia, attname, &atttype, &attlen);
      if ( atttype == DATATYPE_INT )
	{
	  if ( attlen > MAXATT ) attlen = MAXATT;
	  vlistInqAttInt(vlistID, varID, attname, attlen, attint);
	  fprintf(stdout, "  %s=", attname);
	  for ( i = 0; i < attlen; ++i)
	    {
	      if ( i > 0 ) fprintf(stdout, ", ");
	      fprintf(stdout, "%d", attint[i]);
	    }
	  fprintf(stdout, "\n");
	}
      else if ( atttype == DATATYPE_FLT )
	{
	  if ( attlen > MAXATT ) attlen = MAXATT;
	  vlistInqAttFlt(vlistID, varID, attname, MAXATT, attflt);
	  fprintf(stdout, "  %s=", attname);
	  for ( i = 0; i < attlen; ++i)
	    {
	      if ( i > 0 ) fprintf(stdout, ", ");
	      fprintf(stdout, "%g", attflt[i]);
	    }
	  fprintf(stdout, "\n");
	}
      else if ( atttype == DATATYPE_TXT )
	{
	  vlistInqAttTxt(vlistID, varID, attname, sizeof(atttxt), atttxt);
	  atttxt[attlen] = 0;
	  fprintf(stdout, "  %s=\"%s\"\n", attname, atttxt);
	}
    }
}


void *Filedes(void *argument)
{
  int GRIDDES, GRIDDES2, ZAXISDES, VCT, VCT2, PARDES, FILEDES, VLIST, PARTAB, PARTAB2;
  int operatorID;
  int streamID = 0;
  int zaxisID;
  int nvars, ngrids, nzaxis;
  int type, index;
  int vlistID;

  cdoInitialize(argument);

  GRIDDES  = cdoOperatorAdd("griddes",   0, 0, NULL);
  GRIDDES2 = cdoOperatorAdd("griddes2",  0, 0, NULL);
  ZAXISDES = cdoOperatorAdd("zaxisdes",  0, 0, NULL);
  VCT      = cdoOperatorAdd("vct",       0, 0, NULL);
  VCT2     = cdoOperatorAdd("vct2",      0, 0, NULL);
  PARDES   = cdoOperatorAdd("pardes",    0, 0, NULL);
  FILEDES  = cdoOperatorAdd("filedes",   0, 0, NULL);
  VLIST    = cdoOperatorAdd("vlist",     0, 0, NULL);
  PARTAB   = cdoOperatorAdd("partab",    0, 0, NULL);
  PARTAB2  = cdoOperatorAdd("partab2",   0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID = streamOpenRead(cdoStreamName(0));
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));

  vlistID = streamInqVlist(streamID);

  nvars  = vlistNvars(vlistID);
  ngrids = vlistNgrids(vlistID);
  nzaxis = vlistNzaxis(vlistID);

  if ( operatorID == GRIDDES || operatorID == GRIDDES2 )
    {
      int opt = 0;
      if ( operatorID == GRIDDES ) opt = 1;
      for ( index = 0; index < ngrids; index++ )
	gridPrint(vlistGrid(vlistID, index), opt);
    }
  else if ( operatorID == ZAXISDES )
    {
      for ( index = 0; index < nzaxis; index++ )
	zaxisPrint(vlistZaxis(vlistID, index));
    }
  else if ( operatorID == VCT || operatorID == VCT2 )
    {
      for ( index = 0; index < nzaxis; index++)
	{
	  zaxisID = vlistZaxis(vlistID, index);
	  type = zaxisInqType(zaxisID);
	  if ( type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF )
	    {
	      int i, vctsize;
	      const double *vct;

	      vctsize = zaxisInqVctSize(zaxisID);
	      vct     = zaxisInqVctPtr(zaxisID);
		
	      if ( vctsize%2 == 0 )
		{
		  if ( operatorID == VCT )
		    {
		      for ( i = 0; i < vctsize/2; i++ )
			fprintf(stdout, "%5d %25.17f %25.17f\n", i, vct[i], vct[vctsize/2+i]);
		    }
		  else
		    {
		      int nbyte0, nbyte;
		      fprintf(stdout, "vctsize   = %d\n", vctsize);
		      nbyte0 = fprintf(stdout, "vct       = ");
		      nbyte = nbyte0;
		      for ( i = 0; i < vctsize; i++ )
			{
			  if ( nbyte > 70 || i == vctsize/2 )
			    {
			      fprintf(stdout, "\n%*s", nbyte0, "");
			      nbyte = nbyte0;
			    }
			  nbyte += fprintf(stdout, "%.9g ", vct[i]);
			}
		      fprintf(stdout, "\n");
		    }
		}
	      else
		for ( i = 0; i < vctsize; i++ )
		  fprintf(stdout, "%5d %25.17f\n", i, vct[i]);

	      break;
	    }
	}
    }
  else if ( operatorID == VLIST )
    {
      vlistPrint(vlistID);
    }
  else if ( operatorID == PARDES )
    {
      int varID, code;
      char varname[128], varlongname[128], varunits[128];

      for ( varID = 0; varID < nvars; varID++ )
	{
	  varname[0]     = 0;
	  varlongname[0] = 0;
	  varunits[0]    = 0;
	  code     = vlistInqVarCode(vlistID, varID);
	  vlistInqVarName(vlistID, varID, varname);
	  vlistInqVarLongname(vlistID, varID, varlongname);
	  vlistInqVarUnits(vlistID, varID, varunits);
	  fprintf(stdout, "%4d  %-12s", code, varname);
	  if ( strlen(varlongname) )
	    {
	      fprintf(stdout, "  %s", varlongname);
	      if ( strlen(varunits) )
		fprintf(stdout, " [%s]", varunits);
	    }
	  fprintf(stdout, "\n");
	}   
    }
  else if ( operatorID == PARTAB || operatorID == PARTAB2 )
    {
      int varID, code, tabnum, tableID, prec;
      char pstr[4];
      char varname[128], varlongname[128], varstdname[128], varunits[128];
      int natts;
      
      if (  operatorID == PARTAB2 )
	{
	  vlistInqNatts(vlistID, CDI_GLOBAL, &natts);
	  if ( natts > 0 )
	    {
	      fprintf(stdout, "&PARAMETER\n");
	      fprintf(stdout, "  NAME=_GLOBAL_\n");
	      printAtts(vlistID, CDI_GLOBAL);
	      fprintf(stdout, "/\n");
	    }
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  fprintf(stdout, "&PARAMETER\n");

	  varname[0]     = 0;
	  varlongname[0] = 0;
	  varunits[0]    = 0;
	  code     = vlistInqVarCode(vlistID, varID);
	  tableID  = vlistInqVarTable(vlistID, varID);
	  tabnum   = tableInqNum(tableID);
	  vlistInqVarName(vlistID, varID, varname);
	  /* printf("1>%s<\n", varname); */
	  vlistInqVarStdname(vlistID, varID, varstdname);
	  /* printf("2>%s<\n", varname); */
	  vlistInqVarLongname(vlistID, varID, varlongname);
	  /* printf("3>%s<\n", varname); */
	  vlistInqVarUnits(vlistID, varID, varunits);

	  prec = vlistInqVarDatatype(vlistID, varID);
	  if      ( prec == DATATYPE_PACK   ) strcpy(pstr, "P0");
	  else if ( prec > 0 && prec <= 32  ) sprintf(pstr, "P%d", prec);
	  else if ( prec == DATATYPE_FLT32  ) strcpy(pstr, "F32");
	  else if ( prec == DATATYPE_FLT64  ) strcpy(pstr, "F64");
	  else if ( prec == DATATYPE_INT8   ) strcpy(pstr, "I8");
	  else if ( prec == DATATYPE_INT16  ) strcpy(pstr, "I16");
	  else if ( prec == DATATYPE_INT32  ) strcpy(pstr, "I32");
	  else if ( prec == DATATYPE_UINT8  ) strcpy(pstr, "U8");
	  else if ( prec == DATATYPE_UINT16 ) strcpy(pstr, "U16");
	  else if ( prec == DATATYPE_UINT32 ) strcpy(pstr, "U32");
	  else                                strcpy(pstr, "-1");

	  if ( code   > 0 ) fprintf(stdout, "  CODE=%d\n", code);
	  if ( tabnum > 0 ) fprintf(stdout, "  TABLE=%d\n", tabnum);
	  fprintf(stdout, "  NAME=%s\n", varname);
	  if ( strlen(varstdname) )
	    fprintf(stdout, "  STANDARD_NAME=%s\n", varstdname);
	  if ( strlen(varlongname) )
	    fprintf(stdout, "  LONG_NAME=\"%s\"\n", varlongname);
	  if ( strlen(varunits) )
	    fprintf(stdout, "  UNITS=\"%s\"\n", varunits);

	  /* if ( pstr ) fprintf(stdout, "  DATATYPE=%s\n", pstr); */

	  if ( operatorID == PARTAB2 ) printAtts(vlistID, varID);

	  fprintf(stdout, "/\n");
	}   
    }
  else if ( operatorID == FILEDES )
    {
      int filetype;

      printf("\n");
      filetype = streamInqFiletype(streamID);
      switch ( filetype )
	{
	case FILETYPE_GRB:
	  printf("  GRIB data\n");
	  break;
	case FILETYPE_NC:
	  printf("  netCDF data\n");
	  break;
	case FILETYPE_NC2:
	  printf("  netCDF2 data\n");
	  break;
	case FILETYPE_NC4:
	  printf("  netCDF4 data\n");
	  break;
	case FILETYPE_SRV:
	  printf("  SERVICE data\n");
	  switch ( streamInqByteorder(streamID) )
	    {
	    case CDI_BIGENDIAN:
	      printf("  byteorder is BIGENDIAN\n"); break;
	    case CDI_LITTLEENDIAN:
	      printf("  byteorder is LITTLEENDIAN\n"); break;
	    default:
	      printf("  byteorder %d undefined\n", streamInqByteorder(streamID)); break;
	    }
	   break;
	case FILETYPE_EXT:
	  printf("  EXTRA data\n");
	  switch ( streamInqByteorder(streamID) )
	    {
	    case CDI_BIGENDIAN:
	      printf("  byteorder is BIGENDIAN\n"); break;
	    case CDI_LITTLEENDIAN:
	      printf("  byteorder is LITTLEENDIAN\n"); break;
	    default:
	      printf("  byteorder %d undefined\n", streamInqByteorder(streamID)); break;
	    }
	   break;
	case FILETYPE_IEG:
	  printf("  IEG data\n");
	  switch ( streamInqByteorder(streamID) )
	    {
	    case CDI_BIGENDIAN:
	      printf("  byteorder is BIGENDIAN\n"); break;
	    case CDI_LITTLEENDIAN:
	      printf("  byteorder is LITTLEENDIAN\n"); break;
	    default:
	      printf("  byteorder %d undefined\n", streamInqByteorder(streamID)); break;
	    }
	   break;
	default:
	  printf("  unsupported filetype %d\n" , filetype);
	}

      printf("\n");
    }

  streamClose(streamID);

  cdoFinish();

  return (0);
}
