#ifndef _TABLEPAR_H
#define _TABLEPAR_H

typedef struct
{
  int   id;	     /* Parameter number (GRIB) */
  char *name;	     /* Parameter name */
  char *longname;    /* Parameter long name */
  char *units;	     /* Parameter units */
}
PAR;


void tableLink(int tableID, PAR *pars, int npars);
int tableDef(int modelID, int tablegribID, const char *tablename);

#endif
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
