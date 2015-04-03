#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "cdi.h"

char *Progname;

void version(void)
{
  fprintf(stderr, "createtable version 1.00 (11 Nov 2001)\n");
  fprintf(stderr, "\n");
}

void usage(void)
{
  fprintf(stderr, "usage :  createtable  [-options]  ifile  ofile\n");
  fprintf(stderr, "with:\n");

  fprintf(stderr, "   -debug        debugging mode\n");
  fprintf(stderr, "   -version      display version number\n");
}

int main(int argc, char *argv[])
{
  int c, len;
  char  *cstring;
  char *ifile = NULL, *ofile = NULL;
  int debug = 0;
  int tableID;


  Progname = argv[0];

  while ( --argc && (*++argv)[0] == '-' )
    {
      c = *++argv[0];
      cstring = argv[0];
      len = strlen(cstring);
      switch (c)
	{
        case 'd':
	  if ( !strncmp(cstring, "debug", len) )
	    {
	      debug = 1;
	    }
	  break;
        case 'h':
	  if ( !strncmp(cstring, "help", len) )
	    { 
	      usage( );
	      return( 0 );
            }
	  break;
        case 'v':
	  if ( !strncmp(cstring, "version", len) )
	    {
	      version();
	      return (0);
            }
	  break;
        default:
	  usage();
	  fprintf(stderr, "illegal option %s\n", cstring);
	  return (-1);
	  break;
        }
    }

  if ( argc )
    {
      ifile = argv[0];
      argc--;
    }
  if ( argc )
    {
      ofile = (++argv)[0];
      argc--;
    }

  if ( ifile == NULL || ofile == NULL )
    {
      usage();
      fprintf(stderr, "missing filenames\n");
      return (-1);
    }

  if ( debug )
    {
      fprintf(stderr, "\n");
      if ( ifile )
	fprintf(stderr, "  < %s  \n", ifile);
      if ( ofile )
	fprintf(stderr, "  > %s\n\n", ofile);
    }
  /*
  if ( debug ) cdiDebug(debug);
  */
  tableID = tableRead(ifile);
  if ( tableID != CDI_UNDEFID )
    tableWriteC(ofile, tableID);

  return (0);
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
