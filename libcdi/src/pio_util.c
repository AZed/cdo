#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "pio_util.h"
#include "cdi.h"

char commands[][13] = { "FINALIZE\0",
                        "RESOURCES\0",
                        "WINCREATE\0",
                        "WRITETS\0"};


void pcdiAssert   ( bool assumption, const char * filename,
                    const char * functionname, int line )
{
  if ( !assumption )
    {
#ifdef USE_MPI
      int rank;

      MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
      fprintf ( stderr, "ERROR, ASSUMPTION FALSE: PE%d in %s, %s, line %d\n",
                rank, functionname, filename, line );
      fflush ( stderr );
      MPI_Abort ( MPI_COMM_WORLD, 1 );
#else
      fprintf ( stderr, "ERROR ASSUMPTION FALSE, %s, %s, line %d\n",
                functionname, filename, line );
      fflush ( stderr );
      abort();
#endif
    }
}

/****************************************************/

void pcdiAbort ( char * errorString, const char * filename,
		 const char *functionname, int line )
{
#ifdef USE_MPI
  int rank;

  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  fprintf ( stderr, "ERROR, pe%d in %s, %s, line %d, errorString: \"%s\"\n",
	    rank, functionname, filename, line, errorString );
  MPI_Abort ( MPI_COMM_WORLD, 1 );
#else
  fprintf ( stderr, "ERROR, %s, %s, line %d, errorString: \"%s\"\n",
            functionname, filename, line, errorString );
  abort();
#endif
}

/*****************************************************************************/

void * pcdiXmalloc ( size_t size, const char *filename, const char *functionname,
		     int line )
{
  void * value = calloc (1, size );

  if ( value == NULL )
    pcdiAbort ( "malloc failed", filename, functionname, line );

  return value;
}

void * pcdiXcalloc ( size_t nmemb, size_t size, const char *filename,
		     const char *functionname, int line )
{
  void * value = calloc ( nmemb, size );

  if ( value == NULL )
    pcdiAbort ( "calloc failed", filename, functionname, line );

  return value;
}

void * pcdiXrealloc ( void *p, size_t size, const char *functionname,
		      const char *filename, int line )
{
  void * value = realloc ( p, size );

  if ( value == NULL )
    pcdiAbort ( "realloc failed", filename, functionname, line );

  return value;
}

/***************************************************************/

#ifdef USE_MPI
void pcdiXMPI ( int iret, const char *filename, int line )
{
  char errorString1[MPI_MAX_ERROR_STRING + 1];
  char errorString2[MPI_MAX_ERROR_STRING + 1];
  int len, errorClass, rank;

  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  if ( iret != MPI_SUCCESS )
    {
      MPI_Error_class ( iret, &errorClass );
      MPI_Error_string ( errorClass, errorString1, &len );
      errorString1[len] = '\0';
      MPI_Error_string ( iret, errorString2, &len);
      errorString2[len] = '\0';

      fprintf ( stderr, "MPI ERROR, pe%d, %s, line %d,"
                "errorClass: \"%s\""
                "errorString: \"%s\"\n",
                rank, filename, line,
                errorString1, errorString2);

      MPI_Abort ( MPI_COMM_WORLD, iret );
    }
}

/*****************************************************************************/

void pcdiXMPIStat ( int iret, const char *filename, int line, MPI_Status *status )
{
  char errorString[MPI_MAX_ERROR_STRING + 1];
  int len, rank;

  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  if ( iret == MPI_ERR_IN_STATUS )
    {
      switch ( status->MPI_ERROR )
        {
          fprintf ( stderr, "------- checking error in request ----------\n" );
        case MPI_SUCCESS :
          fprintf ( stderr, "-------- mpi_success -----------\n" );
          break;
        case MPI_ERR_PENDING:
          fprintf ( stderr, "-------- mpi_err_pending ----------\n");
          break;
        default:
          MPI_Error_string ( status->MPI_ERROR, errorString, &len );
          errorString[len] = '\0';
          fprintf ( stderr,"MPI ERROR in request, pe%d, %s, line %d,"
                    "return value: %d, error_string: %s\n",
                    rank, filename, line, iret, errorString );
          MPI_Abort ( MPI_COMM_WORLD, iret );
        }
    }
  else
    xmpi ( iret );

  return;
}
#endif

/****************************************************/

#ifdef USE_MPI
void pcdiDebugComm ( const char *filename, const char *functionname, int line, MPI_Comm *comm )
{
  int rank, size, len, rankGlob;
  char *name;

  name = ( char * ) xmalloc ( MPI_MAX_OBJECT_NAME );
  memset ( name, 0, ( MPI_MAX_OBJECT_NAME ) * sizeof ( char ));
  MPI_Comm_get_name ( * comm, name, &len );
  MPI_Comm_size ( * comm, &size );
  MPI_Comm_rank ( * comm, &rank );
  MPI_Comm_rank ( MPI_COMM_WORLD, &rankGlob );
  fprintf ( stdout,
            "pe%d in %s, %s, line %d: comm: name=%s, size=%d, rank=%d\n",
            rankGlob, functionname, filename, line,
            name, size, rank );
  free ( name );

}
#endif

/****************************************************/

#ifdef USE_MPI
void pcdiDebugMsg ( const char * cdiPioDebugString, const char *filename,
                    const char *functionname, int line, int tag, int source,
                    int nfinished )
{
  int rank;

  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  fprintf ( stdout,
            "%s pe%d in %s, %s, line %d: command %s, source %d, finalized=%d\n",
            cdiPioDebugString, rank, functionname, filename, line,
            &commands[tag][0], source, nfinished );
}
#endif
/****************************************************/

#ifdef USE_MPI
void pcdiDebugMsg2 ( const char *filename, const char *functionname, int line,
                   int tag, int source, char * text )
{
  int rank;

  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  fprintf ( stdout,
            "pe%d in %s, %s, line %d: command %s, source %d, %s\n",
            rank, functionname, filename, line,
            &commands[tag][0], source, text );
}
#endif


/****************************************************/


int xmaxInt ( int a, int b )
{
  return a >= b ? a : b;
}


/****************************************************/


int xminInt ( int a, int b )
{
  return a <= b ? a : b;
}


/****************************************************/


int xsum ( int n, int * argarray )
{
  int i, sum = 0;

  for ( i = 0; i < n; i++ )
    sum += * ( argarray + i );

  return sum;
}


/****************************************************/


double xchecksum ( int type, int count, void * buffer )
{
  return 0.0;
}


/****************************************************/

void printArray ( const char * cdiPioDebugString, char * ps, const void * array, int n,
                  int datatype, const char * funname, const char * filename, int line )
{
  int i, rank;
  int * iArray;
  double * dArray;

#ifdef USE_MPI
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  fprintf ( stdout, "%s pe%d in %s, %s, line %d: %s = ",
	    cdiPioDebugString, rank, funname, filename, line, ps );
#else
  fprintf ( stdout, "%s %s, %s, line %d: %s = ",
	    cdiPioDebugString, funname, filename, line, ps );
#endif

  switch ( datatype )
    {
    case DATATYPE_INT:
      iArray = ( int * ) array;
      for ( i = 0; i < n-1; i++ )
	fprintf ( stdout, "%d ", * ( iArray + i ));
      fprintf ( stdout, "%d\n", * ( iArray + n - 1 ));
      break;
    case DATATYPE_FLT:
      dArray = ( double * ) array;
      for ( i = 0; i < n-1; i++ )
	fprintf ( stdout, "%.2f ", * ( dArray + i ));
      fprintf ( stdout, "%.2f\n", * ( dArray + n-1 ));
      break;
    default:
      fprintf ( stdout, " ... no datatype defined\n" );
    }

  return;
}

/****************************************************/
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
