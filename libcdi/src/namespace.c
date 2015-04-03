#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include "cdi.h"
#include "namespace.h"
#include "pio_util.h"


static int nNamespaces = 1;
static int activeNamespace = 0;
static int serialHLF = 1;
static int * hasLocalFiles = &serialHLF;
static int serialRS = STAGE_DEFINITION;
static statusCode * resStatus = (statusCode *) &serialRS;

enum {
  intbits = sizeof(int) * CHAR_BIT,
  nspbits = 4,
  idxbits = intbits - nspbits,
  nspmask = (( 1 << nspbits ) - 1) << idxbits,
  idxmask = ( 1 << idxbits ) - 1,
};

enum {
  NUM_NAMESPACES = 1 << nspbits,
  NUM_IDX = 1 << idxbits,
};


void namespaceShowbits ( int n, char *name )
{
  int i;
  unsigned mask;
  char bitvalues[intbits + 1];

  mask = 1;
  for ( i = 0; i < intbits; i++ )
    {
      bitvalues[i] = ((unsigned)n & mask) ? '1':'0';
      mask <<= 1;
    }
  bitvalues[intbits] = '\0';
  fprintf (stdout, "%s: %s\n", name, bitvalues );
}


int namespaceIdxEncode ( namespaceTuple_t tin )
{
  xassert ( tin.nsp < NUM_NAMESPACES && tin.idx < NUM_IDX);
  return ( tin.nsp << idxbits ) + tin.idx;
}

int namespaceIdxEncode2 ( int nsp, int idx )
{
  xassert(nsp < NUM_NAMESPACES && idx < NUM_IDX);
  return ( nsp << idxbits ) + idx;
}


namespaceTuple_t namespaceResHDecode ( int resH )
{
  namespaceTuple_t tin;

  tin.idx = resH & idxmask;
  tin.nsp = (int)(((unsigned)( resH & nspmask )) >> idxbits);

  return tin;
}


void namespaceInit ( int nspn, int * argHasLocalFile )
{
#ifdef USE_MPI
  int nspID;

  xassert(nspn <= NUM_NAMESPACES && nspn >= 1 );

  nNamespaces = nspn;
  if ( nspn >= 1 )
    {
      hasLocalFiles = xmalloc ( nspn * sizeof ( hasLocalFiles[0] ));
      for ( nspID = 0; nspID < nspn; nspID++ )
	hasLocalFiles[nspID] = argHasLocalFile[nspID];
      resStatus = xmalloc ( nspn * sizeof ( resStatus[0] ));
    }
#endif
}


void namespaceCleanup ( void )
{
  if ( nNamespaces > 1 )
    {
      free ( hasLocalFiles );
      hasLocalFiles = NULL;
      free ( resStatus );
    }
}


int namespaceGetNumber ()
{
  return nNamespaces;
}


void pioNamespaceSetActive ( int nId )
{
#ifdef USE_MPI
  if ( nId >= nNamespaces || nId < 0 )
    abort ();

  activeNamespace = nId;
#endif
}


int namespaceGetActive ()
{
  return activeNamespace;
}


int namespaceHasLocalFile ( int nId )
{
  if ( nId >= nNamespaces || nId < 0 )
    abort ();

  return hasLocalFiles ? hasLocalFiles[nId] : 0;
}


int namespaceAdaptKey ( int key, int nspTarget )
{
  namespaceTuple_t tin;
  int nsp;

  if ( key == CDI_UNDEFID ) return CDI_UNDEFID;

  tin.idx = key & idxmask;
  tin.nsp = (int)(((unsigned)( key & nspmask )) >> idxbits);

  xassert ( tin.nsp == nspTarget );

  nsp = namespaceGetActive ();

  return namespaceIdxEncode2 ( nsp, tin.idx );
}


int namespaceAdaptKey2 ( int key )
{
  namespaceTuple_t tin;
  int nsp;

  if ( key == CDI_UNDEFID ) return CDI_UNDEFID;

  tin.idx = key & idxmask;
  tin.nsp = (int)(((unsigned)( key & nspmask )) >> idxbits);

  nsp = namespaceGetActive ();

  return namespaceIdxEncode2 ( nsp, tin.idx );
}


void namespaceDefResStatus ( statusCode argResStatus )
{
  int nsp = namespaceGetActive ();
  resStatus[nsp] = argResStatus;
}


statusCode namespaceInqResStatus ( void )
{
  int nsp = namespaceGetActive ();
  return resStatus[nsp];
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
