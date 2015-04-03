#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "cdi.h"
#include "dmemory.h"
#include "institution.h"
#include "model.h"
#include "cdi_int.h"
#include "vlist.h"
#include "namespace.h"
#include "serialize.h"
#include "resource_unpack.h"
#include "taxis.h"

extern void   gridUnpack ( char *, int, int *, int, void *context);
extern void  zaxisUnpack ( char *, int, int *, int, void *context);
extern int streamNint;


/*****************************************************************************/



void reshUnpackResources(char * unpackBuffer, int unpackBufferSize,
                         void *context)
{
  int token1, token2, nspTarget;
  int unpackBufferPos = 0;
  int numAssociations = 0, sizeAssociations = 16;
  struct streamAssoc *associations
    = xmalloc(sizeof (associations[0]) * sizeAssociations);

  while ( unpackBufferPos < unpackBufferSize )
    {
      serializeUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                      &token1, 1, DATATYPE_INT, context);

      if (token1 == END)
        break;
      switch (token1)
	{
	case START:
	  serializeUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                          &nspTarget, 1, DATATYPE_INT, context);
	  break;
	case GRID:
	  gridUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                     nspTarget, context );
	  break;
	case ZAXIS:
	  zaxisUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                      nspTarget, context);
	  break;
	case TAXIS:
	  taxisUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                      nspTarget, context, 1);
	  break;
	case INSTITUTE:
          instituteUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                          nspTarget, context);
	  break;
	case MODEL:
          modelUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                      nspTarget, context);
	  break;
	case STREAM:
          if (sizeAssociations == numAssociations)
            associations
              = xrealloc(associations,
                         sizeof (associations[0]) * (sizeAssociations *= 2));
	  associations[numAssociations]
            = streamUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                           nspTarget, context);
          ++numAssociations;
	  break;
	case VLIST:
          vlistUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                      nspTarget, context);
	  break;
	default:
	  xabort ( "TOKEN MAPS NO VALID DATATYPE" );
	}

      serializeUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                       &token2, 1, DATATYPE_INT, context);
      xassert ( token2 == SEPARATOR );
    }
  for (int i = 0; i < numAssociations; ++i)
    {
      cdiStreamSetupVlist(stream_to_pointer(associations[i].streamID),
                          namespaceAdaptKey(associations[i].vlistID,
                                            nspTarget),
                          namespaceAdaptKey(associations[i].vlistIDorig,
                                            nspTarget));
    }
  free(associations);
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
