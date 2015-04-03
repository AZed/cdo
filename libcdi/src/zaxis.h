#ifndef _ZAXIS_H
#define _ZAXIS_H

int zaxisSize(void);

void
zaxisUnpack(char * unpackBuffer, int unpackBufferSize,
            int * unpackBufferPos, int originNamespace, void *context,
            int force_id);

#endif
