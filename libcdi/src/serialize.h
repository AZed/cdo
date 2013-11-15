#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifndef SERIALIZE_H
#define SERIALIZE_H

#include "cdi.h"

/*
 * Generic interfaces for (de-)marshalling
 */
int serializeGetSize(int count, int datatype, void *context);
void serializePack(void *data, int count, int datatype,
                   void *buf, int buf_size, int *position, void *context);
void serializeUnpack(void *buf, int buf_size, int *position,
                     void *data, int count, int datatype, void *context);

/*
 * top-level de-marshalling function
 */

/*
 * Interfaces for marshalling within a single memory domain
 */
int serializeGetSizeInCore(int count, int datatype, void *context);
void serializePackInCore(void *data, int count, int datatype,
                          void *buf, int buf_size, int *position, void *context);
void serializeUnpackInCore(void *buf, int buf_size, int *position,
                            void *data, int count, int datatype, void *context);

#endif
