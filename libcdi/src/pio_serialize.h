#ifndef SERIALIZE_PIO_H
#define SERIALIZE_PIO_H
/*
 * Interfaces for marshalling via MPI
 */
int serializeGetSizeMPI(int count, int datatype, void *context);
void serializePackMPI(void *data, int count, int datatype,
                          void *buf, int buf_size, int *position, void *context);
void serializeUnpackMPI(void *buf, int buf_size, int *position,
                        void *data, int count, int datatype, void *context);
/* switch current namespace to use MPI serialization */
void serializeSetMPI();

#endif
