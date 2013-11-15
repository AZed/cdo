#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <inttypes.h>

void
memcrc_r(uint32_t *state, const unsigned char *block, size_t block_len);

unsigned long
memcrc_finish(uint32_t *state, off_t total_size);

uint32_t
memcrc(const unsigned char *b, size_t n);

