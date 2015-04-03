#include <inttypes.h>
#include <stdlib.h>

#include "cdi_cksum.h"
#include "cksum.h"
#include "error.h"
#include "serialize.h"

uint32_t cdiCheckSum(int type, int count, void *buffer)
{
  uint32_t s = 0U;
  xassert(count >= 0);
  size_t elemSize = (size_t)serializeGetSizeInCore(1, type, NULL);
  memcrc_r_eswap(&s, (const unsigned char*) buffer, count, elemSize);
  s = memcrc_finish(&s, elemSize * (size_t)count);
  return s;
}
