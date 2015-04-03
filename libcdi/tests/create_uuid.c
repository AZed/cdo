#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#define _XOPEN_SOURCE 600

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "create_uuid.h"



#ifdef HAVE_DECL_UUID_GENERATE
#include <sys/time.h>
#include <uuid/uuid.h>
void
create_uuid(unsigned char *uuid)
{
  static int uuid_seeded = 0;
  static char uuid_rand_state[31 * sizeof (long)];
  char *caller_rand_state;
  if (uuid_seeded)
    caller_rand_state = setstate(uuid_rand_state);
  else
    {
      struct timeval tv;
      int status = gettimeofday(&tv, NULL);
      if (status != 0)
        {
          perror("uuid random seed generation failed!");
          exit(1);
        }
      unsigned seed = (unsigned)(tv.tv_sec ^ tv.tv_usec);
      caller_rand_state = initstate(seed, uuid_rand_state,
                                    sizeof (uuid_rand_state));
      uuid_seeded = 1;
    }
  uuid_generate(uuid);
  setstate(caller_rand_state);
}
#elif defined (HAVE_DECL_UUID_CREATE)
typedef uint8_t u_int8_t;
typedef uint16_t u_int16_t;
typedef uint32_t u_int32_t;
#include <uuid.h>
void
create_uuid(unsigned char *uuid)
{
  unsigned32 status;
  uuid_create((uuid_t *)uuid, &status);
  if (status == -1)
    {
      perror("uuid generation failed!");
      exit(1);
    }
}
#else
#include <sys/time.h>
void
create_uuid(unsigned char *uuid)
{
  static int uuid_seeded = 0;
  static char uuid_rand_state[31 * sizeof (long)];
  char *caller_rand_state;
  if (uuid_seeded)
    caller_rand_state = setstate(uuid_rand_state);
  else
    {
      struct timeval tv;
      int status = gettimeofday(&tv, NULL);
      if (status != 0)
        {
          perror("failed seed generation!");
          exit(1);
        }
      unsigned seed = tv.tv_sec ^ tv.tv_usec;
      caller_rand_state = initstate(seed, uuid_rand_state,
                                    sizeof (uuid_rand_state));
      uuid_seeded = 1;
    }
  for (size_t i = 0; i < CDI_UUID_SIZE; ++i)
    uuid[i] = (unsigned char)random();
  /* encode variant into msb of octet 8 */
  uuid[8] = (unsigned char)((uuid[8] & 0x3f) | (1 << 7));
  /* encode version 4 ((pseudo-)random uuid) into msb of octet 7 */
  uuid[7] = (unsigned char)((uuid[7] & 0x0f) | (4 << 4));
  setstate(caller_rand_state);
}
#endif
