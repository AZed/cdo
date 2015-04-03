#ifndef _BASETIME_H
#define _BASETIME_H


typedef struct {
  int   ncvarid;
  int   ncdimid;
  int   ncvarboundsid;
  int   lwrf;     /* TRUE for time axis in WRF format */
}
BaseTime;

void basetimeInit(BaseTime *basetime);

#endif  /* _BASETIME_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
