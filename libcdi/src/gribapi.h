#ifndef _GRIBAPI_H
#define _GRIBAPI_H

#define  GRIBAPI_MISSVAL  -9.E33

/* GRIB2 Level Types */
#define  GRIB2_LTYPE_SURFACE               1
#define  GRIB2_LTYPE_TOA           8
#define  GRIB2_LTYPE_SEA_BOTTOM            9
#define  GRIB2_LTYPE_ATMOSPHERE           10
#define  GRIB2_LTYPE_ISOBARIC            100
#define  GRIB2_LTYPE_MEANSEA             101
#define  GRIB2_LTYPE_ALTITUDE            102
#define  GRIB2_LTYPE_HEIGHT              103
#define  GRIB2_LTYPE_SIGMA               104
#define  GRIB2_LTYPE_HYBRID              105
#define  GRIB2_LTYPE_LANDDEPTH           106
#define  GRIB2_LTYPE_ISENTROPIC          107
#define  GRIB2_LTYPE_REFERENCE           150
#define  GRIB2_LTYPE_SEADEPTH            160

/* GRIB2 Data representation type (Grid Type) */
#define  GRIB2_GTYPE_LATLON                0  /*  latitude/longitude                       */
#define  GRIB2_GTYPE_LATLON_ROT            1  /*  rotated latitude/longitude               */
#define  GRIB2_GTYPE_LATLON_STR            2  /*  stretched latitude/longitude             */
#define  GRIB2_GTYPE_LATLON_ROTSTR         3  /*  rotated and stretched latitude/longitude */
#define  GRIB2_GTYPE_GAUSSIAN             40  /*  gaussian grid                            */
#define  GRIB2_GTYPE_GAUSSIAN_ROT         41  /*  rotated gaussian grid                    */
#define  GRIB2_GTYPE_GAUSSIAN_STR         42  /*  stretched gaussian grid                  */
#define  GRIB2_GTYPE_GAUSSIAN_ROTSTR      43  /*  rotated and stretched gaussian grid      */
#define  GRIB2_GTYPE_LCC                  30  /*  Lambert conformal                        */
#define  GRIB2_GTYPE_SPECTRAL             50  /*  spherical harmonics                      */
#define  GRIB2_GTYPE_GME                 100  /*  hexagonal GME grid                       */
#define  GRIB2_GTYPE_NUMBER              101  /*  General Unstructured Grid                */

const char *gribapiLibraryVersion(void);
void gribContainersNew(int streamID);
void gribContainersDelete(int streamID);
void *gribHandleNew(int editionNumber);
void gribHandleDelete(void *gh);

typedef struct {
  int init;
  void *gribHandle;
}
gribContainer_t;

#endif  /* _GRIBAPI_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
