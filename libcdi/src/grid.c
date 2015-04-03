#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <string.h>

#include "dmemory.h"
#include "cdi.h"
#include "stream_int.h"
#include "grid.h"
#include "gaussgrid.h"


#ifndef  RAD2DEG
#define  RAD2DEG  (180./M_PI)   /* conversion for rad to deg */
#endif

#ifndef  DEG2RAD
#define  DEG2RAD  (M_PI/180.)   /* conversion for deg to rad */
#endif


char *Grids[] = {
  /*  0 */  "undefined",
  /*  1 */  "generic",
  /*  2 */  "gaussian",
  /*  3 */  "gaussian reduced",
  /*  4 */  "lonlat",
  /*  5 */  "spectral",
  /*  6 */  "fourier",
  /*  7 */  "gme",
  /*  8 */  "trajectory",
  /*  9 */  "unstructured",
  /* 10 */  "curvilinear",
  /* 11 */  "lcc",
  /* 12 */  "lcc2",
  /* 13 */  "laea",
  /* 14 */  "sinusoidal",
  /* 15 */  "reference",
  /* 16 */  "projection",
};


static int  GRID_Debug = 0;   /* If set to 1, debugging */

static int _grid_max = MAX_GRIDS;

static void grid_initialize(void);

static int _grid_init = FALSE;

#if  defined  (HAVE_LIBPTHREAD)
#  include <pthread.h>

static pthread_once_t _grid_init_thread = PTHREAD_ONCE_INIT;
static pthread_mutex_t _grid_mutex;

#  define GRID_LOCK()         pthread_mutex_lock(&_grid_mutex)
#  define GRID_UNLOCK()       pthread_mutex_unlock(&_grid_mutex)
#  define GRID_INIT()        \
   if ( _grid_init == FALSE ) pthread_once(&_grid_init_thread, grid_initialize)

#else

#  define GRID_LOCK()
#  define GRID_UNLOCK()
#  define GRID_INIT()        \
   if ( _grid_init == FALSE ) grid_initialize()

#endif


typedef struct _gridPtrToIdx {
  int idx;
  grid_t *ptr;
  struct _gridPtrToIdx *next;
} gridPtrToIdx;


static gridPtrToIdx *_gridList  = NULL;
static gridPtrToIdx *_gridAvail = NULL;


static
void grid_list_new(void)
{
  assert(_gridList == NULL);

  _gridList = (gridPtrToIdx *) malloc(_grid_max*sizeof(gridPtrToIdx));
}

static
void grid_list_delete(void)
{
  if ( _gridList ) free(_gridList);
}

static
void grid_init_pointer(void)
{
  int  i;

  for ( i = 0; i < _grid_max; i++ )
    {
      _gridList[i].next = _gridList + i + 1;
      _gridList[i].idx  = i;
      _gridList[i].ptr  = 0;
    }

  _gridList[_grid_max-1].next = 0;

  _gridAvail = _gridList;
}

static
grid_t *grid_to_pointer(int idx)
{
  grid_t *gridptr = NULL;

  GRID_INIT();

  if ( idx >= 0 && idx < _grid_max )
    {
      GRID_LOCK();

      gridptr = _gridList[idx].ptr;

      GRID_UNLOCK();
    }
  else
    Error("grid index %d undefined!", idx);

  return (gridptr);
}

/* Create an index from a pointer */
static
int grid_from_pointer(grid_t *ptr)
{
  int      idx = -1;
  gridPtrToIdx *newptr;

  if ( ptr )
    {
      GRID_LOCK();

      if ( _gridAvail )
	{
	  newptr       = _gridAvail;
	  _gridAvail   = _gridAvail->next;
	  newptr->next = 0;
	  idx	       = newptr->idx;
	  newptr->ptr  = ptr;

	  if ( GRID_Debug )
	    Message("Pointer %p has idx %d from grid list", ptr, idx);
	}
      else
	Warning("Too many open grids (limit is %d)!", _grid_max);

      GRID_UNLOCK();
    }
  else
    Error("Internal problem (pointer %p undefined)", ptr);

  return (idx);
}


void grid_init(grid_t *gridptr)
{
  gridptr->self         = CDI_UNDEFID;
  gridptr->type         = CDI_UNDEFID;
  gridptr->proj         = CDI_UNDEFID;
  gridptr->mask         = NULL;
  gridptr->mask_gme     = NULL;
  gridptr->xvals        = NULL;
  gridptr->yvals        = NULL;
  gridptr->area         = NULL;
  gridptr->xbounds      = NULL;
  gridptr->ybounds      = NULL;
  gridptr->rowlon       = NULL;
  gridptr->nrowlon      = 0;
  gridptr->xinc         = 0.0;
  gridptr->yinc         = 0.0;
  gridptr->lcc_originLon = 0.0;
  gridptr->lcc_originLat = 0.0;
  gridptr->lcc_lonParY  = 0.0;
  gridptr->lcc_lat1     = 0.0;
  gridptr->lcc_lat2     = 0.0;
  gridptr->lcc_xinc     = 0.0;
  gridptr->lcc_yinc     = 0.0;
  gridptr->lcc_projflag = 0;
  gridptr->lcc_scanflag = 0;
  gridptr->lcc_defined  = FALSE;
  gridptr->lcc2_lon_0   = 0.0;
  gridptr->lcc2_lat_0   = 0.0;
  gridptr->lcc2_lat_1   = 0.0;
  gridptr->lcc2_lat_2   = 0.0;
  gridptr->lcc2_a       = 0.0;
  gridptr->lcc2_defined = FALSE;
  gridptr->laea_lon_0   = 0.0;
  gridptr->laea_lat_0   = 0.0;
  gridptr->laea_a       = 0.0;
  gridptr->laea_defined = FALSE;
  gridptr->trunc        = 0;
  gridptr->nvertex      = 0;
  gridptr->nd           = 0;
  gridptr->ni           = 0;
  gridptr->ni2          = 0;
  gridptr->ni3          = 0;
  gridptr->number       = 0;
  gridptr->position     = 0;
  gridptr->reference    = NULL;
  gridptr->prec         = 0;
  gridptr->size         = 0;
  gridptr->xsize        = 0;
  gridptr->ysize        = 0;
  gridptr->np           = 0;
  gridptr->xdef         = 0;
  gridptr->ydef         = 0;
  gridptr->isCyclic     = CDI_UNDEFID;
  gridptr->isRotated    = FALSE;
  gridptr->xpole        = 0.0;
  gridptr->ypole        = 0.0;
  gridptr->angle        = 0.0;
  gridptr->locked       = FALSE;
  gridptr->lcomplex     = 0;
  gridptr->xname[0]     = 0;
  gridptr->yname[0]     = 0;
  gridptr->xlongname[0] = 0;
  gridptr->ylongname[0] = 0;
  gridptr->xunits[0]    = 0;
  gridptr->yunits[0]    = 0;
  gridptr->xstdname[0]  = 0;
  gridptr->ystdname[0]  = 0;
  gridptr->uuid[0]      = 0;
  gridptr->name         = NULL;
}


void grid_free(grid_t *gridptr)
{
  if ( gridptr->mask      ) free(gridptr->mask);
  if ( gridptr->mask_gme  ) free(gridptr->mask_gme);
  if ( gridptr->xvals     ) free(gridptr->xvals);
  if ( gridptr->yvals     ) free(gridptr->yvals);
  if ( gridptr->area      ) free(gridptr->area);
  if ( gridptr->xbounds   ) free(gridptr->xbounds);
  if ( gridptr->ybounds   ) free(gridptr->ybounds);
  if ( gridptr->rowlon    ) free(gridptr->rowlon);
  if ( gridptr->reference ) free(gridptr->reference);
  if ( gridptr->name      ) free(gridptr->name);

  grid_init(gridptr);
}

static
void grid_init_entry(grid_t *gridptr)
{
  grid_init(gridptr);

  gridptr->self = grid_from_pointer(gridptr);
}

static
grid_t *grid_new_entry(void)
{
  grid_t *gridptr;

  gridptr = (grid_t *) malloc(sizeof(grid_t));

  if ( gridptr ) grid_init_entry(gridptr);

  return (gridptr);
}

static
void grid_delete_entry(grid_t *gridptr)
{
  int idx;

  idx = gridptr->self;

  GRID_LOCK();

  free(gridptr);

  _gridList[idx].next = _gridAvail;
  _gridList[idx].ptr  = 0;
  _gridAvail          = &_gridList[idx];

  GRID_UNLOCK();

  if ( GRID_Debug )
    Message("Removed idx %d from grid list", idx);
}

static
void grid_initialize(void)
{
  char *env;

#if  defined  (HAVE_LIBPTHREAD)
  /* initialize global API mutex lock */
  pthread_mutex_init(&_grid_mutex, NULL);
#endif

  env = getenv("GRID_DEBUG");
  if ( env ) GRID_Debug = atoi(env);

  grid_list_new();
  atexit(grid_list_delete);

  GRID_LOCK();

  grid_init_pointer();

  GRID_UNLOCK();

  _grid_init = TRUE;
}

static
void grid_copy(grid_t *gridptr2, grid_t *gridptr1)
{
  int gridID2;

  gridID2 = gridptr2->self;
  memcpy(gridptr2, gridptr1, sizeof(grid_t));
  gridptr2->self = gridID2;
}

static
void gridCheckPtr(const char *caller, int gridID, grid_t *gridptr)
{
  if ( gridptr == NULL )
    Errorc("grid %d undefined!", gridID);
}

#define  grid_check_ptr(gridID, gridptr)  gridCheckPtr(__func__, gridID, gridptr)

int gridSize(void)
{
  int gridsize = 0;
  long i;
  
  GRID_INIT();

  GRID_LOCK();

  for ( i = 0; i < _grid_max; i++ )
    if ( _gridList[i].ptr ) gridsize++;

  GRID_UNLOCK();

  return (gridsize);
}


void gridGenXvals(int xsize, double xfirst, double xlast, double xinc, double *xvals)
{
  long i;

  if ( (! (fabs(xinc) > 0)) && xsize > 1 )
    {
      if ( xfirst >= xlast )
	{
	  while ( xfirst >= xlast ) xlast += 360;
	  xinc = (xlast-xfirst)/(xsize);
	}
      else
	{
	  xinc = (xlast-xfirst)/(xsize-1);
	}
    }

  for ( i = 0; i < xsize; i++ )
    xvals[i] = xfirst + i*xinc;
}

static
void calc_gaussgrid(double *yvals, int ysize, double yfirst, double ylast)
{
  double *yw;
  long yhsize;
  long i;

  yw = (double *) malloc(ysize*sizeof(double));
  gaussaw(yvals, yw, ysize);
  free(yw);
  for ( i = 0; i < ysize; i++ )
    yvals[i] = asin(yvals[i])/M_PI*180.0;

  if ( yfirst < ylast && yfirst > -90.0 && ylast < 90.0 )
    {
      double ytmp;
      yhsize = ysize/2;
      for ( i = 0; i < yhsize; i++ )
	{
	  ytmp = yvals[i];
	  yvals[i] = yvals[ysize-i-1];
	  yvals[ysize-i-1] = ytmp;
	}
    }
}


void gridGenYvals(int gridtype, int ysize, double yfirst, double ylast, double yinc, double *yvals)
{
  long i;
  double deleps = 0.002;

  if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
    {
      if ( ysize > 2 )
	{
	  calc_gaussgrid(yvals, ysize, yfirst, ylast);

	  if ( ! (IS_EQUAL(yfirst, 0) && IS_EQUAL(ylast, 0)) )
	    if ( fabs(yvals[0] - yfirst) > deleps || fabs(yvals[ysize-1] - ylast) > deleps )
	      {
		double yinc = fabs(ylast-yfirst)/(ysize-1);
		double *ytmp = NULL;
		int nstart, lfound = 0;
		int ny = (int) (180./yinc + 0.5);
		ny -= ny%2;
		/* printf("%g %g %g %g %g %d\n", ylast, yfirst, ylast-yfirst,yinc, 180/yinc, ny); */
		if ( ny > ysize && ny < 4096 )
		  {
		    ytmp = (double *) malloc(ny*sizeof(double));
		    calc_gaussgrid(ytmp, ny, yfirst, ylast);
		    for ( i = 0; i < (ny-ysize); i++ )
		      if ( fabs(ytmp[i] - yfirst) < deleps ) break;

		    nstart = i;

		    if ( (nstart+ysize-1) < ny )
		      if ( fabs(ytmp[nstart+ysize-1] - ylast) < deleps ) lfound = 1;
		  }

		if ( lfound )
		  {
		    for ( i = 0; i < ysize; i++ ) yvals[i] = ytmp[i+nstart];
		  }
		else
		  {
		    Warning("Cannot calculate gaussian latitudes for lat1 = %g latn = %g!", yfirst, ylast);
		    for ( i = 0; i < ysize; i++ ) yvals[i] = 0;
		    yvals[0] = yfirst;
		    yvals[ysize-1] = ylast;
		  }

		if ( ytmp ) free(ytmp);
	      }
	}
      else
	{
	  yvals[0] = yfirst;
	  yvals[ysize-1] = ylast;
	}
    }
  /*     else if ( gridtype == GRID_LONLAT || gridtype == GRID_GENERIC ) */
  else
    {
      if ( (! (fabs(yinc) > 0)) && ysize > 1 )
	{
	  if ( IS_EQUAL(yfirst, ylast) && IS_NOT_EQUAL(yfirst, 0) ) ylast *= -1;

	  if ( yfirst > ylast )
	    yinc = (yfirst-ylast)/(ysize-1);
	  else if ( yfirst < ylast )
	    yinc = (ylast-yfirst)/(ysize-1);
	  else
	    {
	      if ( ysize%2 != 0 )
		{
		  yinc = 180.0/(ysize-1);
		  yfirst = -90;
		}
	      else
		{
		  yinc = 180.0/ysize;
		  yfirst = -90 + yinc/2;
		}
	    }
	}

      if ( yfirst > ylast && yinc > 0 ) yinc = -yinc;

      for ( i = 0; i < ysize; i++ )
	yvals[i] = yfirst + i*yinc;
    }
  /*
    else
    Error("unable to calculate values for %s grid!", gridNamePtr(gridtype));
  */
}


/*
@Function  gridCreate
@Title     Create a horizontal Grid

@Prototype int gridCreate(int gridtype, int size)
@Parameter
    @Item  gridtype  The type of the grid, one of the set of predefined CDI grid types.
                     The valid CDI grid types are @func{GRID_GENERIC}, @func{GRID_GAUSSIAN},
                     @func{GRID_LONLAT}, @func{GRID_LCC}, @func{GRID_SPECTRAL},
                     @func{GRID_GME}, @func{GRID_CURVILINEAR}, @func{GRID_UNSTRUCTURED} and
                     @func{GRID_REFERENCE}.
    @Item  size      Number of gridpoints.

@Description
The function @func{gridCreate} creates a horizontal Grid.

@Result
@func{gridCreate} returns an identifier to the Grid.

@Example
Here is an example using @func{gridCreate} to create a regular lon/lat Grid:

@Source
#include "cdi.h"
   ...
#define  nlon  12
#define  nlat   6
   ...
double lons[nlon] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330};
double lats[nlat] = {-75, -45, -15, 15, 45, 75};
int gridID;
   ...
gridID = gridCreate(GRID_LONLAT, nlon*nlat);
gridDefXsize(gridID, nlon);
gridDefYsize(gridID, nlat);
gridDefXvals(gridID, lons);
gridDefYvals(gridID, lats);
   ...
@EndSource
@EndFunction
*/
int gridCreate(int gridtype, int size)
{
  int gridID;
  grid_t *gridptr;

  if ( CDI_Debug )
    Message("gridtype: %d size: %d", gridtype, size);

  GRID_INIT();

  gridptr = grid_new_entry();
  if ( ! gridptr ) Error("No memory");

  gridID = gridptr->self;

  if ( CDI_Debug ) Message("gridID: %d", gridID);

  gridptr->type = gridtype;
  gridptr->size = size;

  /*  if ( gridtype == GRID_GENERIC )     gridptr->xsize = size; */
  if ( gridtype == GRID_UNSTRUCTURED )  gridptr->xsize = size;
  if ( gridtype == GRID_CURVILINEAR  )  gridptr->nvertex = 4;

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_CURVILINEAR:
    case GRID_TRAJECTORY:
      {
	if ( gridtype == GRID_TRAJECTORY )
	  {
	    gridDefXname(gridID, "tlon");
	    gridDefYname(gridID, "tlat");
	  }
	else
	  {
	    gridDefXname(gridID, "lon");
	    gridDefYname(gridID, "lat");
	  }
	gridDefXlongname(gridID, "longitude");
	gridDefYlongname(gridID, "latitude");
	/*
	if ( gridtype == GRID_CURVILINEAR )
	  {
	    strcpy(gridptr->xstdname, "grid_longitude");
	    strcpy(gridptr->ystdname, "grid_latitude");
	    gridDefXunits(gridID, "degrees");
	    gridDefYunits(gridID, "degrees");
	  }
	else
	*/
	  {
	    strcpy(gridptr->xstdname, "longitude");
	    strcpy(gridptr->ystdname, "latitude");
	    gridDefXunits(gridID, "degrees_east");
	    gridDefYunits(gridID, "degrees_north");
	  }

	break;
      }
    case GRID_GME:
    case GRID_UNSTRUCTURED:
      {
	gridDefXname(gridID, "lon");
	gridDefYname(gridID, "lat");
	strcpy(gridptr->xstdname, "longitude");
	strcpy(gridptr->ystdname, "latitude");
	gridDefXunits(gridID, "degrees_east");
	gridDefYunits(gridID, "degrees_north");
	break;
      }
    case GRID_GENERIC:
      {
	gridDefXname(gridID, "x");
	gridDefYname(gridID, "y");
	strcpy(gridptr->xstdname, "grid_longitude");
	strcpy(gridptr->ystdname, "grid_latitude");
	gridDefXunits(gridID, "degrees");
	gridDefYunits(gridID, "degrees");
	break;
      }
    case GRID_LCC2:
    case GRID_SINUSOIDAL:
    case GRID_LAEA:
      {
	gridDefXname(gridID, "x");
	gridDefYname(gridID, "y");
	strcpy(gridptr->xstdname, "projection_x_coordinate");
	strcpy(gridptr->ystdname, "projection_y_coordinate");
	gridDefXunits(gridID, "m");
	gridDefYunits(gridID, "m");
	break;
      }
    }

  return (gridID);
}


/*
@Function  gridDestroy
@Title     Destroy a horizontal Grid

@Prototype void gridDestroy(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.

@EndFunction
*/
void gridDestroy(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);
  
  grid_check_ptr(gridID, gridptr);

  if ( gridptr->mask      ) free(gridptr->mask);
  if ( gridptr->mask_gme  ) free(gridptr->mask_gme);
  if ( gridptr->xvals     ) free(gridptr->xvals);
  if ( gridptr->yvals     ) free(gridptr->yvals);
  if ( gridptr->area      ) free(gridptr->area);
  if ( gridptr->xbounds   ) free(gridptr->xbounds);
  if ( gridptr->ybounds   ) free(gridptr->ybounds);
  if ( gridptr->rowlon    ) free(gridptr->rowlon);
  if ( gridptr->reference ) free(gridptr->reference);

  grid_delete_entry(gridptr);
}


char *gridNamePtr(int gridtype)
{
  char *name;
  int size = (int) (sizeof(Grids)/sizeof(char *));

  if ( gridtype >= 0 && gridtype < size )
    name = Grids[gridtype];
  else
    name = Grids[GRID_GENERIC];

  return (name);
}


void gridName(int gridtype, char *gridname)
{
  strcpy(gridname, gridNamePtr(gridtype));
}


/*
@Function  gridDefXname
@Title     Define the name of a X-axis

@Prototype void gridDefXname(int gridID, const char *name)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  name     Name of the X-axis.

@Description
The function @func{gridDefXname} defines the name of a X-axis.

@EndFunction
*/
void gridDefXname(int gridID, const char *xname)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( xname )
    strcpy(gridptr->xname, xname);
}


/*
@Function  gridDefXlongname
@Title     Define the longname of a X-axis

@Prototype void gridDefXlongname(int gridID, const char *longname)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  longname Longname of the X-axis.

@Description
The function @func{gridDefXlongname} defines the longname of a X-axis.

@EndFunction
*/
void gridDefXlongname(int gridID, const char *xlongname)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  if ( xlongname )
    strcpy(gridptr->xlongname, xlongname);
}


/*
@Function  gridDefXunits
@Title     Define the units of a X-axis

@Prototype void gridDefXunits(int gridID, const char *units)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  units    Units of the X-axis.

@Description
The function @func{gridDefXunits} defines the units of a X-axis.

@EndFunction
*/
void gridDefXunits(int gridID, const char *xunits)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( xunits )
    strcpy(gridptr->xunits, xunits);
}


/*
@Function  gridDefYname
@Title     Define the name of a Y-axis

@Prototype void gridDefYname(int gridID, const char *name)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  name     Name of the Y-axis.

@Description
The function @func{gridDefYname} defines the name of a Y-axis.

@EndFunction
*/
void gridDefYname(int gridID, const char *yname)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( yname )
    strcpy(gridptr->yname, yname);
}


/*
@Function  gridDefYlongname
@Title     Define the longname of a Y-axis

@Prototype void gridDefYlongname(int gridID, const char *longname)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  longname Longname of the Y-axis.

@Description
The function @func{gridDefYlongname} defines the longname of a Y-axis.

@EndFunction
*/
void gridDefYlongname(int gridID, const char *ylongname)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( ylongname )
    strcpy(gridptr->ylongname, ylongname);
}


/*
@Function  gridDefYunits
@Title     Define the units of a Y-axis

@Prototype void gridDefYunits(int gridID, const char *units)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  units    Units of the Y-axis.

@Description
The function @func{gridDefYunits} defines the units of a Y-axis.

@EndFunction
*/
void gridDefYunits(int gridID, const char *yunits)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( yunits )
    strcpy(gridptr->yunits, yunits);
}


/*
@Function  gridInqXname
@Title     Get the name of a X-axis

@Prototype void gridInqXname(int gridID, char *name)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  name     Name of the X-axis. The caller must allocate space for the 
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqXname} returns the name of a X-axis.

@Result
@func{gridInqXname} returns the name of the X-axis to the parameter name.

@EndFunction
*/
void gridInqXname(int gridID, char *xname)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  strcpy(xname, gridptr->xname);
}


/*
@Function  gridInqXlongname
@Title     Get the longname of a X-axis

@Prototype void gridInqXlongname(int gridID, char *longname)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  longname Longname of the X-axis. The caller must allocate space for the 
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqXlongname} returns the longname of a X-axis.

@Result
@func{gridInqXlongname} returns the longname of the X-axis to the parameter longname.

@EndFunction
*/
void gridInqXlongname(int gridID, char *xlongname)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  strcpy(xlongname, gridptr->xlongname);
}


/*
@Function  gridInqXunits
@Title     Get the units of a X-axis

@Prototype void gridInqXunits(int gridID, char *units)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  units    Units of the X-axis. The caller must allocate space for the 
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqXunits} returns the units of a X-axis.

@Result
@func{gridInqXunits} returns the units of the X-axis to the parameter units.

@EndFunction
*/
void gridInqXunits(int gridID, char *xunits)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  strcpy(xunits, gridptr->xunits);
}


void gridInqXstdname(int gridID, char *xstdname)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  strcpy(xstdname, gridptr->xstdname);
}


/*
@Function  gridInqYname
@Title     Get the name of a Y-axis

@Prototype void gridInqYname(int gridID, char *name)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  name     Name of the Y-axis. The caller must allocate space for the 
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqYname} returns the name of a Y-axis.

@Result
@func{gridInqYname} returns the name of the Y-axis to the parameter name.

@EndFunction
*/
void gridInqYname(int gridID, char *yname)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  strcpy(yname, gridptr->yname);
}


/*
@Function  gridInqYlongname
@Title     Get the longname of a Y-axis

@Prototype void gridInqXlongname(int gridID, char *longname)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  longname Longname of the Y-axis. The caller must allocate space for the 
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqYlongname} returns the longname of a Y-axis.

@Result
@func{gridInqYlongname} returns the longname of the Y-axis to the parameter longname.

@EndFunction
*/
void gridInqYlongname(int gridID, char *ylongname)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  strcpy(ylongname, gridptr->ylongname);
}


/*
@Function  gridInqYunits
@Title     Get the units of a Y-axis

@Prototype void gridInqYunits(int gridID, char *units)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  units    Units of the Y-axis. The caller must allocate space for the 
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqYunits} returns the units of a Y-axis.

@Result
@func{gridInqYunits} returns the units of the Y-axis to the parameter units.

@EndFunction
*/
void gridInqYunits(int gridID, char *yunits)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  strcpy(yunits, gridptr->yunits);
}

void gridInqYstdname(int gridID, char *ystdname)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  strcpy(ystdname, gridptr->ystdname);
}


/*
@Function  gridInqType
@Title     Get the type of a Grid

@Prototype int gridInqType(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.

@Description
The function @func{gridInqType} returns the type of a Grid.

@Result
@func{gridInqType} returns the type of the grid,
one of the set of predefined CDI grid types.
The valid CDI grid types are @func{GRID_GENERIC}, @func{GRID_GAUSSIAN},
@func{GRID_LONLAT}, @func{GRID_LCC}, @func{GRID_SPECTRAL}, @func{GRID_GME},
@func{GRID_CURVILINEAR}, @func{GRID_UNSTRUCTURED} and @func{GRID_REFERENCE}.

@EndFunction
*/
int gridInqType(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->type);
}


/*
@Function  gridInqSize
@Title     Get the size of a Grid

@Prototype int gridInqSize(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.

@Description
The function @func{gridInqSize} returns the size of a Grid.

@Result
@func{gridInqSize} returns the number of grid points of a Grid.

@EndFunction
*/
int gridInqSize(int gridID)
{
  int size = 0;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  size = gridptr->size;

  if ( ! size )
    {
      int xsize, ysize;

      xsize = gridptr->xsize;
      ysize = gridptr->ysize;

      if ( ysize )
	size = xsize *ysize;
      else
	size = xsize;

      gridptr->size = size;  
    }

  return (size);
}

static
int nsp2trunc(int nsp)
{
  /*  nsp = (trunc+1)*(trunc+1)              */
  /*      => trunc^2 + 3*trunc - (x-2) = 0   */
  /*                                         */
  /*  with:  y^2 + p*y + q = 0               */
  /*         y = -p/2 +- sqrt((p/2)^2 - q)   */
  /*         p = 3 and q = - (x-2)           */
  int trunc;

  trunc = (int) (sqrt(nsp*4 + 1.) - 3) / 2;

  return (trunc);
}


int gridInqTrunc(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->trunc == 0 )
    {
      if ( gridptr->type == GRID_SPECTRAL )
	gridptr->trunc = nsp2trunc(gridptr->size);
      /*
      else if      ( gridptr->type == GRID_GAUSSIAN )
	gridptr->trunc = nlat2trunc(gridptr->ysize);
      */
    }

  return (gridptr->trunc);
}


void gridDefTrunc(int gridID, int trunc)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridptr->trunc = trunc;
}


/*
@Function  gridDefXsize
@Title     Define the number of values of a X-axis

@Prototype void gridDefXsize(int gridID, int xsize)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  xsize    Number of values of a X-axis.

@Description
The function @func{gridDefXsize} defines the number of values of a X-axis.

@EndFunction
*/
void gridDefXsize(int gridID, int xsize)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( xsize > gridInqSize(gridID) )
    Error("xsize %d is greater then gridsize %d", xsize, gridInqSize(gridID));

  if ( gridInqType(gridID) == GRID_UNSTRUCTURED && xsize != gridInqSize(gridID) )
    Error("xsize %d must be equal gridsize %d for gridtype: UNSTRUCTURED", xsize, gridInqSize(gridID));

  gridptr->xsize = xsize;

  if ( gridInqType(gridID) != GRID_UNSTRUCTURED )
    {
      long gridsize = gridptr->xsize*gridptr->ysize;
      if ( gridsize > 0 && gridsize != gridInqSize(gridID) )
        Error("Inconsistent grid declaration! (xsize=%d ysize=%d gridsize=%d)",
              gridptr->xsize, gridptr->ysize, gridInqSize(gridID));
    }
}


/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridDefPrec(int gridID, int prec)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridptr->prec = prec;
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
int gridInqPrec(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->prec);
}


/*
@Function  gridInqXsize
@Title     Get the number of values of a X-axis

@Prototype int gridInqXsize(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.

@Description
The function @func{gridInqXsize} returns the number of values of a X-axis.

@Result
@func{gridInqXsize} returns the number of values of a X-axis.

@EndFunction
*/
int gridInqXsize(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->xsize);
}

/*
@Function  gridDefYsize
@Title     Define the number of values of a Y-axis

@Prototype void gridDefYsize(int gridID, int ysize)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  ysize    Number of values of a Y-axis.

@Description
The function @func{gridDefYsize} defines the number of values of a Y-axis.

@EndFunction
*/
void gridDefYsize(int gridID, int ysize)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( ysize > gridInqSize(gridID) )
    Error("ysize %d is greater then gridsize %d", ysize, gridInqSize(gridID));

  if ( gridInqType(gridID) == GRID_UNSTRUCTURED && ysize != gridInqSize(gridID) )
    Error("ysize %d must be equal gridsize %d for gridtype: UNSTRUCTURED", ysize, gridInqSize(gridID));

  gridptr->ysize = ysize;

  if ( gridInqType(gridID) != GRID_UNSTRUCTURED )
    {
      long gridsize = gridptr->xsize*gridptr->ysize;
      if ( gridsize > 0 && gridsize != gridInqSize(gridID) )
        Error("Inconsistent grid declaration! (xsize=%d ysize=%d gridsize=%d)",
              gridptr->xsize, gridptr->ysize, gridInqSize(gridID));
    }
}

/*
@Function  gridInqYsize
@Title     Get the number of values of a Y-axis

@Prototype int gridInqYsize(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.

@Description
The function @func{gridInqYsize} returns the number of values of a Y-axis.

@Result
@func{gridInqYsize} returns the number of values of a Y-axis.

@EndFunction
*/
int gridInqYsize(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->ysize);
}

/*
@Function  gridDefNP
@Title     Define the number of parallels between a pole and the equator

@Prototype void gridDefNP(int gridID, int np)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  np       Number of parallels between a pole and the equator.

@Description
The function @func{gridDefNP} defines the number of parallels between a pole and the equator
of a Gaussian grid.

@EndFunction
*/
void gridDefNP(int gridID, int np)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridptr->np = np;
}

/*
@Function  gridInqNP
@Title     Get the number of parallels between a pole and the equator

@Prototype int gridInqNP(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.

@Description
The function @func{gridInqNP} returns the number of parallels between a pole and the equator
of a Gaussian grid.

@Result
@func{gridInqNP} returns the number of parallels between a pole and the equator.

@EndFunction
*/
int gridInqNP(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->np);
}

/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridDefRowlon(int gridID, int nrowlon, const int *rowlon)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);
  
  grid_check_ptr(gridID, gridptr);

  gridptr->rowlon = (int *) malloc(nrowlon*sizeof(int));
  gridptr->nrowlon = nrowlon;

  memcpy(gridptr->rowlon, rowlon, nrowlon*sizeof(int));
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridInqRowlon(int gridID, int *rowlon)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->rowlon == 0 )  Error("undefined pointer!");

  memcpy(rowlon, gridptr->rowlon, gridptr->nrowlon*sizeof(int));
}


int gridInqMask(int gridID, int *mask)
{
  long size, i;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  size = gridptr->size;

  if ( CDI_Debug && size == 0 )
    Warning("Size undefined for gridID = %d", gridID);
    
  if ( mask && gridptr->mask )
    for ( i = 0; i < size; ++i )
      mask[i] = gridptr->mask[i];

  if ( gridptr->mask == NULL ) size = 0;

  return (size);
}


void gridDefMask(int gridID, const int *mask)
{
  long size, i;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  size = gridptr->size;

  if ( size == 0 )
    Error("Size undefined for gridID = %d", gridID);

  if ( mask == NULL )
    {
      if ( gridptr->mask )
	{
	  free(gridptr->mask);
	  gridptr->mask = NULL;
	}
    }
  else
    {
      if ( gridptr->mask == NULL )
	gridptr->mask = (mask_t *) malloc(size*sizeof(mask_t));
      else if ( CDI_Debug )
	Warning("grid mask already defined!");

      for ( i = 0; i < size; ++i )
	gridptr->mask[i] = mask[i];
    }
}


int gridInqMaskGME(int gridID, int *mask)
{
  long size, i;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  size = gridptr->size;

  if ( CDI_Debug && size == 0 )
    Warning("Size undefined for gridID = %d", gridID);
    
  if ( mask && gridptr->mask_gme )
    for ( i = 0; i < size; ++i )
      mask[i] = gridptr->mask_gme[i];

  if ( gridptr->mask_gme == NULL ) size = 0;

  return (size);
}


void gridDefMaskGME(int gridID, const int *mask)
{
  long size, i;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  size = gridptr->size;

  if ( size == 0 )
    Error("Size undefined for gridID = %d", gridID);
    
  if ( gridptr->mask_gme == NULL )
    gridptr->mask_gme = (mask_t *) malloc(size*sizeof(mask_t));
  else if ( CDI_Debug )
    Warning("mask already defined!");

  for ( i = 0; i < size; ++i )
    gridptr->mask_gme[i] = mask[i];
}


/*
@Function  gridInqXvals
@Title     Get all values of a X-axis

@Prototype int gridInqXvals(int gridID, double *xvals)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  xvals    Pointer to the location into which the X-values are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{gridInqXvals} returns all values of the X-axis.

@Result
Upon successful completion @func{gridInqXvals} returns the number of values and
the values are stored in @func{xvals}.
Otherwise, 0 is returned and @func{xvals} is empty.

@EndFunction
*/
int gridInqXvals(int gridID, double *xvals)
{
  long size;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->type == GRID_CURVILINEAR || gridptr->type == GRID_UNSTRUCTURED )
    size = gridptr->size;
  else
    size = gridptr->xsize;

  if ( CDI_Debug && size == 0 )
    Warning("Size undefined for gridID = %d", gridID);
    
  if ( xvals && gridptr->xvals )
    memcpy(xvals, gridptr->xvals, size*sizeof(double));

  if ( gridptr->xvals == NULL ) size = 0;

  return (size);
}


/*
@Function  gridDefXvals
@Title     Define the values of a X-axis

@Prototype void gridDefXvals(int gridID, const double *xvals)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  xvals    X-values of the grid.

@Description
The function @func{gridDefXvals} defines all values of the X-axis.

@EndFunction
*/
void gridDefXvals(int gridID, const double *xvals)
{
  int gridtype;
  long size;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridtype = gridptr->type;

  if ( gridtype == GRID_UNSTRUCTURED || gridtype == GRID_CURVILINEAR )
    size = gridptr->size;
  else
    size = gridptr->xsize;

  if ( size == 0 )
    Error("Size undefined for gridID = %d", gridID);
    
  if ( gridptr->xvals == NULL )
    gridptr->xvals = (double *) malloc(size*sizeof(double));
  else if ( CDI_Debug )
    Warning("values already defined!");

  memcpy(gridptr->xvals, xvals, size*sizeof(double));
}


/*
@Function  gridInqYvals
@Title     Get all values of a Y-axis

@Prototype int gridInqYvals(int gridID, double *yvals)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  yvals    Pointer to the location into which the Y-values are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{gridInqYvals} returns all values of the Y-axis.

@Result
Upon successful completion @func{gridInqYvals} returns the number of values and
the values are stored in @func{yvals}.
Otherwise, 0 is returned and @func{yvals} is empty.

@EndFunction
*/
int gridInqYvals(int gridID, double *yvals)
{
  long size;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->type == GRID_CURVILINEAR || gridptr->type == GRID_UNSTRUCTURED )
    size = gridptr->size;
  else
    size = gridptr->ysize;

  if ( CDI_Debug && size == 0 )
    Warning("Size undefined for gridID = %d!", gridID);
    
  if ( yvals && gridptr->yvals )
    memcpy(yvals, gridptr->yvals, size*sizeof(double));

  if ( gridptr->yvals == NULL ) size = 0;

  return (size);
}


/*
@Function  gridDefYvals
@Title     Define the values of a Y-axis

@Prototype void gridDefYvals(int gridID, const double *yvals)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  yvals    Y-values of the grid.

@Description
The function @func{gridDefYvals} defines all values of the Y-axis.

@EndFunction
*/
void gridDefYvals(int gridID, const double *yvals)
{
  int gridtype;
  long size;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridtype = gridptr->type;

  if ( gridtype == GRID_UNSTRUCTURED || gridtype == GRID_CURVILINEAR )
    size = gridptr->size;
  else
    size = gridptr->ysize;

  if ( size == 0 )
    Error("Size undefined for gridID = %d!", gridID);
    
  if ( gridptr->yvals == NULL )
    gridptr->yvals = (double *) malloc(size*sizeof(double));
  else if ( CDI_Debug )
    Warning("Values already defined!");
    
  memcpy(gridptr->yvals, yvals, size*sizeof(double));
}


double gridInqXval(int gridID, int index)
{
  double xval = 0;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->xvals )
    xval = gridptr->xvals[index];

  return (xval);
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
double gridInqYval(int gridID, int index)
{
  double yval = 0;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->yvals )
    yval = gridptr->yvals[index];

  return (yval);
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
double gridInqXinc(int gridID)
{
  double xinc;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  xinc = gridptr->xinc;

  if ( (! (fabs(xinc) > 0)) && gridptr->xvals )
    {
      int xsize;
      double *xvals;

      xsize = gridptr->xsize;
      xvals = gridptr->xvals;

      if ( xsize > 1 )
	{
	  long i;
	  xinc = fabs(xvals[xsize-1] - xvals[0])/(xsize-1);
	  for ( i = 2; i < xsize; i++ )
	    if ( fabs(fabs(xvals[i-1] - xvals[i]) - xinc) > 0.01*xinc ) break;
		  
	  if ( i < xsize ) xinc = 0;
	}
    }

  return (xinc);
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
double gridInqYinc(int gridID)
{
  double yinc;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  yinc = gridptr->yinc;

  if ( (! (fabs(yinc) > 0)) && gridptr->yvals )
    {
      int ysize;
      double *yvals;

      ysize = gridptr->ysize;
      yvals = gridptr->yvals;

      if ( ysize > 1 )
	{
	  long i;
	  yinc = fabs(yvals[1] - yvals[0]);
	  for ( i = 2; i < ysize; i++ )
	    if ( fabs(fabs(yvals[i] - yvals[i-1]) - yinc) > (yinc/1000) ) break;

	  if ( i < ysize ) yinc = 0;
	  else             yinc = yvals[1] - yvals[0];
	}
    }

  return (yinc);
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
double gridInqXpole(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->xpole);
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridDefXpole(int gridID, double xpole)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( memcmp(gridptr->xstdname, "grid", 4) != 0 )
    strcpy(gridptr->xstdname, "grid_longitude");

  gridptr->isRotated = TRUE;
  gridptr->xpole = xpole;
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
double gridInqYpole(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->ypole);
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridDefYpole(int gridID, double ypole)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( memcmp(gridptr->ystdname, "grid", 4) != 0 )
    strcpy(gridptr->ystdname, "grid_latitude");

  gridptr->isRotated = TRUE;
  gridptr->ypole = ypole;
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
double gridInqAngle(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->angle);
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridDefAngle(int gridID, double angle)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridptr->isRotated = TRUE;
  gridptr->angle = angle;
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
int gridInqGMEnd(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->nd);
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridDefGMEnd(int gridID, int nd)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridptr->nd = nd;
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
int gridInqGMEni(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->ni);
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridDefGMEni(int gridID, int ni)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridptr->ni = ni;
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
int gridInqGMEni2(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->ni2);
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridDefGMEni2(int gridID, int ni2)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridptr->ni2 = ni2;
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
int gridInqGMEni3(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->ni3);
}

void gridDefGMEni3(int gridID, int ni3)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridptr->ni3 = ni3;
}


/*
@Function  
@Title     

@Prototype 
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridChangeType(int gridID, int gridtype)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  Message("Changed grid type from %s to %s",
	  gridNamePtr(gridptr->type),
	  gridNamePtr(gridtype));
  
  gridptr->type = gridtype;
}


static 
void grid_check_cyclic(grid_t *gridptr)
{
  int xsize, ysize, gridsize;
  long i1, i2, in, j, k1, k2, nc;
  double xinc, x0;
  const double *xvals, *xbounds;

  gridptr->isCyclic = FALSE;

  gridsize = gridptr->size;
  xsize = gridptr->xsize;
  ysize = gridptr->ysize;
  xvals = gridptr->xvals;
  xbounds = gridptr->xbounds;

  if ( gridptr->type == GRID_GAUSSIAN || gridptr->type == GRID_LONLAT )
    {
      if ( xvals && xsize > 1 )
	{
	  xinc = xvals[1] - xvals[0];
	  if ( IS_EQUAL(xinc, 0) )
	    xinc = (xvals[xsize-1] - xvals[0])/(xsize-1);
	  x0 = 2*xvals[xsize-1]-xvals[xsize-2]-360;
	  if ( IS_NOT_EQUAL(xvals[0], xvals[xsize-1]) )
	    if ( fabs(x0 - xvals[0]) < 0.01*xinc ) gridptr->isCyclic = TRUE;
	}
    }
  else if ( gridptr->type == GRID_CURVILINEAR )
    {
      if ( xvals && xsize > 1 )
	{
	  double val1, val2, valn;

	  nc = 0;
	  gridptr->isCyclic = FALSE;
	  for ( j = 0; j < ysize; ++j )
	    {
	      i1 = j*xsize;
	      i2 = j*xsize+1;
	      in = j*xsize+(xsize-1);
	      val1 = xvals[i1];
	      val2 = xvals[i2];
	      valn = xvals[in];

	      xinc = fabs(val2-val1);

	      if ( val1 <    1 && valn > 300 ) val1 += 360;
	      if ( valn <    1 && val1 > 300 ) valn += 360;
	      if ( val1 < -179 && valn > 120 ) val1 += 360;
	      if ( valn < -179 && val1 > 120 ) valn += 360;

	      if ( valn > val1 ) x0 = valn - xinc;
	      else               x0 = valn + xinc;

	      if ( fabs(x0-val1) < 0.5*xinc ) nc++;
	    }

	  if ( nc > 0.5*ysize ) gridptr->isCyclic = TRUE;
	}

      if ( xbounds && xsize > 1 )
	{
	  double val1, val2;

	  gridptr->isCyclic = TRUE;
	  for ( j = 0; j < ysize; ++j )
	    {
	      i1 = j*xsize*4;
	      i2 = j*xsize*4+(xsize-1)*4;
	      nc = 0;
	      for ( k1 = 0; k1 < 4; ++k1 )
		{
		  val1 = xbounds[i1+k1];
		  for ( k2 = 0; k2 < 4; ++k2 )
		    {
		      val2 = xbounds[i2+k2];

		      if ( val1 <    1 && val2 > 300 ) val1 += 360;
		      if ( val2 <    1 && val1 > 300 ) val2 += 360;
		      if ( val1 < -179 && val2 > 120 ) val1 += 360;
		      if ( val2 < -179 && val1 > 120 ) val2 += 360;

		      if ( fabs(val1-val2) < 0.001 )
			{
			  nc++;
			  break;
			}
		    }
		}

	      if ( nc < 1 )
		{
		  gridptr->isCyclic = FALSE;
		  break;
		}
	    }
	}
    }
}


int gridIsCircular(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->isCyclic == CDI_UNDEFID ) grid_check_cyclic(gridptr);

  return ( gridptr->isCyclic );
}


int gridIsRotated(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return ( gridptr->isRotated );
}

static
int compareXYvals(int gridID, long xsize, long ysize, double *xvals0, double *yvals0)
{
  long i;
  int differ = 0;

  if ( !differ && xsize == gridInqXvals(gridID, NULL) )
    {
      double *xvals;

      xvals = (double *) malloc(xsize*sizeof(double));

      gridInqXvals(gridID, xvals);

      for ( i = 0; i < xsize; ++i )
	if ( fabs(xvals0[i] - xvals[i]) > 1.e-10 )
	  {
	    differ = 1;
	    break;
	  }
  
      free(xvals);
    }

  if ( !differ && ysize == gridInqYvals(gridID, NULL) )
    {
      double *yvals;

      yvals = (double *) malloc(ysize*sizeof(double));

      gridInqYvals(gridID, yvals);

      for ( i = 0; i < ysize; ++i )
	if ( fabs(yvals0[i] - yvals[i]) > 1.e-10 )
	  {
	    differ = 1;
	    break;
	  }
      
      free(yvals);
    }

  return (differ);
}

static
int compareXYvals2(int gridID, long gridsize, double *xvals, double *yvals)
{
  int differ = 0;

  if ( !differ && xvals && gridInqXvalsPtr(gridID) )
    {
      if ( fabs(xvals[0] - gridInqXval(gridID, 0)) > 1.e-9 ||
	   fabs(xvals[gridsize-1] - gridInqXval(gridID, gridsize-1)) > 1.e-9 )
	differ = 1;
    }

  if ( !differ && yvals && gridInqYvalsPtr(gridID) )
    {
      if ( fabs(yvals[0] - gridInqYval(gridID, 0)) > 1.e-9 ||
	   fabs(yvals[gridsize-1] - gridInqYval(gridID, gridsize-1)) > 1.e-9 )
	differ = 1;
    }

  return (differ);
}


int gridCompare(int gridID, grid_t grid)
{
  int differ = 1;

  if ( grid.type == gridInqType(gridID) || grid.type == GRID_GENERIC )
    {
      if ( grid.size == gridInqSize(gridID) )
	{
	  differ = 0;
	  if ( grid.type == GRID_LONLAT )
	    {
	      /*
	      printf("gridID      %d\n", gridID);
	      printf("grid.xdef   %d\n", grid.xdef);
	      printf("grid.ydef   %d\n", grid.ydef);
	      printf("grid.xsize  %d\n", grid.xsize);
	      printf("grid.ysize  %d\n", grid.ysize);
	      printf("grid.xfirst %f\n", grid.xfirst);
	      printf("grid.yfirst %f\n", grid.yfirst);
	      printf("grid.xfirst %f\n", gridInqXval(gridID, 0));
	      printf("grid.yfirst %f\n", gridInqYval(gridID, 0));
	      printf("grid.xinc   %f\n", grid.xinc);
	      printf("grid.yinc   %f\n", grid.yinc);
	      printf("grid.xinc   %f\n", gridInqXinc(gridID));
	      printf("grid.yinc   %f\n", gridInqYinc(gridID));
	      */
	      if ( grid.xsize == gridInqXsize(gridID) && grid.ysize == gridInqYsize(gridID) )
		{
		  if ( grid.xdef == 2 && grid.ydef == 2 )
		    {
		      if ( ! (IS_EQUAL(grid.xfirst, 0) && IS_EQUAL(grid.xlast, 0) && IS_EQUAL(grid.xinc, 0)) &&
			   ! (IS_EQUAL(grid.yfirst, 0) && IS_EQUAL(grid.ylast, 0) && IS_EQUAL(grid.yinc, 0)) &&
			   IS_NOT_EQUAL(grid.xfirst, grid.xlast) && IS_NOT_EQUAL(grid.yfirst, grid.ylast) )
			{
			  if ( IS_NOT_EQUAL(grid.xfirst, gridInqXval(gridID, 0)) ||
			       IS_NOT_EQUAL(grid.yfirst, gridInqYval(gridID, 0)))
			    {
			      differ = 1;
			    }
			  if ( !differ && fabs(grid.xinc) > 0 &&
			       fabs(fabs(grid.xinc) - fabs(gridInqXinc(gridID))) > fabs(grid.xinc/1000))
			    {
			      differ = 1;
			    }
			  if ( !differ && fabs(grid.yinc) > 0 &&
			       fabs(fabs(grid.yinc) - fabs(gridInqYinc(gridID))) > fabs(grid.yinc/1000))
			    {
			      differ = 1;
			    }
			}
		    }
		  else
		    {
		      if ( grid.xvals && grid.yvals )
			differ = compareXYvals(gridID, grid.xsize, grid.ysize, grid.xvals, grid.yvals);
		    }
		}
	      else
		differ = 1;
	    }
	  else if ( grid.type == GRID_GENERIC )
	    {
	      if ( grid.xsize == gridInqXsize(gridID) && grid.ysize == gridInqYsize(gridID) )
		{
		  if ( grid.xdef == 1 && grid.ydef == 1 )
		    {
		      if ( grid.xvals && grid.yvals )
			differ = compareXYvals(gridID, grid.xsize, grid.ysize, grid.xvals, grid.yvals);
		    }
		}
	      else if ( (grid.ysize == 0 || grid.ysize == 1) &&
			grid.xsize == gridInqXsize(gridID)*gridInqYsize(gridID) )
		{
		}
	      else
		differ = 1;		
	    }
	  else if ( grid.type == GRID_GAUSSIAN )
	    {
	      if ( grid.xsize == gridInqXsize(gridID) && grid.ysize == gridInqYsize(gridID) )
		{
		  if ( grid.xdef == 2 && grid.ydef == 2 )
		    {
		      if ( ! (IS_EQUAL(grid.xfirst, 0) && IS_EQUAL(grid.xlast, 0) && IS_EQUAL(grid.xinc, 0)) &&
			   ! (IS_EQUAL(grid.yfirst, 0) && IS_EQUAL(grid.ylast, 0)) )
			if ( fabs(grid.xfirst - gridInqXval(gridID, 0)) > 0.001 ||
			     fabs(grid.yfirst - gridInqYval(gridID, 0)) > 0.001 ||
			     (fabs(grid.xinc)>0 && fabs(fabs(grid.xinc) - fabs(gridInqXinc(gridID))) > fabs(grid.xinc/1000)) )
			  {
			    differ = 1;
			  }
		    }
		  else
		    {
		      if ( grid.xvals && grid.yvals )
			differ = compareXYvals(gridID, grid.xsize, grid.ysize, grid.xvals, grid.yvals);
		    }
		}
	      else
		differ = 1;		
	    }
	  else if ( grid.type == GRID_CURVILINEAR )
	    {
	      /*
	      printf("gridID      %d\n", gridID);
	      printf("grid.xsize  %d\n", grid.xsize);
	      printf("grid.ysize  %d\n", grid.ysize);
	      printf("grid.xfirst %f\n", grid.xvals[0]);
	      printf("grid.yfirst %f\n", grid.yvals[0]);
	      printf("grid xfirst %f\n", gridInqXval(gridID, 0));
	      printf("grid yfirst %f\n", gridInqYval(gridID, 0));
	      printf("grid.xlast  %f\n", grid.xvals[grid.size-1]);
	      printf("grid.ylast  %f\n", grid.yvals[grid.size-1]);
	      printf("grid xlast  %f\n", gridInqXval(gridID, grid.size-1));
	      printf("grid ylast  %f\n", gridInqYval(gridID, grid.size-1));
	      printf("grid.nv     %d\n", grid.nvertex);
	      printf("grid nv     %d\n", gridInqNvertex(gridID));
	      */
	      if ( grid.xsize == gridInqXsize(gridID) && grid.ysize == gridInqYsize(gridID) )
		differ = compareXYvals2(gridID, grid.size, grid.xvals, grid.yvals);
	    }
	  else if ( grid.type == GRID_UNSTRUCTURED )
	    {
	      if ( grid.nvertex == gridInqNvertex(gridID) )
		differ = compareXYvals2(gridID, grid.size, grid.xvals, grid.yvals);
	    }
	}
    }

  return (differ);
}


int gridGenerate(grid_t grid)
{
  int gridID;
  grid_t *gridptr;

  gridID = gridCreate(grid.type, grid.size);

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridDefPrec(gridID, grid.prec);

  switch (grid.type)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_UNSTRUCTURED:
    case GRID_CURVILINEAR:
    case GRID_GENERIC:
    case GRID_LCC:
    case GRID_LCC2:
    case GRID_SINUSOIDAL:
    case GRID_LAEA:
    case GRID_PROJECTION:
      {
	if ( grid.xsize > 0 ) gridDefXsize(gridID, grid.xsize);
	if ( grid.ysize > 0 ) gridDefYsize(gridID, grid.ysize);

        if ( grid.type == GRID_GAUSSIAN ) gridDefNP(gridID, grid.np);

	if ( grid.nvertex > 0 )
	  gridDefNvertex(gridID, grid.nvertex);

	if ( grid.xdef == 1 )
	  {
	    gridDefXvals(gridID, grid.xvals);
	    if ( grid.xbounds )
	      gridDefXbounds(gridID, grid.xbounds);
	  }
	else if ( grid.xdef == 2 )
	  {
	    double *xvals = (double *) malloc(grid.xsize*sizeof(double));
	    gridGenXvals(grid.xsize, grid.xfirst, grid.xlast, grid.xinc, xvals);
	    gridDefXvals(gridID, xvals);
	    free(xvals);
	    /*
	    gridDefXinc(gridID, grid.xinc);
	    */
	  }

	if ( grid.ydef == 1 )
	  {
	    gridDefYvals(gridID, grid.yvals);
	    if ( grid.ybounds && grid.nvertex )
	      gridDefYbounds(gridID, grid.ybounds);
	  }
	else if ( grid.ydef == 2 )
	  {
	    double *yvals = (double *) malloc(grid.ysize*sizeof(double));
	    gridGenYvals(grid.type, grid.ysize, grid.yfirst, grid.ylast, grid.yinc, yvals);
	    gridDefYvals(gridID, yvals);
	    free(yvals);
	    /*
	    gridDefYinc(gridID, grid.yinc);
	    */
	  }

	if ( grid.isRotated )
	  {
	    gridDefXname(gridID, "rlon");
	    gridDefYname(gridID, "rlat");
	    gridDefXlongname(gridID, "longitude in rotated pole grid");
	    gridDefYlongname(gridID, "latitude in rotated pole grid");
	    strcpy(gridptr->xstdname, "grid_longitude");
	    strcpy(gridptr->ystdname, "grid_latitude");
	    gridDefXunits(gridID, "degrees");
	    gridDefYunits(gridID, "degrees");

	    gridDefXpole(gridID, grid.xpole);
	    gridDefYpole(gridID, grid.ypole);
	    gridDefAngle(gridID, grid.angle);
	  }

	if ( grid.area )
	  {
	    gridDefArea(gridID, grid.area);
	  }

	if ( grid.type == GRID_LAEA )
	  gridDefLaea(gridID, grid.laea_a, grid.laea_lon_0, grid.laea_lat_0);

	if ( grid.type == GRID_LCC2 )
	  gridDefLcc2(gridID, grid.lcc2_a, grid.lcc2_lon_0, grid.lcc2_lat_0, grid.lcc2_lat_1, grid.lcc2_lat_2);

	if ( grid.type == GRID_LCC )
	  gridDefLCC(gridID, grid.lcc_originLon, grid.lcc_originLat, grid.lcc_lonParY,
		     grid.lcc_lat1, grid.lcc_lat2, grid.lcc_xinc, grid.lcc_yinc,
		     grid.lcc_projflag, grid.lcc_scanflag);

	if ( grid.type == GRID_PROJECTION )
	  {
	    gridptr->name = strdup(grid.name);
	  }

	break;
      }
    case GRID_GAUSSIAN_REDUCED:
      {
	gridDefNP(gridID, grid.np);
	gridDefYsize(gridID, grid.ysize);
	gridDefRowlon(gridID, grid.ysize, grid.rowlon);

	if ( grid.ydef == 1 )
	  {
	    gridDefYvals(gridID, grid.yvals);
	    if ( grid.ybounds && grid.nvertex )
	      gridDefYbounds(gridID, grid.ybounds);
	  }
	else if ( grid.ydef == 2 )
	  {
	    double *yvals = (double *) malloc(grid.ysize*sizeof(double));
	    gridGenYvals(grid.type, grid.ysize, grid.yfirst, grid.ylast, grid.yinc, yvals);
	    gridDefYvals(gridID, yvals);
	    free(yvals);
	    /*
	    gridDefYinc(gridID, grid.yinc);
	    */
	  }
	break;
      }
    case GRID_SPECTRAL:
      {
	gridDefTrunc(gridID, grid.trunc);
	if ( grid.lcomplex ) gridDefComplexPacking(gridID, 1);
	break;
      }
    case GRID_FOURIER:
      {
	gridDefTrunc(gridID, grid.trunc);
	break;
      }
    case GRID_GME:
      {
	gridDefGMEnd(gridID, grid.nd);
	gridDefGMEni(gridID, grid.ni);
	gridDefGMEni2(gridID, grid.ni2);
	gridDefGMEni3(gridID, grid.ni3);
	break;
      }
    case GRID_REFERENCE:
      {
	gridDefNumber(gridID, grid.number);
	gridDefPosition(gridID, grid.position);
        gridDefUUID(gridID, grid.uuid);
	if ( grid.reference ) gridDefReference(gridID, grid.reference);
	break;
      }
      /*
    case GRID_GENERIC:
      {
	if ( grid.xsize > 0 && grid.ysize > 0 )
	  {
	    gridDefXsize(gridID, grid.xsize);
	    gridDefYsize(gridID, grid.ysize);
	    if ( grid.xvals ) gridDefXvals(gridID, grid.xvals);
	    if ( grid.yvals ) gridDefYvals(gridID, grid.yvals);
	  }
	break;
      }
      */
    case GRID_TRAJECTORY:
      {
	gridDefXsize(gridID, 1);
	gridDefYsize(gridID, 1);
	break;
      }
    default:
      {
	Error("Gridtype %s unsupported!", gridNamePtr(grid.type));
	break;
      }
    }

  if ( grid.xname[0]     ) gridDefXname(gridID, grid.xname);
  if ( grid.xlongname[0] ) gridDefXlongname(gridID, grid.xlongname);
  if ( grid.xunits[0]    ) gridDefXunits(gridID, grid.xunits);
  if ( grid.yname[0]     ) gridDefYname(gridID, grid.yname);
  if ( grid.ylongname[0] ) gridDefYlongname(gridID, grid.ylongname);
  if ( grid.yunits[0]    ) gridDefYunits(gridID, grid.yunits);

  return (gridID);
}


/*
@Function  gridDuplicate
@Title     Duplicate a horizontal Grid

@Prototype int gridDuplicate(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate},
                    @fref{gridDuplicate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridDuplicate} duplicates a horizontal Grid.

@Result
@func{gridDuplicate} returns an identifier to the duplicated Grid.

@EndFunction
*/
int gridDuplicate(int gridID)
{
  int gridIDnew;
  int gridtype, gridsize;
  int nrowlon;
  int size;
  grid_t *gridptr, *gridptrnew;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridtype = gridInqType(gridID);
  gridsize = gridInqSize(gridID);

  gridIDnew = gridCreate(gridtype, gridsize);
  gridptrnew = grid_to_pointer(gridIDnew);

  grid_copy(gridptrnew, gridptr);

  strcpy(gridptrnew->xname, gridptr->xname);
  strcpy(gridptrnew->yname, gridptr->yname);
  strcpy(gridptrnew->xlongname, gridptr->xlongname);
  strcpy(gridptrnew->ylongname, gridptr->ylongname);
  strcpy(gridptrnew->xunits, gridptr->xunits);
  strcpy(gridptrnew->yunits, gridptr->yunits);
  strcpy(gridptrnew->xstdname, gridptr->xstdname);
  strcpy(gridptrnew->ystdname, gridptr->ystdname);

  nrowlon = gridptr->nrowlon;
  if ( nrowlon )
    {
      gridptrnew->rowlon = (int *) malloc(nrowlon*sizeof(int));
      memcpy(gridptrnew->rowlon, gridptr->rowlon, nrowlon*sizeof(int));
    }

  if ( gridptr->xvals != NULL )
    {
      if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED )
	size = gridsize;
      else
	size = gridptr->xsize;

      gridptrnew->xvals = (double *) malloc(size*sizeof(double));
      memcpy(gridptrnew->xvals, gridptr->xvals, size*sizeof(double));
    }

  if ( gridptr->yvals != NULL )
    {
      if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED )
	size = gridsize;
      else
	size = gridptr->ysize;

      gridptrnew->yvals = (double *) malloc(size*sizeof(double));
      memcpy(gridptrnew->yvals, gridptr->yvals, size*sizeof(double));
    }

  if ( gridptr->xbounds != NULL )
    {
      if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED )
	size = gridsize;
      else
	size = gridptr->xsize;

      size *= gridptr->nvertex;

      gridptrnew->xbounds = (double *) malloc(size*sizeof(double));
      memcpy(gridptrnew->xbounds, gridptr->xbounds, size*sizeof(double));
    }

  if ( gridptr->ybounds != NULL )
    {
      if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED )
	size = gridsize;
      else
	size = gridptr->ysize;

      size *= gridptr->nvertex;

      gridptrnew->ybounds = (double *) malloc(size*sizeof(double));
      memcpy(gridptrnew->ybounds, gridptr->ybounds, size*sizeof(double));
    }

  if ( gridptr->area != NULL )
    {
      size = gridsize;

      gridptrnew->area = (double *) malloc(size*sizeof(double));
      memcpy(gridptrnew->area, gridptr->area, size*sizeof(double));
    }

  if ( gridptr->mask != NULL )
    {
      size = gridsize;

      gridptrnew->mask = (mask_t *) malloc(size*sizeof(mask_t));
      memcpy(gridptrnew->mask, gridptr->mask, size*sizeof(mask_t));
    }

  if ( gridptr->mask_gme != NULL )
    {
      size = gridsize;

      gridptrnew->mask_gme = (mask_t *) malloc(size*sizeof(mask_t));
      memcpy(gridptrnew->mask_gme, gridptr->mask_gme, size*sizeof(mask_t));
    }

  return (gridIDnew);
}


void gridCompress(int gridID)
{
  int gridtype, gridsize;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridtype = gridInqType(gridID);
  gridsize = gridInqSize(gridID);

  if ( gridtype == GRID_UNSTRUCTURED )
    {
      if ( gridptr->mask_gme != NULL )
	{
	  long i, j, iv, nv;

	  nv = gridptr->nvertex;

	  j = 0;
	  for ( i = 0; i < gridsize; i++ )
	    {
	      if ( gridptr->mask_gme[i] )
		{
		  if ( gridptr->xvals != NULL ) gridptr->xvals[j] = gridptr->xvals[i];
		  if ( gridptr->yvals != NULL ) gridptr->yvals[j] = gridptr->yvals[i];
		  if ( gridptr->area  != NULL ) gridptr->area[j]  = gridptr->area[i];
		  if ( gridptr->xbounds != NULL )
		    for ( iv = 0; iv < nv; iv++ )
		      gridptr->xbounds[j*nv+iv] = gridptr->xbounds[i*nv+iv];
		  if ( gridptr->ybounds != NULL )
		    for ( iv = 0; iv < nv; iv++ )
		      gridptr->ybounds[j*nv+iv] = gridptr->ybounds[i*nv+iv];

		  j++;
		}
	    }

	  /* fprintf(stderr, "grid compress %d %d %d\n", i, j, gridsize); */
	  gridsize = j;
	  gridptr->size  = gridsize;
	  gridptr->xsize = gridsize;
	  gridptr->ysize = gridsize;

	  if ( gridptr->xvals )
	    gridptr->xvals = (double *) realloc(gridptr->xvals, gridsize*sizeof(double));

	  if ( gridptr->yvals )
	    gridptr->yvals = (double *) realloc(gridptr->yvals, gridsize*sizeof(double));

	  if ( gridptr->area )
	    gridptr->area  = (double *) realloc(gridptr->area, gridsize*sizeof(double));

	  if ( gridptr->xbounds )
	    gridptr->xbounds = (double *) realloc(gridptr->xbounds, nv*gridsize*sizeof(double));

	  if ( gridptr->ybounds )
	    gridptr->ybounds = (double *) realloc(gridptr->ybounds, nv*gridsize*sizeof(double));

	  free(gridptr->mask_gme);
	  gridptr->mask_gme = NULL;
	}
    }
  else
    Warning("Unsupported grid type: %s", gridNamePtr(gridtype));
    
}


void gridDefArea(int gridID, const double *area)
{
  long size;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  size = gridptr->size;

  if ( size == 0 )
    Error("size undefined for gridID = %d", gridID);
    
  if ( gridptr->area == NULL )
    gridptr->area = (double *) malloc(size*sizeof(double));
  else if ( CDI_Debug )
    Warning("values already defined!");

  memcpy(gridptr->area, area, size*sizeof(double));
}


void gridInqArea(int gridID, double *area)
{
  long size;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  size = gridptr->size;

  if ( gridptr->area )
    memcpy(area, gridptr->area, size*sizeof(double));
}


int gridHasArea(int gridID)
{
  int hasArea = FALSE;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->area != NULL ) hasArea = TRUE;

  return (hasArea);
}


const double *gridInqAreaPtr(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->area);
}


void gridDefNvertex(int gridID, int nvertex)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridptr->nvertex = nvertex;
}


int gridInqNvertex(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->nvertex);
}


/*
@Function  gridDefXbounds
@Title     Define the bounds of a X-axis

@Prototype void gridDefXbounds(int gridID, const double *xbounds)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  xbounds  X-bounds of the grid.

@Description
The function @func{gridDefXbounds} defines all bounds of the X-axis.

@EndFunction
*/
void gridDefXbounds(int gridID, const double *xbounds)
{
  long size;
  long nvertex;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  nvertex = gridptr->nvertex;
  if ( nvertex == 0 )
    {
      Warning("nvertex undefined for gridID = %d. Cannot define bounds!", gridID);
      return;
    }

  if ( gridptr->type == GRID_CURVILINEAR || gridptr->type == GRID_UNSTRUCTURED )
    size = nvertex*gridptr->size;
  else
    size = nvertex*gridptr->xsize;

  if ( size == 0 )
    Error("size undefined for gridID = %d", gridID);
    
  if ( gridptr->xbounds == NULL )
    gridptr->xbounds = (double *) malloc(size*sizeof(double));
  else if ( CDI_Debug )
    Warning("values already defined!");

  memcpy(gridptr->xbounds, xbounds, size*sizeof(double));
}


/*
@Function  gridInqXbounds
@Title     Get the bounds of a X-axis

@Prototype int gridInqXbounds(int gridID, double *xbounds)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  xbounds  Pointer to the location into which the X-bounds are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{gridInqXbounds} returns the bounds of the X-axis.

@Result
Upon successful completion @func{gridInqXbounds} returns the number of bounds and
the bounds are stored in @func{xbounds}.
Otherwise, 0 is returned and @func{xbounds} is empty.

@EndFunction
*/
int gridInqXbounds(int gridID, double *xbounds)
{
  long size;
  long nvertex;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  nvertex = gridptr->nvertex;
  if ( CDI_Debug && nvertex == 0 )
    Warning("nvertex undefined for gridID = %d", gridID);

  if ( gridptr->type == GRID_CURVILINEAR || gridptr->type == GRID_UNSTRUCTURED )
    size = nvertex*gridptr->size;
  else
    size = nvertex*gridptr->xsize;

  if ( CDI_Debug && size == 0 )
    Warning("size undefined for gridID = %d", gridID);

  if ( xbounds && gridptr->xbounds )
    memcpy(xbounds, gridptr->xbounds, size*sizeof(double));

  if ( gridptr->xbounds == NULL ) size = 0;

  return ((int)size);
}


double *gridInqXboundsPtr(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->xbounds);
}


/*
@Function  gridDefYbounds
@Title     Define the bounds of a Y-axis

@Prototype void gridDefYbounds(int gridID, const double *ybounds)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  ybounds  Y-bounds of the grid.

@Description
The function @func{gridDefYbounds} defines all bounds of the Y-axis.

@EndFunction
*/
void gridDefYbounds(int gridID, const double *ybounds)
{
  long size;
  long nvertex;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  nvertex = gridptr->nvertex;
  if ( nvertex == 0 )
    {
      Warning("nvertex undefined for gridID = %d. Cannot define bounds!", gridID);
      return;
    }

  if ( gridptr->type == GRID_CURVILINEAR || gridptr->type == GRID_UNSTRUCTURED )
    size = nvertex*gridptr->size;
  else
    size = nvertex*gridptr->ysize;

  if ( size == 0 )
    Error("size undefined for gridID = %d", gridID);
    
  if ( gridptr->ybounds == NULL )
    gridptr->ybounds = (double *) malloc(size*sizeof(double));
  else if ( CDI_Debug )
    Warning("values already defined!");

  memcpy(gridptr->ybounds, ybounds, size*sizeof(double));
}


/*
@Function  gridInqYbounds
@Title     Get the bounds of a Y-axis

@Prototype int gridInqYbounds(int gridID, double *ybounds)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  ybounds  Pointer to the location into which the Y-bounds are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{gridInqYbounds} returns the bounds of the Y-axis.

@Result
Upon successful completion @func{gridInqYbounds} returns the number of bounds and
the bounds are stored in @func{ybounds}.
Otherwise, 0 is returned and @func{ybounds} is empty.

@EndFunction
*/
int gridInqYbounds(int gridID, double *ybounds)
{
  long size;
  long nvertex;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  nvertex = gridptr->nvertex;
  if ( CDI_Debug && nvertex == 0 )
    Warning("nvertex undefined for gridID = %d", gridID);

  if ( gridptr->type == GRID_CURVILINEAR || gridptr->type == GRID_UNSTRUCTURED )
    size = nvertex*gridptr->size;
  else
    size = nvertex*gridptr->ysize;

  if ( CDI_Debug && size == 0 )
    Warning("size undefined for gridID = %d", gridID);

  if ( ybounds && gridptr->ybounds )
    memcpy(ybounds, gridptr->ybounds, size*sizeof(double));

  if ( gridptr->ybounds == NULL ) size = 0;

  return ((int)size);
}


double *gridInqYboundsPtr(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->ybounds);
}


void gridPrint(int gridID, int opt)
{
  FILE *fp = stdout;
  int type;
  int gridsize, xsize, ysize, xdim, ydim;
  int trunc;
  int nbyte0, nbyte;
  int i;
  int nvertex, iv;
  char uuid[17];
  const double *area    = gridInqAreaPtr(gridID);
  const double *xvals   = gridInqXvalsPtr(gridID);
  const double *yvals   = gridInqYvalsPtr(gridID);
  const double *xbounds = gridInqXboundsPtr(gridID);
  const double *ybounds = gridInqYboundsPtr(gridID);
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  type     = gridInqType(gridID);
  trunc    = gridInqTrunc(gridID);
  gridsize = gridInqSize(gridID);
  xsize    = gridInqXsize(gridID);
  ysize    = gridInqYsize(gridID);
  nvertex  = gridInqNvertex(gridID);

  nbyte0 = 0;
  fprintf(fp, "#\n");
  fprintf(fp, "# gridID %d\n", gridID);
  fprintf(fp, "#\n");
  fprintf(fp, "gridtype  = %s\n", gridNamePtr(type));
  fprintf(fp, "gridsize  = %d\n", gridsize);

  if ( type != GRID_GME )
    {
      if ( gridptr->xname[0]     )     fprintf(fp, "xname     = %s\n", gridptr->xname);
      if ( gridptr->xlongname[0] )     fprintf(fp, "xlongname = %s\n", gridptr->xlongname);
      if ( gridptr->xunits[0]    )     fprintf(fp, "xunits    = %s\n", gridptr->xunits);
      if ( gridptr->yname[0]     )     fprintf(fp, "yname     = %s\n", gridptr->yname);
      if ( gridptr->ylongname[0] )     fprintf(fp, "ylongname = %s\n", gridptr->ylongname);
      if ( gridptr->yunits[0]    )     fprintf(fp, "yunits    = %s\n", gridptr->yunits);
      if ( type == GRID_UNSTRUCTURED ) fprintf(fp, "nvertex   = %d\n", nvertex);
    }

  switch (type)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GENERIC:
    case GRID_LCC2:
    case GRID_SINUSOIDAL:
    case GRID_LAEA:
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
      {
        if ( type == GRID_GAUSSIAN ) fprintf(fp, "np        = %d\n", gridptr->np);

	if ( type == GRID_CURVILINEAR || type == GRID_UNSTRUCTURED )
	  {
	    xdim = gridsize;
	    ydim = gridsize;
	  }
	else
	  {
	    xdim = xsize;
	    ydim = ysize;
	  }

	if ( type != GRID_UNSTRUCTURED )
	  {
	    if ( xsize > 0 ) fprintf(fp, "xsize     = %d\n", xsize);
	    if ( ysize > 0 ) fprintf(fp, "ysize     = %d\n", ysize);
	  }

	if ( type == GRID_LAEA )
	  {
	    double a, lon_0, lat_0;
	    gridInqLaea(gridID, &a, &lon_0, &lat_0);
	    fprintf(fp, "a         = %g\n", a);
	    fprintf(fp, "lon_0     = %g\n", lon_0);
	    fprintf(fp, "lat_0     = %g\n", lat_0);
	  }

	if ( type == GRID_LCC2 )
	  {
	    double a, lon_0, lat_0, lat_1, lat_2;
	    gridInqLcc2(gridID, &a, &lon_0, &lat_0, &lat_1, &lat_2);
	    fprintf(fp, "a         = %g\n", a);
	    fprintf(fp, "lon_0     = %g\n", lon_0);
	    fprintf(fp, "lat_0     = %g\n", lat_0);
	    fprintf(fp, "lat_1     = %g\n", lat_1);
	    fprintf(fp, "lat_2     = %g\n", lat_2);
	  }

	if ( gridptr->isRotated )
	  {
	    if ( xsize > 0 ) fprintf(fp, "xnpole    = %g\n", gridptr->xpole);
	    if ( ysize > 0 ) fprintf(fp, "ynpole    = %g\n", gridptr->ypole);
	    if ( gridptr->angle > 0 ) fprintf(fp, "angle     = %g\n", gridptr->angle);
 	  }

	if ( xvals )
	  {
	    double xfirst = 0.0, xinc = 0.0;

	    if ( type == GRID_LONLAT     || type == GRID_GAUSSIAN || 
		 type == GRID_GENERIC    || type == GRID_LCC2     || 
                 type == GRID_SINUSOIDAL || type == GRID_LAEA )
	      {
		xfirst = gridInqXval(gridID, 0);
		xinc   = gridInqXinc(gridID);
	      }

	    if ( IS_NOT_EQUAL(xinc, 0) && opt )
	      {
	  	fprintf(fp, "xfirst    = %g\n", xfirst);
		fprintf(fp, "xinc      = %g\n", xinc);
	      }
	    else
	      {
		nbyte0 = fprintf(fp, "xvals     = ");
		nbyte = nbyte0;
		for ( i = 0; i < xdim; i++ )
		  {
		    if ( nbyte > 80 )
		      {
			fprintf(fp, "\n");
			fprintf(fp, "%*s", nbyte0, "");
			nbyte = nbyte0;
		      }
		    nbyte += fprintf(fp, "%.9g ", xvals[i]);
		  }
		fprintf(fp, "\n");
	      }
	  }

	if ( xbounds )
	  {
	    nbyte0 = fprintf(fp, "xbounds   = ");
	    for ( i = 0; i < xdim; i++ )
	      {
		if ( i ) fprintf(fp, "%*s", nbyte0, "");

		for ( iv = 0; iv < nvertex; iv++ )
		  fprintf(fp, "%.9g ", xbounds[i*nvertex+iv]);
		fprintf(fp, "\n");
	      }
	  }

	if ( yvals )
	  {
	    double yfirst = 0.0, yinc = 0.0;

	    if ( type == GRID_LONLAT || type == GRID_GENERIC || type == GRID_LCC2 ||
		 type == GRID_SINUSOIDAL || type == GRID_LAEA )
	      {
		yfirst = gridInqYval(gridID, 0);
		yinc   = gridInqYinc(gridID);
	      }

	    if ( IS_NOT_EQUAL(yinc, 0) && opt )
	      {
	  	fprintf(fp, "yfirst    = %g\n", yfirst);
		fprintf(fp, "yinc      = %g\n", yinc);
	      }
	    else
	      {
		nbyte0 = fprintf(fp, "yvals     = ");
		nbyte = nbyte0;
		for ( i = 0; i < ydim; i++ )
		  {
		    if ( nbyte > 80 )
		      {
			fprintf(fp, "\n");
			fprintf(fp, "%*s", nbyte0, "");
			nbyte = nbyte0;
		      }
		    nbyte += fprintf(fp, "%.9g ", yvals[i]);
		  }
		fprintf(fp, "\n");
	      }
	  }

	if ( ybounds )
	  {
	    nbyte0 = fprintf(fp, "ybounds   = ");
	    for ( i = 0; i < ydim; i++ )
	      {
		if ( i ) fprintf(fp, "%*s", nbyte0, "");

		for ( iv = 0; iv < nvertex; iv++ )
		  fprintf(fp, "%.9g ", ybounds[i*nvertex+iv]);
		fprintf(fp, "\n");
	      }
	  }

	if ( area )
	  {
	    nbyte0 = fprintf(fp, "area      = ");
	    nbyte  = nbyte0;
	    for ( i = 0; i < gridsize; i++ )
	      {
		if ( nbyte > 80 )
		  {
		    fprintf(fp, "\n");
		    fprintf(fp, "%*s", nbyte0, "");
		    nbyte = nbyte0;
		  }
		nbyte += fprintf(fp, "%.9g ", area[i]);
	      }
	    fprintf(fp, "\n");
	  }
	break;
      }
   case GRID_GAUSSIAN_REDUCED:
      {
	int *rowlon;
	fprintf(fp, "ysize = %d\n", ysize);
	nbyte0 = fprintf(fp, "rowlon = %d  ", ysize);
	nbyte  = nbyte0;
	rowlon = (int *) malloc(ysize*sizeof(int));
	gridInqRowlon(gridID, rowlon);
	for ( i = 0; i < ysize; i++ )
	  {
	    if ( nbyte > 80 )
	      {
		fprintf(fp, "\n");
		fprintf(fp, "%*s", nbyte0, "");
		nbyte = nbyte0;
	      }
	    nbyte += fprintf(fp, "%d ", rowlon[i]);
	  }
	fprintf(fp, "\n");
	free(rowlon);
	break;
      }
    case GRID_LCC:
      {
	double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
	int projflag, scanflag;
	gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
		   &projflag, &scanflag);

	fprintf(fp, "xsize     = %d\n", xsize);
	fprintf(fp, "ysize     = %d\n", ysize);

	fprintf(fp, "originLon = %g\n", originLon);
	fprintf(fp, "originLat = %g\n", originLat);
	fprintf(fp, "lonParY   = %g\n", lonParY);
	fprintf(fp, "lat1      = %g\n", lat1);
	fprintf(fp, "lat2      = %g\n", lat2);
	fprintf(fp, "xinc      = %g\n", xincm);
	fprintf(fp, "yinc      = %g\n", yincm);
	if ( (projflag & 128) == 0 )
	  fprintf(fp, "projection = northpole\n");
	else
	  fprintf(fp, "projection = southpole\n");

	break;
      }
    case GRID_SPECTRAL:
      {
	fprintf(fp, "truncation = %d\n", trunc);
	fprintf(fp, "complexpacking = %d\n", gridInqComplexPacking(gridID));
	break;
      }
    case GRID_FOURIER:
      {
	fprintf(fp, "truncation = %d\n", trunc);
	break;
      }
    case GRID_GME:
      {
	fprintf(fp, "ni        = %d\n", gridInqGMEni(gridID));
	break;
      }
    case GRID_REFERENCE:
      {
        const unsigned char *d;
	fprintf(fp, "number    = %d\n", gridInqNumber(gridID));
	fprintf(fp, "position  = %d\n", gridInqPosition(gridID));
        gridInqUUID(gridID, uuid);
        d = (unsigned char *) &uuid;
	fprintf(fp, "uuid      = %02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x\n", 
                d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7],
                d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15]);
	if ( gridInqReference(gridID, NULL) )
	  {
	    char reference_link[8192];
	    gridInqReference(gridID, reference_link);
	    fprintf(fp, "path      = %s\n", reference_link);
	  }
	break;
      }
   default:
      {
	fprintf(stderr, "Unsupported grid type: %s\n", gridNamePtr(type));
      }
    }

  if ( gridptr->mask )
    {
      nbyte0 = fprintf(fp, "mask      = ");
      nbyte  = nbyte0;
      for ( i = 0; i < gridsize; i++ )
	{
	  if ( nbyte > 80 )
	    {
	      fprintf(fp, "\n");
	      fprintf(fp, "%*s", nbyte0, "");
	      nbyte = nbyte0;
	    }
	  nbyte += fprintf(fp, "%d ", (int) gridptr->mask[i]);
	}
      fprintf(fp, "\n");
    }
}


const double *gridInqXvalsPtr(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return ( gridptr->xvals );
}


const double *gridInqYvalsPtr(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return ( gridptr->yvals );
}


/*
@Function  gridDefLCC
@Title     Define the parameter of a Lambert Conformal Conic grid

@Prototype void gridDefLCC(int gridID, double originLon, double originLat, double lonParY, double lat1, double lat2, double xinc, double yinc, int projflag, int scanflag)
@Parameter
    @Item  gridID    Grid ID, from a previous call to @fref{gridCreate}.
    @Item  originLon Longitude of the first grid point.
    @Item  originLat Latitude of the first grid point.
    @Item  lonParY   The East longitude of the meridian which is parallel to the Y-axis.
    @Item  lat1      First latitude from the pole at which the secant cone cuts the sphere.
    @Item  lat2      Second latitude at which the secant cone cuts the sphere.
    @Item  xinc      X-direction grid lenght in meter.
    @Item  yinc      Y-direction grid lenght in meter.
    @Item  projflag  Projection centre flag.
    @Item  scanflag  Scanning mode flag.
 
@Description
The function @func{gridDefLCC} defines the parameter of a Lambert Conformal Conic grid.

@EndFunction
*/
void gridDefLCC(int gridID, double originLon, double originLat, double lonParY,
		double lat1, double lat2, double xinc, double yinc,
		int projflag, int scanflag)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->type != GRID_LCC )
    Warning("Definition of LCC grid for %s grid not allowed!",
	    gridNamePtr(gridptr->type));
  else
    {
      gridptr->lcc_originLon = originLon;
      gridptr->lcc_originLat = originLat;
      gridptr->lcc_lonParY   = lonParY;
      gridptr->lcc_lat1      = lat1;
      gridptr->lcc_lat2      = lat2;
      gridptr->lcc_xinc      = xinc;
      gridptr->lcc_yinc      = yinc;
      gridptr->lcc_projflag  = projflag;
      gridptr->lcc_scanflag  = scanflag;
      gridptr->lcc_defined   = TRUE;
    }
}


/*
@Function  gridInqLCC
@Title     Get the parameter of a Lambert Conformal Conic grid

@Prototype void gridInqLCC(int gridID, double *originLon, double *originLat, double *lonParY, double *lat1, double *lat2, double *xinc, double *yinc, int *projflag, int *scanflag)
@Parameter
    @Item  gridID    Grid ID, from a previous call to @fref{gridCreate}.
    @Item  originLon Longitude of the first grid point.
    @Item  originLat Latitude of the first grid point.
    @Item  lonParY   The East longitude of the meridian which is parallel to the Y-axis.
    @Item  lat1      First latitude from the pole at which the secant cone cuts the sphere.
    @Item  lat2      Second latitude at which the secant cone cuts the sphere.
    @Item  xinc      X-direction grid lenght in meter.
    @Item  yinc      Y-direction grid lenght in meter.
    @Item  projflag  Projection centre flag.
    @Item  scanflag  Scanning mode flag.
 
@Description
The function @func{gridInqLCC} returns the parameter of a Lambert Conformal Conic grid.

@EndFunction
*/
void gridInqLCC(int gridID, double *originLon, double *originLat, double *lonParY,
		double *lat1, double *lat2, double *xinc, double *yinc,
		int *projflag, int *scanflag)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->type != GRID_LCC )
    Warning("Inquire of LCC grid definition for %s grid not allowed!",
	    gridNamePtr(gridptr->type));
  else
    {
      if ( gridptr->lcc_defined )
	{
	  *originLon = gridptr->lcc_originLon;
	  *originLat = gridptr->lcc_originLat;
	  *lonParY   = gridptr->lcc_lonParY;
	  *lat1      = gridptr->lcc_lat1;
	  *lat2      = gridptr->lcc_lat2;
	  *xinc      = gridptr->lcc_xinc;
	  *yinc      = gridptr->lcc_yinc;
	  *projflag  = gridptr->lcc_projflag;
	  *scanflag  = gridptr->lcc_scanflag;
	}
      else
	Warning("Lambert Conformal grid undefined (gridID = %d)", gridID);
    }
}

void gridDefLcc2(int gridID, double earth_radius, double lon_0, double lat_0, double lat_1, double lat_2)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->type != GRID_LCC2 )
    Warning("Definition of LCC2 grid for %s grid not allowed!",
	    gridNamePtr(gridptr->type));
  else
    {
      gridptr->lcc2_a       = earth_radius;
      gridptr->lcc2_lon_0   = lon_0;
      gridptr->lcc2_lat_0   = lat_0;
      gridptr->lcc2_lat_1   = lat_1;
      gridptr->lcc2_lat_2   = lat_2;
      gridptr->lcc2_defined = TRUE;
    }
}


void gridInqLcc2(int gridID, double *earth_radius, double *lon_0, double *lat_0, double *lat_1, double *lat_2)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->type != GRID_LCC2 )
    Warning("Inquire of LCC2 grid definition for %s grid not allowed!",
	    gridNamePtr(gridptr->type));
  else
    {
      if ( gridptr->lcc2_defined )
	{
	  *earth_radius = gridptr->lcc2_a;
	  *lon_0        = gridptr->lcc2_lon_0;
	  *lat_0        = gridptr->lcc2_lat_0;
	  *lat_1        = gridptr->lcc2_lat_1;
	  *lat_2        = gridptr->lcc2_lat_2;
	}
      else
	Warning("LCC2 grid undefined (gridID = %d)", gridID);
    }
}

void gridDefLaea(int gridID, double earth_radius, double lon_0, double lat_0)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->type != GRID_LAEA )
    Warning("Definition of LAEA grid for %s grid not allowed!",
	    gridNamePtr(gridptr->type));
  else
    {
      gridptr->laea_a       = earth_radius;
      gridptr->laea_lon_0   = lon_0;
      gridptr->laea_lat_0   = lat_0;
      gridptr->laea_defined = TRUE;
    }
}


void gridInqLaea(int gridID, double *earth_radius, double *lon_0, double *lat_0)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->type != GRID_LAEA )
    Warning("Inquire of LAEA grid definition for %s grid not allowed!",
	    gridNamePtr(gridptr->type));
  else
    {
      if ( gridptr->laea_defined )
	{
	  *earth_radius = gridptr->laea_a;
	  *lon_0        = gridptr->laea_lon_0;
	  *lat_0        = gridptr->laea_lat_0;
	}
      else
	Warning("LAEA grid undefined (gridID = %d)", gridID);
    }
}


void gridDefComplexPacking(int gridID, int lcomplex)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridptr->lcomplex = lcomplex;
}


int gridInqComplexPacking(int gridID)
{
  int lcomplex;
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  lcomplex = gridptr->lcomplex;

  return (lcomplex);
}


int gridInqNumber(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->number);
}


void gridDefNumber(int gridID, int number)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridptr->number = number;
}


int gridInqPosition(int gridID)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  return (gridptr->position);
}


void gridDefPosition(int gridID, int position)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  gridptr->position = position;
}


int gridInqReference(int gridID, char *reference)
{
  grid_t *gridptr;
  int len = 0;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( gridptr->reference )
    {
      len = (int) strlen(gridptr->reference);

      if ( reference )
	strcpy(reference, gridptr->reference);
    }

  return (len);
}


void gridDefReference(int gridID, const char *reference)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  if ( reference )
    {
      if ( gridptr->reference )
	{
	  free(gridptr->reference);
	  gridptr->reference = NULL;
	}

      gridptr->reference = strdupx(reference);
    }
}


char *gridInqUUID(int gridID, char *uuid)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  strncpy(uuid, gridptr->uuid, 16);

  return (uuid);
}


void gridDefUUID(int gridID, const char *uuid)
{
  grid_t *gridptr;

  gridptr = grid_to_pointer(gridID);

  grid_check_ptr(gridID, gridptr);

  strncpy(gridptr->uuid, uuid, 16);

  return;
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
