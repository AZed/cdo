#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif


#include <ctype.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#include "gradsdeslib.h"

static
void get_dim_vals(dsets_t *pfi, double *vals, int dimlen, int dim)
{
  gadouble (*conv) (gadouble *, gadouble);
  gadouble *cvals;
  int i;

  assert( dimlen == pfi->dnum[dim] );

  if ( pfi->linear[dim] == 0 )
    {
      for ( i = 0; i < dimlen; ++i )
	{
	  vals[i] = pfi->grvals[dim][i+1];
	  /* printf("%d %g\n", i, vals[i]); */
	}
    }
  else if ( pfi->linear[dim] == 1 )
    {
      /*
      for ( i = 0; i < 3; ++i )
	printf("%d %g %g\n", i, pfi->grvals[dim][i] , pfi->abvals[dim][i]);
      */
      conv = pfi->gr2ab[dim];
      cvals = pfi->grvals[dim];
      for ( i = 0; i < dimlen; ++i )
	{
	  vals[i] = conv(cvals, i+1);
	  /* printf("%d %g\n", i, vals[i]); */
	}
    }
  
}


static
void rev_yvals(double *yvals, int ny)
{
  int i;
  double dum;

  for ( i = 0; i < ny/2; ++i )
    {
      dum = yvals[i];
      yvals[i] = yvals[ny-1-i];
      yvals[ny-1-i] = dum;
    }
}


static
int y_is_gauss(double *gridyvals, int ysize)
{
  static char func[] = "y_is_gauss";
  int lgauss = FALSE;
  int i;

  if ( ysize > 2 )
    {
      double *yvals, *yw;
      yvals = (double *) malloc(ysize*sizeof(double));
      yw    = (double *) malloc(ysize*sizeof(double));
      gaussaw(yvals, yw, ysize);
      free(yw);
      for ( i = 0; i < (int) ysize; i++ )
	yvals[i] = asin(yvals[i])/M_PI*180.0;

      for ( i = 0; i < (int) ysize; i++ )
	if ( fabs(yvals[i] - gridyvals[i]) >
	     ((yvals[0] - yvals[1])/500) ) break;
		      
      if ( i == (int) ysize ) lgauss = TRUE;

      /* check S->N */
      if ( lgauss == FALSE )
	{		  
	  for ( i = 0; i < (int) ysize; i++ )
	    if ( fabs(yvals[i] - gridyvals[ysize-i-1]) >
		 ((yvals[0] - yvals[1])/500) ) break;
		      
	  if ( i == (int) ysize ) lgauss = TRUE;
	}

      free(yvals);
    }

  return (lgauss);
}


static
int define_grid(dsets_t *pfi)
{
  static char func[] = "define_grid";
  int gridID, gridtype;
  int nx, ny;
  double *xvals, *yvals;
  int lgauss;

  nx = pfi->dnum[0];
  ny = pfi->dnum[1];

  xvals = (double *) malloc(nx*sizeof(double));
  yvals = (double *) malloc(ny*sizeof(double));

  get_dim_vals(pfi, xvals, nx, 0);
  get_dim_vals(pfi, yvals, ny, 1);

  if ( pfi->yrflg ) rev_yvals(yvals, ny);

  lgauss = y_is_gauss(yvals, ny);

  if ( lgauss ) gridtype = GRID_GAUSSIAN;
  else          gridtype = GRID_LONLAT;

  gridID = gridCreate(gridtype, nx*ny);
  gridDefXsize(gridID, nx);
  gridDefYsize(gridID, ny);

  gridDefXvals(gridID, xvals);
  gridDefYvals(gridID, yvals);

  free(xvals);
  free(yvals);
  
  return (gridID);
}


static
int define_level(dsets_t *pfi)
{
  static char func[] = "define_level";
  int zaxisID = -1;
  int nz;

  nz = pfi->dnum[2];

  if ( nz )
    {
      double *zvals = NULL;

      zvals = (double *) malloc(nz*sizeof(double));

      get_dim_vals(pfi, zvals, nz, 2);

      if ( nz == 1 && IS_EQUAL(zvals[0], 0) )
	zaxisID = zaxisCreate(ZAXIS_SURFACE, nz);
      else
	zaxisID = zaxisCreate(ZAXIS_GENERIC, nz);
      zaxisDefLevels(zaxisID, zvals);

      free(zvals);
    }
  else
    {
      double level = 0;
      nz = 1;

      zaxisID = zaxisCreate(ZAXIS_SURFACE, nz);
      zaxisDefLevels(zaxisID, &level);
    }

  
  return (zaxisID);
}


void *Importgrads(void *argument)
{
  static char func[] = "Importgrads";
  int streamID;
  int gridID = -1, zaxisID, zaxisIDsfc, taxisID, vlistID;
  int i;
  int nmiss, n_nan;
  int ivar;
  int varID, levelID, tsID;
  int gridsize;
  int  status;
  dsets_t pfi;
  int vdate, vtime;
  int tcur, told,fnum;
  int tmin=0,tmax=0;
  char *ch=NULL;
  off_t flen;
  int nvars, nlevels, nrecs;
  int recID;
  int e, flag;
  size_t rc, recsize;
  int recoffset;
  char *rec = NULL;
  struct gavar *pvar;
  struct dt dtim, dtimi;
  double fmin, fmax;
  float *farray;
  double *array;
  double sfclevel = 0;
  int *recVarID, *recLevelID;

  cdoInitialize(argument);

  dsets_init(&pfi);

  status = read_gradsdes((char *)cdoStreamName(0), &pfi);
  if ( cdoVerbose ) fprintf(stderr, "status %d\n", status);
  if ( status ) cdoAbort("Open failed on %s!", pfi.name);

  nrecs = pfi.trecs;
  nvars = pfi.vnum;
  pvar  = pfi.pvar1;

  if ( nvars == 0 ) cdoAbort("No variables found!");

  gridID = define_grid(&pfi);
  if ( cdoVerbose ) gridPrint(gridID, 1);

  zaxisID = define_level(&pfi);
  if ( cdoVerbose ) zaxisPrint(zaxisID);

  zaxisIDsfc = zaxisCreate(ZAXIS_SURFACE, 1);
  zaxisDefLevels(zaxisIDsfc, &sfclevel);

  vlistID = vlistCreate();

  recVarID   = (int *) malloc(nrecs*sizeof(int));
  recLevelID = (int *) malloc(nrecs*sizeof(int));

  recID = 0;
  for ( ivar = 0; ivar < nvars; ++ivar )
    {
      /*
      if ( cdoVerbose )
	fprintf(stderr, "1:%s 2:%s %d %d %d %d 3:%s %d \n", 
		pvar->abbrv, pvar->longnm, pvar->offset, pvar->recoff, pvar->levels, 
		pvar->nvardims, pvar->varnm, pvar->var_t);
      */
      nlevels = pvar->levels;
      
      if ( nlevels == 0 )
	{
	  nlevels = 1;
	  varID = vlistDefVar(vlistID, gridID, zaxisIDsfc, TIME_VARIABLE);
	}
      else
	{
	  if ( nlevels != zaxisInqSize(zaxisID) ) cdoAbort("Number of levels differ!");
	  varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
	}

      vlistDefVarName(vlistID, varID, pvar->abbrv);
      vlistDefVarLongname(vlistID, varID, pvar->varnm);
      vlistDefVarDatatype(vlistID, varID, DATATYPE_FLT32);
      vlistDefVarMissval(vlistID, varID, pfi.undef);

      for ( levelID = 0; levelID < nlevels; ++levelID )
	{
	  if ( recID >= nrecs ) cdoAbort("Internal problem with number of records!");
	  recVarID[recID]   = varID;
	  recLevelID[recID] = levelID;
	  recID++;
	}

      pvar++;
    }

  taxisID = taxisCreate(TAXIS_RELATIVE);

  taxisDefCalendar(taxisID, CALENDAR_STANDARD);

  vlistDefTaxis(vlistID, taxisID);

  streamID = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID, vlistID);


  gridsize = pfi.dnum[0]*pfi.dnum[1];
  recoffset = pfi.xyhdr*4;
  if ( pfi.seqflg ) recoffset += 4;

  recsize = pfi.gsiz*4;
  rec = (char *) malloc(recsize);

  array = (double *) malloc(gridsize*sizeof(double));

  /*
  if (pfi.tmplat)
    for ( i = 0; i <  pfi.dnum[3]; ++i )
      printf("%d %d\n", i, pfi.fnums[i]);
  */

  pfi.infile = NULL;
  tcur = 0;
  e = 1;
  while (1)
    {    /* loop over all times for this ensemble */
      if (pfi.tmplat)
	{
	  /* make sure no file is open */
	  if (pfi.infile!=NULL) {
	    fclose(pfi.infile);
	    pfi.infile=NULL;
	  }
	  /* advance to first valid time step for this ensemble */
	  if (tcur==0) {
	    told = 0;
	    tcur = 1;
	    while (pfi.fnums[tcur-1] == -1) tcur++;  
	  }
	  else {  /* tcur!=0 */
	    told = pfi.fnums[tcur-1];
	    /* increment time step until fnums changes */
	    while (told==pfi.fnums[tcur-1] && tcur<=pfi.dnum[3]) {
	      tcur++;
	      if ( tcur > pfi.dnum[3] ) break;
	    }
	  }

	  /* make sure we haven't advanced past end of time axis */
	  if (tcur>pfi.dnum[3]) break;

	  /* check if we're past all valid time steps for this ensemble */
	  if ((told != -1) && (pfi.fnums[tcur-1] == -1)) break;

	  /* Find the range of t indexes that have the same fnums value.
	     These are the times that are contained in this particular file */
	  tmin = tcur;
	  tmax = tcur-1;
	  fnum = pfi.fnums[tcur-1];
	  if (fnum != -1) {
	    while (fnum == pfi.fnums[tmax])
	      {
		tmax++; 
		if (tmax == pfi.dnum[3]) break;
	      }
	    gr2t(pfi.grvals[3], (gadouble)tcur, &dtim); 
	    gr2t(pfi.grvals[3], (gadouble)1, &dtimi);
	    ch = gafndt(pfi.name, &dtim, &dtimi, pfi.abvals[3], pfi.pchsub1, NULL,tcur,e,&flag);
	    if (ch==NULL) cdoAbort(" couldn't determine data file name for e=%d t=%d\n",e,tcur);
	  }
	}
      else { 
	/* Data set is not templated */
	ch = pfi.name;
	tmin = 1;
	tmax = pfi.dnum[3];
      }
       
      /* Open this file and position to start of first record */
      if ( cdoVerbose) cdoPrint("Opening file: %s", ch);
      pfi.infile = fopen(ch,"rb");
      if (pfi.infile==NULL) {
	if (pfi.tmplat) {
	  if ( cdoVerbose ) cdoPrint("Could not open file: %s",ch);
	  break;
	} else {
	  cdoAbort("Could not open file: %s",ch);
	}
      }
      if (pfi.tmplat) gree(ch,"312");
       
      /* Get file size */
      /*
      fseeko(pfi.infile,0L,2);
      flen = ftello(pfi.infile);

      printf("flen %d tsiz %d\n", flen, pfi.tsiz);
       
      fseeko (pfi.infile,0,0);
      */
      for ( tsID = tmin-1; tsID < tmax; ++tsID )
	{
	  gr2t(pfi.grvals[3], (gadouble)(tsID+1), &dtim); 
	  if ( cdoVerbose )
	    cdoPrint(" Reading timestep: %3d  %4.4d-%2.2d-%2.2d %2.2d:%2.2d",
		     tsID+1, dtim.yr, dtim.mo, dtim.dy, dtim.hr, dtim.mn);

	  vdate = encode_date(dtim.yr, dtim.mo, dtim.dy);
	  vtime = encode_time(dtim.hr, dtim.mn);
	  taxisDefVdate(taxisID, vdate);
	  taxisDefVtime(taxisID, vtime);
	  streamDefTimestep(streamID, tsID);

	  for ( recID = 0; recID < nrecs; ++recID )
	    {
	      rc = fread (rec, 1, recsize, pfi.infile);
	      if ( rc < recsize ) cdoAbort("I/O error reading record!");

	      if ( pfi.bswap ) gabswp(rec+recoffset, gridsize);
	      farray = (float *) (rec+recoffset);
	      fmin =  1.e99;
	      fmax = -1.e99;
	      nmiss = 0;
	      n_nan = 0;
	      for ( i = 0; i < gridsize; ++i )
		{
		  array[i] = (double) farray[i];
		  if ( array[i] > pfi.ulow && array[i] < pfi.uhi )
		    {
		      array[i] = pfi.undef;
		      nmiss++;
		    }
		  else if ( DBL_IS_NAN(array[i]) )
		    {
		      array[i] = pfi.undef;
		      nmiss++;
		      n_nan++;
		      /* printf("Nan at %d\n", i); */
		    }
		  else
		    {
		      if ( array[i] < fmin ) fmin = array[i];
		      if ( array[i] > fmax ) fmax = array[i];
		    }
		}
	      /*
	      if ( cdoVerbose )
		printf("%d %d %g %g %d %d %d\n", tsID, recID, fmin, fmax, recoffset, nmiss, n_nan);
	      */
	      varID   = recVarID[recID];
	      levelID = recLevelID[recID];
	      streamDefRecord(streamID,  varID,  levelID);
              streamWriteRecord(streamID, array, nmiss);
 	    }
	}

      /* break out if not templating */
      if (!pfi.tmplat) break;
      
    } /* end of while (1) loop */


  processDefVarNum(vlistNvars(vlistID), streamID);

  streamClose(streamID);

  vlistDestroy(vlistID);
  gridDestroy(gridID);
  zaxisDestroy(zaxisID);
  taxisDestroy(taxisID);

  free(array);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  cdoFinish();

  return (0);
}
