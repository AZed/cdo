#define DATE_FORMAT "%5.4d-%2.2d-%2.2d"
#define TIME_FORMAT "%2.2d:%2.2d:%2.2d"

void uuid2str(const char *uuid, char *uuidstr);

void date2str(int date, char *datestr, int maxlen)
{
  int year, month, day;
  int len;

  cdiDecodeDate(date, &year, &month, &day);

  len = sprintf(datestr, DATE_FORMAT, year, month, day);

  if ( len > ( maxlen-1) )
    fprintf(stderr, "Internal problem (%s): sizeof input string is too small!\n", __func__);
}


void time2str(int time, char *timestr, int maxlen)
{
  int hour, minute, second;
  int len;

  cdiDecodeTime(time, &hour, &minute, &second);

  len = sprintf(timestr, TIME_FORMAT, hour, minute, second);

  if ( len > ( maxlen-1) )
    fprintf(stderr, "Internal problem (%s): sizeof input string is too small!\n", __func__);
}


void printFiletype(int streamID, int vlistID)
{
  int filetype;

  filetype = streamInqFiletype(streamID);

  switch ( filetype )
    {
    case FILETYPE_GRB:
      printf("GRIB");
      break;
    case FILETYPE_GRB2:
      printf("GRIB2");
      break;
    case FILETYPE_NC:
      printf("netCDF");
      break;
    case FILETYPE_NC2:
      printf("netCDF2");
      break;
    case FILETYPE_NC4:
      printf("netCDF4");
      break;
    case FILETYPE_NC4C:
      printf("netCDF4 classic");
      break;
    case FILETYPE_SRV:
      printf("SERVICE");
      break;
    case FILETYPE_EXT:
      printf("EXTRA");
      break;
    case FILETYPE_IEG:
      printf("IEG");
      break;
    default:
      printf("  File format: unsupported filetype %d" , filetype);
      break;
    }

  if ( filetype == FILETYPE_SRV || filetype == FILETYPE_EXT || filetype == FILETYPE_IEG )
    {
      switch ( streamInqByteorder(streamID) )
	{
	case CDI_BIGENDIAN:
	  printf("  BIGENDIAN"); break;
	case CDI_LITTLEENDIAN:
	  printf("  LITTLEENDIAN"); break;
	default:
	  printf("  byteorder: %d undefined", streamInqByteorder(streamID)); break;
	}
    }

  if ( filetype == FILETYPE_GRB || filetype == FILETYPE_NC4 || filetype == FILETYPE_NC4C )
    {
      int nvars, varID;
      int comptype;

      nvars = vlistNvars(vlistID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  comptype = vlistInqVarCompType(vlistID, varID);
	  if ( comptype )
	    {
	      if ( comptype == COMPRESS_SZIP )
		printf(" SZIP");
	      else if ( comptype == COMPRESS_ZIP )
		printf(" ZIP");

	      break;
	    }
	}
    }

  if ( filetype == FILETYPE_GRB2 )
    {
      int nvars, varID;
      int comptype;

      nvars = vlistNvars(vlistID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  comptype = vlistInqVarCompType(vlistID, varID);
	  if ( comptype )
	    {
	      if ( comptype == COMPRESS_JPEG )
		printf(" JPEG");

	      break;
	    }
	}
    }

  printf("\n");
}

static
void printGridInfo(int vlistID)
{
  int ngrids, index;
  int gridID, gridtype, trunc, gridsize, xsize, ysize;
  int nbyte0;
  char xname[CDI_MAX_NAME], yname[CDI_MAX_NAME], xunits[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
  char uuidOfHGrid[17];

  ngrids = vlistNgrids(vlistID);
  for ( index = 0; index < ngrids; index++ )
    {
      gridID   = vlistGrid(vlistID, index);
      gridtype = gridInqType(gridID);
      trunc    = gridInqTrunc(gridID);
      gridsize = gridInqSize(gridID);
      xsize    = gridInqXsize(gridID);
      ysize    = gridInqYsize(gridID);
      gridInqXname(gridID, xname);
      gridInqYname(gridID, yname);
      gridInqXunits(gridID, xunits);
      gridInqYunits(gridID, yunits);

      nbyte0   = fprintf(stdout, "  %4d : %-12s > ", index+1, gridNamePtr(gridtype));

      if ( gridtype == GRID_LONLAT   ||
	   gridtype == GRID_LCC2 ||
	   gridtype == GRID_LAEA ||
	   gridtype == GRID_SINUSOIDAL ||
	   gridtype == GRID_GENERIC ||
	   gridtype == GRID_GAUSSIAN ||
	   gridtype == GRID_GAUSSIAN_REDUCED )
	{
          int lxcoord = 1, lycoord = 1;
	  double xfirst = 0.0, xlast = 0.0;
	  double yfirst = 0.0, ylast = 0.0;
	  double xinc = 0.0, yinc = 0.0;

	  yfirst = gridInqYval(gridID, 0);
	  ylast  = gridInqYval(gridID, ysize-1);
	  yinc   = gridInqYinc(gridID);

	  if ( gridtype == GRID_GAUSSIAN_REDUCED )
	    fprintf(stdout, "size : dim = %d  nlat = %d", gridsize, ysize);
	  else
	    fprintf(stdout, "size      : dim = %d  nx = %d  ny = %d", gridsize, xsize, ysize);

	  if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
	    fprintf(stdout, "  np = %d", gridInqNP(gridID));

	  fprintf(stdout, "\n");

          if ( gridInqXvals(gridID, NULL) == 0 ) lxcoord = 0;
          if ( gridInqYvals(gridID, NULL) == 0 ) lycoord = 0;

	  if ( xsize > 0 && lxcoord )
	    {
	      if ( gridtype == GRID_GAUSSIAN_REDUCED )
		{
		  fprintf(stdout, "size : dim = %d  nlat = %d\n", gridsize, ysize);
		  fprintf(stdout, "%*s", nbyte0, "");
		  fprintf(stdout, "longitude : reduced\n");
		}
	      else
		{
		  xfirst = gridInqXval(gridID, 0);
		  xlast  = gridInqXval(gridID, xsize-1);
		  xinc   = gridInqXinc(gridID);
		  fprintf(stdout, "%*s", nbyte0, "");
		  fprintf(stdout, "%-9s : first = %.9g", xname, xfirst);
		  if ( xsize > 1 ) fprintf(stdout, "  last = %.9g", xlast);
		  if ( IS_NOT_EQUAL(xinc, 0) )
		    fprintf(stdout, "  inc = %.9g", xinc);
		  fprintf(stdout, "  %s", xunits);
		  if ( gridIsCircular(gridID) ) fprintf(stdout, "  circular");
		  fprintf(stdout, "\n");
		}
	    }

	  if ( ysize > 0 && lycoord )
	    {
	      fprintf(stdout, "%*s", nbyte0, "");
	      fprintf(stdout, "%-9s : first = %.9g", yname, yfirst);
	      if ( ysize > 1 ) fprintf(stdout, "  last = %.9g", ylast);
	      if ( IS_NOT_EQUAL(yinc, 0) && gridtype != GRID_GAUSSIAN && gridtype != GRID_GAUSSIAN_REDUCED )
		fprintf(stdout, "  inc = %.9g", yinc);
	      fprintf(stdout, "  %s", yunits);
	      fprintf(stdout, "\n");
	    }

	  if ( gridIsRotated(gridID) )
	    {
	      double lonpole, latpole, angle;
	      lonpole = gridInqXpole(gridID);
	      latpole = gridInqYpole(gridID);
	      angle   = gridInqAngle(gridID);
	      fprintf(stdout, "%*s", nbyte0, "");
	      fprintf(stdout, "northpole : lon = %.9g  lat = %.9g", lonpole, latpole);
	      if ( angle > 0 ) fprintf(stdout, "  angle = %.9g", angle);
	      fprintf(stdout, "\n");
	    }

	  if ( gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL) )
	    {
	      fprintf(stdout, "%*s", nbyte0, "");
	      fprintf(stdout, "available :");
	      if ( gridInqXbounds(gridID, NULL) ) fprintf(stdout, " xbounds");
	      if ( gridInqYbounds(gridID, NULL) ) fprintf(stdout, " ybounds");
	      if ( gridHasArea(gridID) )          fprintf(stdout, " area");
	      if ( gridInqMask(gridID, NULL) )    fprintf(stdout, " mask");
	      fprintf(stdout, "\n");
	    }

	  if ( gridtype == GRID_LAEA )
	    {
	      double a, lon_0, lat_0;
	      gridInqLaea(gridID, &a, &lon_0, &lat_0);
	      fprintf(stdout, "%*s", nbyte0, "");
	      fprintf(stdout, "projpar   : a = %g  lon_0 = %g  lat_0 = %g\n", a, lon_0, lat_0);
	    }

	  if ( gridtype == GRID_LCC2 )
	    {
	      double a, lon_0, lat_0, lat_1, lat_2;
	      gridInqLcc2(gridID, &a, &lon_0, &lat_0, &lat_1, &lat_2);
	      fprintf(stdout, "%*s", nbyte0, "");
	      fprintf(stdout, "projpar   : a = %7.0f  lon_0 = %g  lat_0 = %g  lat_1 = %g  lat_2 = %g\n",
		      a, lon_0, lat_0, lat_1, lat_2);
	    }
	}
      else if ( gridtype == GRID_SPECTRAL )
	{
	  fprintf(stdout, "size      : dim = %d  nsp = %d  truncation = %d\n", gridsize, gridsize/2, trunc);
	  fprintf(stdout, "%*s", nbyte0, "");
	  fprintf(stdout, "            complexPacking = %d\n", gridInqComplexPacking(gridID));
	}
      else if ( gridtype == GRID_FOURIER )
	{
	  fprintf(stdout, "size      : dim = %d  nfc = %d  truncation = %d\n", gridsize, gridsize/2, trunc);
	}
      else if ( gridtype == GRID_GME )
	{
	  int ni, nd;
	  ni = gridInqGMEni(gridID);
	  nd = gridInqGMEnd(gridID);
	  fprintf(stdout, "size      : dim = %d  nd = %d  ni = %d\n", gridsize, nd, ni);
	}
      else if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED )
	{
	  if ( gridtype == GRID_CURVILINEAR )
	    fprintf(stdout, "size      : dim = %d  nx = %d  ny = %d", gridsize, xsize, ysize);
	  else
	    fprintf(stdout, "size      : dim = %d", gridsize);

          if ( gridtype == GRID_UNSTRUCTURED && gridInqNvertex(gridID) > 0 )
	    fprintf(stdout, "  nvertex = %d", gridInqNvertex(gridID));

          fprintf(stdout, "\n");

          if ( gridtype == GRID_UNSTRUCTURED )
            {
              int number   = gridInqNumber(gridID);
              int position = gridInqPosition(gridID);

              if ( number > 0 )
                {
                  fprintf(stdout, "%*s", nbyte0, "");
                  fprintf(stdout, "grid      : number = %d  position = %d\n", number, position);
                }

              if ( gridInqReference(gridID, NULL) )
                {
                  char reference_link[8192];
                  gridInqReference(gridID, reference_link);
                  fprintf(stdout, "%*s", nbyte0, "");
                  fprintf(stdout, "uri       : %s\n", reference_link);
                }
            }

	  if ( gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL) )
	    {
	      int i;
	      double *xvals, *yvals;
	      double xfirst, xlast, yfirst, ylast;
	      xvals = (double *) malloc(gridsize*sizeof(double));
	      yvals = (double *) malloc(gridsize*sizeof(double));

	      gridInqXvals(gridID, xvals);
	      gridInqYvals(gridID, yvals);

	      xfirst = xvals[0];
	      xlast  = xvals[0];
	      yfirst = yvals[0];
	      ylast  = yvals[0];
	      for ( i = 1; i < gridsize; i++ )
		{
		  if ( xvals[i] < xfirst ) xfirst = xvals[i];
		  if ( xvals[i] > xlast )  xlast  = xvals[i];
		  if ( yvals[i] < yfirst ) yfirst = yvals[i];
		  if ( yvals[i] > ylast )  ylast  = yvals[i];
		}

	      fprintf(stdout, "%*s", nbyte0, "");
	      fprintf(stdout, "%-9s : min = %.9g  max = %.9g  %s", xname, xfirst, xlast, xunits);
	      if ( gridIsCircular(gridID) ) fprintf(stdout, "  circular");
	      fprintf(stdout, "\n");
	      fprintf(stdout, "%*s", nbyte0, "");
	      fprintf(stdout, "%-9s : min = %.9g  max = %.9g  %s\n", yname, yfirst, ylast, yunits);

	      free(xvals);
	      free(yvals);
	    }
	}
      else if ( gridtype == GRID_LCC )
	{
	  double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
	  int projflag, scanflag;

	  gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
		     &projflag, &scanflag);

	  fprintf(stdout, "size      : dim = %d  nx = %d  ny = %d  ", gridsize, xsize, ysize);
	  if ( (projflag&128) == 0 )
	    fprintf(stdout, "North Pole\n");
	  else
	    fprintf(stdout, "South Pole\n");
	  fprintf(stdout, "%*s", nbyte0, "");
	  fprintf(stdout, "            originLon = %g  originLat = %g  lonParY = %g\n",
		  originLon, originLat, lonParY);
	  fprintf(stdout, "%*s", nbyte0, "");
	  fprintf(stdout, "            lat1 = %g  lat2 = %g  xinc = %g m  yinc = %g m\n",
		  lat1, lat2, xincm, yincm);
	}
      else /* if ( gridtype == GRID_GENERIC ) */
	{
	  if ( ysize == 0 )
	    fprintf(stdout, "size      : dim = %d\n", gridsize);
	  else
            fprintf(stdout, "size      : dim = %d  nx = %d  ny = %d\n", gridsize, xsize, ysize);
	}

      if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED || gridtype == GRID_LCC )
	{
	  if ( gridInqXvals(gridID, NULL) || gridInqYvals(gridID, NULL) || gridHasArea(gridID) ||
	       gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL) )
	    {
	      fprintf(stdout, "%*s", nbyte0, "");
	      fprintf(stdout, "available :");
	      if ( gridInqXvals(gridID, NULL) )   fprintf(stdout, " xvals");
	      if ( gridInqYvals(gridID, NULL) )   fprintf(stdout, " yvals");
	      if ( gridInqXbounds(gridID, NULL) ) fprintf(stdout, " xbounds");
	      if ( gridInqYbounds(gridID, NULL) ) fprintf(stdout, " ybounds");
	      if ( gridHasArea(gridID) )          fprintf(stdout, " area");
	      if ( gridInqMask(gridID, NULL) )    fprintf(stdout, " mask");
	      fprintf(stdout, "\n");
	    }
	}

      gridInqUUID(gridID, uuidOfHGrid);
      if ( uuidOfHGrid[0] != 0 )
        {
          char uuidOfHGridStr[37];
          uuid2str(uuidOfHGrid, uuidOfHGridStr);
          if ( uuidOfHGridStr[0] != 0  && strlen(uuidOfHGridStr) == 36 )
            {
	      fprintf(stdout, "%*s", nbyte0, "");
	      fprintf(stdout, "uuid      : %s\n", uuidOfHGridStr);
            }
        }
    }
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
