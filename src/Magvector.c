#if  defined  (HAVE_CONFIG_H)
#  include "config.h" /* HAVE_LIBMAGICS */
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"

#if  defined  (HAVE_LIBMAGICS)
#include "magics_api.h"
#endif


#if  defined  (HAVE_LIBXML)

#include<libxml/parser.h>
#include<libxml/tree.h>
#include "template_parser.h"
#include "magics_template_parser.h"
#include "results_template_parser.h"

//xmlDoc *param_doc = NULL;
//extern xmlNode *root_node, *magics_node, *results_node;
extern xmlNode  *magics_node;

#endif

int VECTOR, STREAM;

static
void magvector( const char *plotfile, int operatorID, const char *varname, long nlon, long nlat, double *grid_center_lon, double *grid_center_lat, double *uarray, double *varray )

{
        long i;
        double dlon = 0, dlat = 0;
        double thin_fac;
	char plotfilename[4096];

	if( uarray == NULL && varray == NULL )
	{
          	fprintf( stderr," No Velocity Components in input file, cannot creaate Vector PLOT!\n" );
		return ;
	}

	if( uarray == NULL || varray == NULL )
	{
          	fprintf( stderr," Found only one Velocity Component in input file, cannot create Vector PLOT!\n" );
		return ;
	}

        if ( nlon > 1 )
        {
           for ( i = 1; i < nlon; ++i ) dlon += (grid_center_lon[i] - grid_center_lon[i-1]);
               dlon /= (nlon-1);
        }  

        if ( nlat > 1 )
        {
           for ( i = 1; i < nlat; ++i ) dlat += (grid_center_lat[nlon*i] - grid_center_lat[nlon*(i-1)]);
               dlat /= (nlat-1);
        }


#if  defined  (HAVE_LIBMAGICS)
        
        magics_template_parser( magics_node );

        /* results_template_parser(results_node, varname ); */

        sprintf(plotfilename, "%s", plotfile);

        mag_setc ("output_name",      plotfilename);
	/* Set the input data */

        mag_setr("input_field_initial_latitude", grid_center_lat[0]);
        mag_setr("input_field_latitude_step", dlat);

        mag_setr("input_field_initial_longitude", grid_center_lon[0]);
        mag_setr("input_field_longitude_step", dlon);

	mag_set2r("input_wind_u_component", uarray, nlon, nlat);
	mag_set2r("input_wind_v_component", varray, nlon, nlat);

        if ( operatorID == VECTOR ) 
	{
		/* Magics functions for performing vector operation */
		/*
		
		mag_setc("wind_legend_only", "on" );
		mag_setc("wind_legend_text", "on" );
		*/
		mag_setc( "legend", "on" );
		mag_setc( "wind_flag_cross_boundary", "on" );
		//mag_setr( "wind_arrow_unit_velocity", 1.0);
		mag_seti( "wind_arrow_thickness",1 );
		mag_coast();
		mag_text();
                mag_enqr("wind_thinning_factor",&thin_fac);
                printf( " %g \n", thin_fac );
		mag_wind();
	}
#else
        cdoAbort("MAGICS support not compiled in!");

#endif

}


#if  defined  (HAVE_LIBMAGICS)

static
void init_MAGICS( )

{
	mag_open();
}

static
void quit_MAGICS( )

{

  mag_close ();
  fprintf( stdout,"Exiting From MAGICS\n" );

}

#endif


void *Magvector(void *argument)

{
  int operatorID;
  int varID, recID;
  int gridsize;
  int gridID;
  int nrecs;
  int levelID;
  int tsID;
  int streamID;
  int vlistID;
  int nmiss;
  int nlon, nlat;
  int nlev;
  int zaxisID, taxisID;
  int vdate, vtime;
  int found;
  char varname[CDI_MAX_NAME];
  double missval;
  double *uarray = NULL;
  double *varray = NULL;
  double *grid_center_lat = NULL, *grid_center_lon = NULL;
  char units[CDI_MAX_NAME];
  char vdatestr[32], vtimestr[32];

  char  *Filename = "combined.xml";

  cdoInitialize(argument);


  VECTOR  = cdoOperatorAdd("vector", 0, 0, NULL);
  STREAM  = cdoOperatorAdd("stream", 0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID = streamOpenRead(cdoStreamName(0));

  vlistID = streamInqVlist(streamID);
  taxisID = vlistInqTaxis(vlistID);

  found = 0;
  varID = 0;
  gridID  = vlistInqVarGrid(vlistID, varID);
  zaxisID = vlistInqVarZaxis(vlistID, varID);
  missval = vlistInqVarMissval(vlistID, varID);

  if ( gridInqType(gridID) == GRID_GME          ) cdoAbort("GME grid unspported!");
  if ( gridInqType(gridID) == GRID_UNSTRUCTURED ) cdoAbort("Unstructured grid unspported!");

  if ( gridInqType(gridID) != GRID_CURVILINEAR )
    gridID = gridToCurvilinear(gridID, 1);

  gridsize = gridInqSize(gridID);
  nlon     = gridInqXsize(gridID);
  nlat     = gridInqYsize(gridID);
  nlev     = zaxisInqSize(zaxisID);

  uarray           = (double *) malloc(gridsize*sizeof(double));
  varray           = (double *) malloc(gridsize*sizeof(double));
  grid_center_lat = (double *) malloc(gridsize*sizeof(double));
  grid_center_lon = (double *) malloc(gridsize*sizeof(double));

  gridInqYvals(gridID, grid_center_lat);
  gridInqXvals(gridID, grid_center_lon);

  /* Convert lat/lon units if required */
  gridInqXunits(gridID, units);
  gridToDegree(units, "grid center lon", gridsize, grid_center_lon);
  gridInqYunits(gridID, units);
  gridToDegree(units, "grid center lat", gridsize, grid_center_lat);
					
  tsID = 0;

#if  defined  (HAVE_LIBXML)
  /* HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR */
  init_XMLtemplate_parser( Filename );
  updatemagics_and_results_nodes( );
#endif


#if  defined  (HAVE_LIBMAGICS)
  init_MAGICS( );
#endif

  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);
	      
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID, &varID, &levelID);

	  vlistInqVarName(vlistID, varID, varname);

          if ( operatorID == VECTOR )
	  {
	  	if( !strcmp( varname, "var131" ) || !strcmp( varname, "u" ) ) /* U Velocity as per GRIB is 'var131, as per NC 'u' */
	  	{
          		fprintf( stderr,"Found U VEL in Varname %s\n",varname );
			streamReadRecord(streamID, uarray, &nmiss);
			found ++;
	  	}

	  	if( !strcmp( varname, "var132" ) || !strcmp( varname, "v" ) ) /* V Velocity as per GRIB  is 'var132, as per NC 'v'*/
	  	{
          		fprintf( stderr,"Found V VEL in Varname %s\n",varname );
			streamReadRecord(streamID, varray, &nmiss);
			found ++;
	  	}	

	  	if( found == 2 )
	    		break;
	  }
	  else if ( operatorID == STREAM )
          	fprintf( stderr," Stream Operator Un-Supported!\n" );
	  else 
          	fprintf( stderr," Operator Un-Supported!\n" );
		
        }
         
        if ( operatorID == VECTOR )
	{
	  	if( found == 2 )
	  	{
          		fprintf( stderr,"Found Both U & V VEL, Creating vector fields! \n" );
			magvector(cdoStreamName(1), operatorID, varname, nlon, nlat, grid_center_lon, grid_center_lat, uarray, varray );
	  	}
	  	else if( found == 1 )
          		fprintf( stderr,"Found only one Velocity Component in input file, cannot creaate Vector PLOT!\n" );
	  	else if( found == 0 )
          		fprintf( stderr,"No Velocity Components in input file, cannot create Vector PLOT!\n" );
	}
        break;

        tsID++;
    }

  streamClose(streamID);

  if ( uarray  ) free(uarray);
  if ( varray  ) free(varray);
  if ( grid_center_lon ) free(grid_center_lon);
  if ( grid_center_lat ) free(grid_center_lat);

#if  defined  (HAVE_LIBXML)
  quit_XMLtemplate_parser( );
#endif

#if  defined  (HAVE_LIBMAGICS)
  quit_MAGICS( );
#endif

  cdoFinish();

  return (0);

}
