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

xmlDoc *param_doc = NULL;
xmlNode *root_node = NULL, *magics_node = NULL, *results_node = NULL;

#endif

#define DBG 0


int CONTOUR, SHADED, GRFILL;

static
void magplot( const char *plotfile, int operatorID, const char *varname, long nlon, long nlat, double *grid_center_lon, double *grid_center_lat, double *array )

{
  static int once = 1;
  long i;
  double dlon = 0, dlat = 0;
  char plotfilename[4096];
  char *titlename;

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

  sprintf(plotfilename, "%s_%s", plotfile, varname);

  titlename = strdup( plotfilename );

#if  defined  (HAVE_LIBMAGICS)

  mag_setc ("output_name",      plotfilename);

  // Set the input data arrays to magics++
   
  mag_set2r("input_field", array, nlon, nlat);

  /*
  	mag_setc("input_field_organization", "REGULAR");
  	mag_set2r("input_field_latitudes", grid_center_lat, nlon, nlat);
  	mag_set2r("input_field_longitudes", grid_center_lon, nlon, nlat);
  */

  mag_setr("input_field_initial_latitude", grid_center_lat[0]);
  mag_setr("input_field_latitude_step", dlat);

  mag_setr("input_field_initial_longitude", grid_center_lon[0]);
  mag_setr("input_field_longitude_step", dlon);
 

#if 0
  if( once )
  {
	  magics_template_parser( magics_node );
	  once = 0;
  }
#endif

  magics_template_parser( magics_node );
  results_template_parser(results_node, varname );


  /* set up the coastline attributes */
  /* mag_setc ("map_coastline_colour", "khaki"); */
  /* mag_setc ("map_grid_colour",      "grey");  */ 

  /* define the contouring parameters */
  if ( operatorID == SHADED )
    {

      mag_setc ( "contour", "off" );
      mag_setc ( "contour_shade", "on" );
      mag_setc ( "contour_shade_method", "area_fill" );
      mag_setc ( "contour_label", "off" );

      /* Adjust Set The page slightly to fit the legend */
      mag_setr ( "subpage_x_length", 24. );
      mag_setr ( "subpage_y_length", 30. );

      /* Legend Settings */
      mag_setc ( "legend", "on" );
      mag_setc ( "legend_display_type", "continuous" );
      mag_setc ( "legend_entry_plot_direction", "column" );
      mag_setc ( "legend_box_mode", "positional" );
      mag_setr ( "legend_box_x_position", 26.5 );
      mag_setr ( "legend_box_y_position", 0.39 );
      mag_setr ( "legend_box_x_length", 2.0 );
      mag_setr ( "legend_box_y_length", 12.69 );

    }
  else if ( operatorID == CONTOUR )
    {

      mag_setc ("contour",                  "on");
      // mag_setc ("contour_line_colour",      "sky");
      // mag_setc ("CONTOUR_HIGHLIGHT_COLOUR", "GREEN");
      mag_setc ("contour_shade",            "off");
      mag_setc ("contour_label",            "on");

    }
  else if ( operatorID == GRFILL )
    {

      mag_setc ( "contour", "off" );
      mag_setc ( "contour_shade", "on" );

      mag_setc ( "contour_shade_technique", "cell_shading" );

      mag_setc ( "contour_shade_method", "area_fill" );
      mag_setc ( "contour_label", "off" );

      /* Adjust Set The page slightly to fit the legend */
      mag_setr ( "subpage_x_length", 24. );
      mag_setr ( "subpage_y_length", 30. );

      /* Legend Settings */
      mag_setc ( "legend", "on" );
      mag_setc ( "legend_display_type", "continuous" );
      mag_setc ( "legend_entry_plot_direction", "column" );
      mag_setc ( "legend_box_mode", "positional" );
      mag_setr ( "legend_box_x_position", 26.5 );
      mag_setr ( "legend_box_y_position", 0.39 );
      mag_setr ( "legend_box_x_length", 2.0 );
      mag_setr ( "legend_box_y_length", 12.69 );

      fprintf( stderr, " GrFILL Done!\n");
    }

  /* plot the title text and the coastlines */
  mag_cont ();
  mag_coast ();


  mag_set1c("text_lines", (const char **) &titlename, 1);
  mag_setc("text_colour", "black");

/*
  mag_setr("text_font_size", 0.6);
  mag_setc("text_mode", "positional");
  mag_setr("text_box_x_position", 1.5);
  mag_setr("text_box_y_position", 16.5);
  mag_setr("text_box_x_length", 20.);
  mag_setr("text_box_y_length", 2.5);
  mag_setc("text_border", "off");
*/

  mag_setc("text_justification", "left");
  mag_text();

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
  if( DBG )
    fprintf( stderr,"Exiting From MAGICS\n" );

}

#endif


void *Magplot(void *argument)
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
  char varname[CDI_MAX_NAME];
  double missval;
  double *array = NULL;
  double *grid_center_lat = NULL, *grid_center_lon = NULL;
  char units[CDI_MAX_NAME];
  char vdatestr[32], vtimestr[32];

  char  *Filename = "combined.xml";

  cdoInitialize(argument);

  CONTOUR = cdoOperatorAdd("contour", 0, 0, NULL);
  SHADED  = cdoOperatorAdd("shaded", 0, 0, NULL);
  GRFILL  = cdoOperatorAdd("grfill", 0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID = streamOpenRead(cdoStreamName(0));

  vlistID = streamInqVlist(streamID);
  taxisID = vlistInqTaxis(vlistID);

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

  array           = (double *) malloc(gridsize*sizeof(double));
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



  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);
	      
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));

      for ( recID = 0; recID < nrecs; recID++ )
	{

#if  defined  (HAVE_LIBMAGICS)
	  init_MAGICS( );
#endif
	  streamInqRecord(streamID, &varID, &levelID);
	  streamReadRecord(streamID, array, &nmiss);

	  vlistInqVarName(vlistID, varID, varname);


	  if ( operatorID == SHADED || operatorID == CONTOUR || operatorID == GRFILL )
          {

                if( operatorID == SHADED )
                   fprintf( stderr," Creating SHADED PLOT for %s\n",varname );
                else if( operatorID == CONTOUR )
                   fprintf( stderr," Creating CONTOUR PLOT for %s\n",varname );
                else if( operatorID == GRFILL )
                   fprintf( stderr," Creating GRFILL PLOT for %s\n",varname );

	  	magplot(cdoStreamName(1), operatorID, varname, nlon, nlat, grid_center_lon, grid_center_lat, array);

          }
	  else
	  	fprintf(stderr,"operator not implemented\n");

	  //	  break;

#if  defined  (HAVE_LIBMAGICS)
	  quit_MAGICS( );
#endif
	}

      break;

      tsID++;
    }

  streamClose(streamID);

  if ( array  ) free(array);
  if ( grid_center_lon ) free(grid_center_lon);
  if ( grid_center_lat ) free(grid_center_lat);

#if  defined  (HAVE_LIBXML)
  quit_XMLtemplate_parser( );
#endif

  cdoFinish();

  return (0);

}
