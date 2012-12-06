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
#include <ctype.h>


extern xmlNode  *magics_node;

#endif

#define DBG 0


char *line_colours[] = {

     "RGB(1.,0.,0.)",
     "RGB(0.,1.,0.)",
     "RGB(0.,0.,1.)",
     "RGB(0.,1.,1.)",
     "RGB(1.,1.,0.)",
     "RGB(1.,0.,1.)",
     "RGB(0.,0.,0.)",

};

char  *graph_params[] = {"ymin","ymax","sigma","stat","obsv","xml"};

int graph_param_count = sizeof(graph_params)/sizeof(char*);
int num_colours = sizeof( line_colours )/sizeof( char* );

void VerifyGraphParameters( int num_param, char **param_names );

extern int IsNumeric();
extern void StrToUpperCase();
extern int StringSplitWithSeperator();


static
void maggraph(const char *plotfile, const char *varname,const char *varunits, long nfiles, long nts, int *vdate, int *vtime, double **datatab, int nparam, char **params)
{
  
  char *lines[1];
  char *temp_str;
  char **split_str = NULL;
  char *sep_char = "=";
  char *date_time_str[nts];
  char vdatestr[32], vtimestr[32], legend_text_data[256];
  int num_sigma = 2;
  int stat = FALSE, obsv = FALSE;
  int split_str_count;
  int file_begin = 0;
  int count ;
  int num_years = 0, num_months = 0, num_days = 0;
  long tsID, fileID, i, j, k;
  double *date_time;
  double min_val = 1.0e+200, max_val = -1.0e+200;
  double suppress_min_val = 1.0e+200, suppress_max_val = -1.0e+200;
  double *mean_val, *std_dev_val;
  double *spread_min, *spread_max;
  double y_min_val = 1.0e+200, y_max_val = -1.0e+200;

  
  
  if( DBG )
    {
      fprintf(stderr, "Num params %d\n", nparam);
  
      for( i = 0; i< nparam; i++ )
	fprintf(stderr, "Param %s\n", params[i]);
    }
  
  for( i = 0; i < nparam; ++i )
    {
      split_str_count = 0;
      sep_char = "=";
      split_str_count = StringSplitWithSeperator( params[i], sep_char, &split_str );
      
      if( !strcmp( split_str[0],"obsv" ) ) 
	{  
	  temp_str = strdup( split_str[1] );    
	  StrToUpperCase( temp_str );
	  if( !strcmp( temp_str, "TRUE" ) )
	    {
	      obsv = TRUE;
	      file_begin = 1;
	      if( DBG )
		fprintf( stderr,"OBSV TRUE\n" );
	    }
	}
	
      if( !strcmp( split_str[0],"stat" ) ) 
	{  
	  temp_str = strdup( split_str[1] );    
	  StrToUpperCase( temp_str );
	  
	  if( !strcmp( temp_str, "TRUE" ) )
	    {
	      stat = TRUE;
	      if( DBG )
		fprintf(stderr,"STAT TRUE\n");
	    }
	}
	
      if( !strcmp( split_str[0],"ymin" ) )
	{
	  y_min_val = atof( split_str[1] );
	  if( DBG )
	    fprintf(stderr,"Y min Val %g\n",y_min_val);
	}
	
      if( !strcmp( split_str[0],"ymax" ) )
	{
	  y_max_val = atof( split_str[1] );
	  if( DBG )
	    fprintf(stderr,"Y max Val %g\n",y_max_val);
	}
	
      if( !strcmp( split_str[0],"sigma" ) )
	{
	  num_sigma = atof( split_str[1] );
	  if( DBG )
	    fprintf(stderr,"SIGMA %d\n",num_sigma);
	}   
      
      free( split_str );
      /*for( k = 0; k < split_str_count; ++k )
	{
	  free( split_str[k] );
	}
      */	  
    }
    

  date_time = (double *) malloc( nts*sizeof(double) );
  mean_val  = (double *) malloc( nts*sizeof(double) );
  std_dev_val = (double *) malloc( nts*sizeof(double) );
  spread_min = (double *) malloc( nts*sizeof(double) );
  spread_max = (double *) malloc( nts*sizeof(double) );


  if( DBG )
    {
      fprintf(stderr," %6d %6d\n", nfiles, nts );
      fprintf(stderr,"\n");
    }

  for ( tsID = 0; tsID < nts; ++tsID )
    {
      date_time[tsID] = tsID+1;
      date2str(vdate[tsID], vdatestr, sizeof(vdatestr));
      time2str(vtime[tsID], vtimestr, sizeof(vtimestr));
      date_time_str[tsID] = (char *)malloc(256);
      sprintf(date_time_str[tsID], "%s %s", vdatestr, vtimestr);
      mean_val[tsID] = 0.;
      std_dev_val[tsID] = 0.;

      if( DBG )
	{
	  fprintf(stderr,"%d: %s\n", tsID, date_time_str[tsID]);
	  fprintf(stderr,"%6d %6d", vdate[tsID], vtime[tsID]);
	}
	
      for ( fileID = 0; fileID < nfiles; ++fileID )
	{
	  if( DBG )
	    printf("%d\n", fileID );
	  if( datatab[fileID][tsID] < min_val )
	    min_val = datatab[ fileID ][ tsID ];	
	  if( datatab[fileID][tsID] > max_val )
	    max_val = datatab[ fileID ][ tsID ];	
	  
	  mean_val[tsID] += datatab[fileID][tsID];
	  std_dev_val[tsID] = 0.;
	  spread_min[tsID] = 0.;
	  spread_max[tsID] = 0.;

	  if( DBG )
	    {
	      fprintf(stderr," %6g", datatab[fileID][tsID]);
	      fprintf(stderr,"\n");
	    }
	}
    }

  for ( tsID = 0; tsID < nts; ++tsID )
    {
    	mean_val[tsID] /= ( double )nfiles;
        spread_min[tsID] = mean_val[tsID];
        spread_max[tsID] = mean_val[tsID];

        for ( fileID = 0; fileID < nfiles; ++fileID )
          {
             std_dev_val[tsID] += ( datatab[fileID][tsID]-mean_val[tsID] ) * ( datatab[fileID][tsID]-mean_val[tsID] );
          }
        std_dev_val[tsID] /= ( double )nfiles;
        std_dev_val[tsID] = pow( std_dev_val[tsID], 0.5 ); 
	
	if( DBG )
	  fprintf(stderr," Mean : %g Std Dev: %g\n",mean_val[tsID],std_dev_val[tsID] ); 

        spread_min[tsID] = mean_val[tsID] - num_sigma * std_dev_val[tsID];
        spread_max[tsID] = mean_val[tsID] + num_sigma * std_dev_val[tsID];
	
	if( DBG )
	  fprintf(stderr," Min : %g Max: %g\n",spread_min[tsID],spread_max[tsID] ); 

    }
  
  for ( tsID = 0; tsID < nts; ++tsID )
    {
         if( spread_min[tsID] < min_val )
            min_val = spread_min[ tsID ];	
         if( spread_max[tsID] > max_val )
            max_val = spread_max[ tsID ];	
    }

  
    if( DBG )
    {
      fprintf(stderr," %6g %6g\n", min_val, max_val );
      fprintf(stderr," %s %s\n", date_time_str[0], date_time_str[ nts-1 ] );
      fprintf(stderr,"\n");
    }


    //fprintf(stderr,"HERE 1\n");
    split_str_count = 0;
    sep_char = "-";
    split_str_count = StringSplitWithSeperator( date_time_str[ nts - 1 ], sep_char, &split_str );
    num_years = atoi( split_str[0] );
    num_months = atoi( split_str[1] );
    num_days   = atoi( split_str[2] );
    free( split_str  );
    /*
    for( k = 0; k < split_str_count; ++k )
      {
	free( split_str[k] );
      }
   */	
    
    
    split_str_count = StringSplitWithSeperator( date_time_str[0], sep_char, &split_str );
    num_years -= atoi( split_str[0] );
    //fprintf(stderr,"HERE 2\n");
    if( num_years <= 1 )
      {
	if( num_years == 1 )
	  num_months += ( 12 - atoi( split_str[1] ) );
	else
	  num_months -= ( atoi( split_str[1] ) );
	
	if( !num_months )
	  num_days -= atoi( split_str[2] );
	else if( num_months == 1 )
	  num_days += ( 31- atoi( split_str[2] ) );
      }
      
    /*  
    for( k = 0; k < split_str_count; ++k )
      {
	free( split_str[k] );
      }
    */  
    free( split_str );
    if( DBG )
      fprintf(stderr," %d %d\n", num_years, num_months );
    

  /* 
	1. Loop over the Files
	2. Loop over the number of time steps 
	3. Set the attributes for the magics data and plot
  */  
   
#if  defined  (HAVE_LIBMAGICS)

  magics_template_parser( magics_node );

  mag_setc("output_name", plotfile);
  mag_setc("subpage_map_projection", "cartesian"); 
  mag_setr("subpage_y_length", 14.);
  mag_setr("subpage_y_position", 1.5);


  /* Horizontal Axis attributes */
  mag_setc("axis_orientation","horizontal");
  mag_setc("axis_grid", "on");
  mag_setc("axis_grid_colour", "grey");
  mag_seti("axis_grid_thickness", 1);
  mag_setc("axis_grid_line_style", "dot");
  mag_setc("axis_type", "date");
  
  if( num_years > 1 )
    mag_setc("axis_date_type", "years");
  else if( num_years <= 1 )
    {
      if( num_months > 1 )
	mag_setc("axis_date_type", "months");
      else
	{
	  if( num_months == 1 )
	    mag_setc("axis_date_type", "days");
	  else
	    {
	      if( num_days )
		mag_setc("axis_date_type", "days");
	      else
		mag_setc("axis_date_type", "hours");
	    }
	}
    }
  
  
  mag_setc("axis_date_min_value", date_time_str[0]);
  mag_setc("axis_date_max_value", date_time_str[nts-1]);
  mag_setc("axis_title_text","Time");
  mag_setc("axis_title_orientation","horizontal");
  mag_axis();

  /* Vertical Axis attributes */
  mag_setc("axis_orientation", "vertical");
  mag_setc("axis_grid", "on");
  mag_setc("axis_type", "regular");
  mag_setc("axis_grid_colour", "grey");
  mag_seti("axis_grid_thickness", 1);
  mag_setc("axis_grid_line_style", "dot");

  /*  To redefine the y- axis scale based on user input in .xml file */

  /*
  mag_enqr("graph_y_suppress_above",&suppress_max_val);
  mag_enqr("graph_y_suppress_below",&suppress_min_val);
  
  if( DBG )
    fprintf(stderr," %6g %6g\n", suppress_min_val, suppress_max_val );
  
  if( min_val < suppress_min_val )
      min_val = suppress_min_val;
  
  if( max_val > suppress_max_val )
      max_val = suppress_max_val;
  */
 
  mag_setr("axis_min_value", min_val);
  mag_setr("axis_max_value", max_val);
  
  if( y_min_val < 1.0e+200 )
    mag_setr("axis_min_value", y_min_val);
  
  if( y_max_val > -1.0e+200)
    mag_setr("axis_max_value", y_max_val);
  
  mag_setc("axis_title_text",varname);
  
  mag_setc("axis_title_orientation","vertical");
  mag_axis();

  

  /* Legend */
  mag_setc("legend", "on");
  mag_setc("legend_text_colour", "black");

  mag_setc("graph_symbol","off");
  mag_seti("graph_line_thickness", 8 );
  
  for ( i = file_begin; i < nfiles; ++i )
    {
      count = i; 
      if( obsv == TRUE )
	count = i -1;
      
      sprintf(legend_text_data, "ens_%d", count + 1);
      mag_setc("graph_line_colour", line_colours[ count%num_colours ]);
      mag_setc("legend_user_text", legend_text_data);
      mag_set1c("graph_curve_date_x_values",(const char**)date_time_str, nts);
      mag_set1r("graph_curve_y_values", datatab[i], nts);
      mag_graph ();
    }
  
    
  if( obsv == TRUE )
    {
      mag_setc("graph_line_colour", line_colours[ num_colours - 1 ]);
      sprintf(legend_text_data, "%s","Obsv" );
      mag_setc("legend_user_text", legend_text_data);
      mag_set1c("graph_curve_date_x_values",(const char**)date_time_str, nts);
      mag_set1r("graph_curve_y_values", datatab[0], nts);
      mag_setc("graph_line_style", "dot" );
      mag_seti("graph_line_thickness", 10 );
      mag_graph ();
    }
    

  if( stat == TRUE )
    {
      mag_seti("graph_line_thickness", 8 );
      mag_setc("graph_line_colour", "grey" );
      mag_setc("graph_line_style", "dash" );
      mag_set1c("graph_curve_date_x_values", (const char**)date_time_str, nts);
      mag_set1r("graph_curve_y_values",mean_val, nts);
      sprintf(legend_text_data, "Mean");
      mag_setc("legend_user_text", legend_text_data);
      mag_graph ();

      mag_reset("graph_type");
      mag_setc("graph_type", "area");
      mag_seti("graph_line_thickness", 1 );
      mag_setc("graph_shade_style", "dot");
      mag_setr("graph_shade_dot_size",1.);
      mag_set1c("graph_curve2_date_x_values", (const char**)date_time_str, nts);
      mag_set1r("graph_curve2_y_values",spread_max, nts);
      mag_set1c("graph_curve_date_x_values", (const char**)date_time_str, nts);
      mag_set1r("graph_curve_y_values",spread_min, nts);
      mag_setc("graph_shade_colour", "grey");
      sprintf(legend_text_data, "%dSigma", num_sigma);
      mag_setc("legend_user_text", legend_text_data);
      mag_graph ();
    }
  
  
  lines[0] = (char *)malloc(1024);
  sprintf(lines[0],"%s","ExpID : ");/* To be obtained from Meta Data */
  sprintf( lines[0],"%sxxxx  Variable : %s[%s]",lines[0], varname, varunits );
  sprintf( lines[0],"%s  Date : %s --%s",lines[0], date_time_str[0], date_time_str[ nts-1 ] );
  mag_set1c("text_lines", (const char**)lines, 1);
  
  if( DBG )
    fprintf(stderr, "%s\n",lines[0]);
  
  mag_setc("text_html", "true");
  mag_setc("text_colour", "black");
  mag_setr("text_font_size", 0.6);
  mag_setc("text_mode", "positional");
  mag_setr("text_box_x_position", 1.5);
  mag_setr("text_box_y_position", 16.5);
  mag_setr("text_box_x_length", 20.);
  mag_setr("text_box_y_length", 2.5);
  mag_setc("text_border", "off");
  mag_setc("text_justification", "left");
  mag_text();

  free( date_time );
  free( mean_val );
  free( std_dev_val );
  free( spread_min );
  free( spread_max );

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
    fprintf( stdout,"Exiting From MAGICS\n" );

}

#endif

#define NINC_ALLOC 1024

void *Maggraph(void *argument)
{
  const char *ofilename;
  char varname[CDI_MAX_NAME], units[CDI_MAX_NAME];
  char  *Filename = "combined.xml";
  /* char  *Filename = "combined.xml"; */
  char **split_str = NULL;
  char sep_char = '=';
  char *temp_str;
  char **pnames = NULL;
  int operatorID;
  int varID, levelID, recID;
  int gridID;
  int nrecs;
  int tsID;
  int streamID;
  int vlistID, vlistID0 = -1;
  int nmiss;
  int zaxisID, taxisID;
  int *vdate = NULL, *vtime = NULL;
  int fileID, nfiles;
  int nts = 0, nts_alloc = 0;
  int nparam = 0;
  int  found = FALSE, syntax = TRUE, halt_flag = FALSE, split_str_count;
  double missval;
  double **datatab = NULL;
  double val;
  int i;
  
  cdoInitialize(argument);

  nparam = operatorArgc();
  pnames = operatorArgv();
  
  if( nparam )
    VerifyGraphParameters(nparam,pnames);
  
  nfiles = cdoStreamCnt() - 1;
  ofilename = cdoStreamName(nfiles);

  datatab = (double **) malloc(nfiles*sizeof(double *));

  for ( fileID = 0; fileID < nfiles; fileID++ )
    datatab[fileID] = NULL;

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID = streamOpenRead(cdoStreamName(fileID));

      vlistID = streamInqVlist(streamID);
      taxisID = vlistInqTaxis(vlistID);

      vlistInqVarUnits(vlistID, 0, units);
      if( DBG )
	fprintf(stderr," %s\n", units );
      if ( fileID == 0 )
	{
	  vlistInqVarName(vlistID, 0, varname);
	  
	  
	  gridID = vlistInqVarGrid(vlistID, 0);

	  if ( gridInqSize(gridID) != 1 ) cdoAbort("Variable has more than one grid point!");

	  vlistID0 = vlistDuplicate(vlistID);
	}
      else
	{
	  vlistCompare(vlistID0, vlistID, CMP_ALL);
	}

      tsID = 0;
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{
	  if ( nrecs != 1 ) cdoAbort("Input streams have more than one record!\n");
	  if ( fileID == 0 )
	    {
	      nts++;

	      if ( nts > nts_alloc )
		{
		  nts_alloc += NINC_ALLOC;
		  datatab[fileID] = (double *) realloc(datatab[fileID], nts_alloc*sizeof(double));
		  vdate = (int *) realloc(vdate, nts_alloc*sizeof(int));
		  vtime = (int *) realloc(vtime, nts_alloc*sizeof(int));
		}
	      vdate[tsID] = taxisInqVdate(taxisID);
	      vtime[tsID] = taxisInqVtime(taxisID);
	    }
	  else
	    {
	      if ( (tsID+1) > nts ) cdoAbort("Too many timesteps in stream %s", cdoStreamName(fileID));

	      if ( tsID == 0 )
		{
		  datatab[fileID] = (double *) malloc(nts*sizeof(double));
		}
	    }
	  
	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID, &varID, &levelID);
	      streamReadRecord(streamID, &val, &nmiss);	
	      datatab[fileID][tsID] = val;
	    }

	  tsID++;
	}

      streamClose(streamID);
    }
  
#if  defined  (HAVE_LIBXML)
  /* HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR */
  init_XMLtemplate_parser( Filename );
  updatemagics_and_results_nodes( );
#endif


#if  defined  (HAVE_LIBMAGICS)
  init_MAGICS( );
#endif

  cdoPrint(" Creating PLOT for %s", varname);
  if( DBG )
    {
      fprintf(stderr, "Num params %d\n", nparam);
  
      for( i = 0; i< nparam; i++ )
	fprintf(stderr, "Param %s\n", pnames[i]);
    }
  maggraph(ofilename, varname, units, nfiles, nts, vdate, vtime, datatab, nparam, pnames);

#if  defined  (HAVE_LIBXML)
  quit_XMLtemplate_parser( );
#endif

#if  defined  (HAVE_LIBMAGICS)
  quit_MAGICS( );
#endif

  if ( vlistID0 != -1 ) vlistDestroy(vlistID0);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      if ( datatab[fileID] ) free(datatab[fileID]);
    }

  free(datatab);


  if ( vdate ) free(vdate);
  if ( vtime ) free(vtime);

  cdoFinish();

  return (0);
}


void VerifyGraphParameters( int num_param, char **param_names )

{
  int i, j, k;
  int  found = FALSE, syntax = TRUE, halt_flag = FALSE, file_found = TRUE, split_str_count;
  char **split_str = NULL;
  char *sep_char = "=";
  char *temp_str;
  FILE *fp;
  
  
  
  for ( i = 0; i < num_param; ++i )
    {
      split_str_count = 0;
      found = FALSE;
      syntax = TRUE;
      split_str_count = StringSplitWithSeperator( param_names[i], sep_char, &split_str );
      if( split_str_count > 1 ) 
	{
	  for ( j = 0; j < graph_param_count; ++j )
	    {
	      if( !strcmp( split_str[0], graph_params[j] ) )
		{
		  found = TRUE;
		  if( !strcmp( split_str[0],"obsv" ) ||  !strcmp( split_str[0],"stat" ) )
		    {  
		      if( IsNumeric( split_str[1] ) )
			syntax = FALSE;
		      else 
			{			
			  temp_str = strdup( split_str[1] );    
			  StrToUpperCase( temp_str );
			  if( strcmp( temp_str,"TRUE" ) && strcmp( temp_str,"FALSE" ) )
			    syntax = FALSE;			      
			}
		    }	 
		      
		  if( !strcmp( split_str[0],"ymin" ) ||  !strcmp( split_str[0],"ymax" ) || !strcmp( split_str[0],"sigma" )  )
		    {
		      if( !IsNumeric( split_str[1] ) )
			syntax = FALSE;       
		    }
		    
		  if( !strcmp( split_str[0],"xml" ) )
		    {
		      if( ( fp = fopen( split_str[1],"r") ) == NULL )
			{
			  fprintf( stderr,"Input XML File not found in specified path '%s'\n", split_str[1] );
			  halt_flag = TRUE;
			}
		      else
			{
#if  defined  (HAVE_LIBXML)
			  /* HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR */
			  fclose(fp);
			  init_XMLtemplate_parser( split_str[1] );
			  updatemagics_and_results_nodes( );
#endif			
			}
		    }
		}
	    }
	}
      else
	{
	  syntax = FALSE;
	}
	
      if( found == FALSE )
	{
	  halt_flag = TRUE;
	  fprintf( stderr,"Invalid parameter  '%s'\n", param_names[i] );
	} 
      if( found == TRUE && syntax == FALSE )
	{
	  halt_flag = TRUE;
	  fprintf( stderr,"Invalid parameter specification  '%s'\n", param_names[i] );
	}
	free( split_str );
      
	/*  
	  for ( k = 0; k < split_str_count; ++k )
	    {  
	      free( split_str[k] );
	    }
	*/    
    }
      
    if( halt_flag == TRUE )
    {
      exit(0);
    }
    
}
