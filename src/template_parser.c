#include "template_parser.h"
#include "magics_template_parser.h"
#include "results_template_parser.h"

#define DBG_MSG 0 


#if  defined  (HAVE_LIBXML)
extern xmlNode *root_node, *magics_node, *results_node;
extern xmlDoc *param_doc;

extern int magics_template_parser();
extern int results_template_parser();
#endif


int template_parser(  char *Filename, const char *varname )

{

#if  defined  (HAVE_LIBXML)
        xmlDoc         *doc = NULL;
        xmlNode        *root_element = NULL;

        doc = xmlReadFile( Filename, NULL, 0 );
        if ( doc == NULL )
        {
                  printf( "Error: Could not parse the file \"%s\"\n", Filename );
        	  return (1);
        }
        else
        {
                  /* 
		     Get the name of the root element node 
		     If "magics" , call "magics" parser
		     If "results", call "results" parser
		  */                      

                  root_element = xmlDocGetRootElement( doc );
		  
		  if( !strcmp( root_element->name, "magics" ) )
		  {
                  	if ( magics_template_parser( root_element ) == 1 )
		  	{
				printf( "Un-Supported version of Magics++! \n" );
	        		return (2);
			}
		  }
		  else if( !strcmp( root_element->name, "results" ) )
		  {
                  	results_template_parser( root_element, varname );
			 /* Needs some error handling */
		  }

                  /*** free the document ***/
                  xmlFreeDoc( doc );
        }

        /*** Free the global variables that may
         *   have been allocated by the parser. 
        ***/

        xmlCleanupParser();

        return (0);
#else
	fprintf(stderr, "XML support not compiled in!");
	return (1);
#endif

}


int init_XMLtemplate_parser( char *Filename )

{

#if  defined  (HAVE_LIBXML)
        param_doc = xmlReadFile( Filename, NULL, 0 );
        if ( param_doc == NULL )
        {
                  printf( "Error: Could not parse the file \"%s\"\n", Filename );
        	  return (1);
        }
        else
        {
		  fprintf( stderr, "XML file %s being parsed \n", Filename );
                  root_node = xmlDocGetRootElement( param_doc );
        }
        return 0;
#else
	fprintf(stderr, "XML support not compiled in!");
	return (1);
#endif

}


int updatemagics_and_results_nodes(  )

{

#if  defined  (HAVE_LIBXML)
    int param_set_flag;
    xmlNode *cur_node = NULL;
	
    if( root_node == NULL )
    {
        printf( "Invalid Root Node\n" );
    	return 0;
    }

    for ( cur_node = root_node->children; cur_node; cur_node = cur_node->next )
    {   
        if ( cur_node->type == XML_ELEMENT_NODE )
        {   
    
#if 0
            fprintf( stdout, "Node Name: %s \n", cur_node->name );
#endif
            if( !strcmp( cur_node->name, "magics" ) ) 
            {
		magics_node = cur_node;
#if 0
                fprintf( stdout, "Node Name: %s \n", cur_node->name );
#endif
	    }  

            if( !strcmp( cur_node->name, "results" ) ) 
            {
		results_node = cur_node;
#if 0
                fprintf( stdout, "Node Name: %s \n", cur_node->name );
#endif
	    }  
	}
    }
    return 0;
#else
    fprintf(stderr, "XML support not compiled in!");
    return (1);
#endif

}


int quit_XMLtemplate_parser( )

{

#if  defined  (HAVE_LIBXML)
        xmlFreeDoc( param_doc );
        xmlCleanupParser( );
	if( param_doc == NULL )
		printf( "Cleaned XML parser\n" );
#if 0
        fprintf( stdout, "Cleaned XML parser\n" );
#endif
	return 0;
#else
	fprintf(stderr, "XML support not compiled in!");
	return (1);
#endif

}
