#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#define ANIM_DIR "./anim"
#define POV_DIR "./pov"


char *program_name;
void usage(void);
int isNumber( char c );

int main(argc, argv)
     int argc;
     char *argv[];
{
  
  int i, j, k; /*loop counters */
  mode_t dir_mode = S_IRWXU;
  
  char pov_flags[120] = ""; /*flags to be placed on the povray command line */
  char *prefix = NULL; /*the prefix for the input and output files */
  int start_frame = 0; /*the start frame */
  int in_frame; /*the current in frame we're working with. */
  int out_frame; /*the current out frame were working with */
  int end_frame; /*the end frame */
  short is_end_frame = 0; /*tells us to stop before we run out of pov files */
  int in_file_exists = 0; /*lets us know if the in_file exists */
  int compression = 0; /* boolean for toggling the compression flag */
  int duplicate = 0; /* boolean to turn on duplicate frames */
  int n_duplicates = 0; /* the number of duplicate frames */
  int done;
  
  char current_flag;

  char command[1000];
  char pov_dir[500];
  char anim_dir[500];
  char out_format[500];
  char in_format[500];
  char in_name[5000];
  char out_name[500];
  char out_name2[5000];

  DIR *theDir;
  struct dirent *dirEntry;
  int nFrames;
  size_t prefixLength;
  double count;
  int out_count;
  int in_count;
  char test_c;
  

  program_name = argv[0]; /*save the program name in case we need it*/
  


  for( i = 1; i < argc; i++){
    
    if(argv[i][0] =='-'){

      // parse the option
      
      if(argv[i][1] == '-' ){

	// parse long word options
	
	fprintf( stderr, 
		 "Invalid option \"%s\"\n", argv[i] );
	usage();
	
      }
      
      else{
	
	// parse single character options
	
	done = 0;
	j = 1;
	current_flag = argv[i][j];
	while( (current_flag != '\0') && (!done) ){
	  
	  switch(current_flag){

	  case 'c':	
	    // -c => turn on the compression flag for x-povray

	    compression = 1;
	    strcat(pov_flags, " +FP ");
	    break;
	    
	  case 'd':
	    // -d <#> => duplicate the frames 

	    duplicate = 1;
	    i++;
	    n_duplicates = atoi(argv[i]);
	    done =1;
	    break;


	
	  case 'f':
	    // -f <#> => tells us the final frame

	    is_end_frame = 1;
	    i++;
	    end_frame = atoi(argv[i]);
	    done=1;
	    break;

	  case 's':
	    // -s <#> => tells the start frame

	    i++;
	    start_frame = atoi(argv[i]);
	    done=1;
	    break;

	  case 'F':
	    // -F <flags> => The flag arguments to be passed to x-povray

	    i++;
	    strcat(pov_flags, argv[i]);
	    done=1;
	    break;
	    
	  default:

	    (void)fprintf(stderr, "Bad option \"-%c\"\n", current_flag);
	    usage();
	  }
	  j++;
	  current_flag = argv[i][j];
	}
      }
    }

    else{
      
      if( prefix != NULL ){
	fprintf( stderr, 
		 "Error at \"%s\", program does not currently support\n"
		 "more than one prefix.\n"
		 "\n",
		 argv[i]);
	usage();
      }
      
      prefix = argv[i];
    }
  }


  if(prefix == NULL){
    usage();
  }
  

  if(access(ANIM_DIR, F_OK)){
    /*create the anim directory*/
    mkdir(ANIM_DIR, dir_mode);
  }

  if(!compression) strcat(pov_flags, " +FT ");

  strcpy(anim_dir, ANIM_DIR); strcat(anim_dir, "/");
  strcpy(pov_dir, POV_DIR); strcat(pov_dir, "/");

  // find the number of digits to output
  
  nFrames = 0;
  if( is_end_frame ){
    nFrames = end_frame - start_frame;
  }
  else{

    theDir = opendir( POV_DIR );
    dirEntry = readdir( theDir );
    prefixLength = strlen( prefix );
    while( dirEntry != NULL ){
      
      if( !strncmp( dirEntry->d_name, prefix, prefixLength ) ) nFrames++;
      
      dirEntry = readdir( theDir );
    }
    closedir( theDir );
  }
  
  out_count = 1;
  count = (double)( nFrames * (n_duplicates+1) );
  while( count >= 10.0 ){
    count /= 10.0;
    out_count++;
  }

  // how many padded digits are there for the input
  
  theDir = opendir( POV_DIR );
  dirEntry = readdir( theDir );
  prefixLength = strlen( prefix );
  while( dirEntry != NULL ){
    
    in_count = 0;
    if( !strncmp( dirEntry->d_name, prefix, prefixLength ) ){
      
      i = prefixLength;
      test_c = dirEntry->d_name[i];
      while( isNumber( test_c ) ){
	
	in_count++;
	i++;
	test_c = dirEntry->d_name[i];
      }
      break;
    }

    dirEntry = readdir( theDir );
  }
  closedir( theDir );


  strcpy(in_name, pov_dir);
  strcpy(out_name, anim_dir);
  
  sprintf( in_format, "%s%s%%0%dd.pov", pov_dir, prefix, in_count );
  
  sprintf( out_format, "%s%s%%0%dd", anim_dir, prefix, out_count );
  if(compression) strcat(out_format, ".ppm");
  else strcat(out_format, ".tga");

  in_frame = start_frame;
  out_frame = 0;

  if(access(ANIM_DIR, F_OK)){
    /*create the anim directory*/
    mkdir(ANIM_DIR, dir_mode);
  }


  sprintf( in_name, in_format, in_frame );
  sprintf( out_name, out_format, out_frame );
    
  in_file_exists = !access(in_name, F_OK);
  
  if( is_end_frame ){
    
    while( (in_frame <= end_frame) && in_file_exists ){  
      
      
      sprintf(command, "x-povray +I %s +O %s %s", 
	      in_name, out_name, pov_flags);
      system(command);
      
      if(duplicate){
	for(k = 0; k < n_duplicates; k++){
	  out_frame++;
	  
	  sprintf( out_name2, out_format, out_frame );
	  sprintf(command, "cp %s %s", out_name, out_name2);
	  system(command);
	}
      }
      
      in_frame++;
      out_frame++;
      
      sprintf( in_name, in_format, in_frame );
      sprintf( out_name, out_format, out_frame );
      
      in_file_exists = !access(in_name, F_OK);
    }
  }
  else{
    
    while( in_file_exists ){  
      
      
      sprintf(command, "x-povray +I %s +O %s %s", 
	      in_name, out_name, pov_flags);
      system(command);
      
      if(duplicate){
	for(k = 0; k < n_duplicates; k++){
	  out_frame++;
	  
	  sprintf( out_name2, out_format, out_frame );
	  sprintf(command, "cp %s %s", out_name, out_name2);
	  system(command);
	}
      }
      
      in_frame++;
      out_frame++;
      
      sprintf( in_name, in_format, in_frame );
      sprintf( out_name, out_format, out_frame );
      
      in_file_exists = !access(in_name, F_OK);
    }
  }
  
  return 0;    
}





/***************************************************************************
 * prints out the usage for the command line arguments, then exits.
 ***************************************************************************/

void usage(){
  (void)fprintf(stderr, 
		"The proper usage is: %s [options] <pov_prefix>\n\n"
		"Options:\n"
		"   -c         turns on compression for x-povray\n"
		"   -d <#>     turns on and duplicates the number of frames\n"
		"   -s <#>     the starting frame number\n"
		"   -f <#>     the final frame number\n"
		"   -F \"<povray flags>\"\n"
		"              flags to pass to x-povray\n" 
		
		"\n",
		program_name);
  exit(8);
}


int isNumber( char c ){

  if( c == '0') return 1;
  else if( c == '1') return 1;
  else if( c == '2') return 1;
  else if( c == '3') return 1;
  else if( c == '4') return 1;
  else if( c == '5') return 1;
  else if( c == '6') return 1;
  else if( c == '7') return 1;
  else if( c == '8') return 1;
  else if( c == '9') return 1;
  
  return 0;
}
