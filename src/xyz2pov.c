#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "frameCount.h"
#include "atom_parser.h"
#include "pov_writer.h"


#define POV_DIR "./pov"


struct linked_xyz{
  struct coords *r;
  double boxX, boxY, boxZ;
  struct linked_xyz *next;
};

char *program_name; /*the name of the program */
int draw_bonds = 0; /* boolean to draw bonds or not */
int draw_hydrogens = 0; /*boolean to draw hydrogens */
int draw_atoms = 0; /*boolean to draw atoms */
int draw_box = 0;  // boolean to draw the periodic Box
int regenerateBonds = 0; // boolean to regenearate bonds each frame

void usage(void);

int main(argc, argv)
     int argc;
     char *argv[];
{


  struct coords *out_coords;

  int i,j; /* loop counters */
  mode_t dir_mode = S_IRWXU;

  int generate_header = 0; /* boolean for generating the pov ray header */
  double big_x = 0;
  double big_y = 0; /* lets me know the biggest x y and z */
  double big_z = 0;
  double small_x = 0;
  double small_y = 0; /* lets me know the smallest x, y, z */
  double small_z = 0;
  double rsqr; /* the square of the diagonal */
  double diagonal; /* the diagonal length of the sim box */

  unsigned int n_atoms; /*the number of atoms in each time step */
  char read_buffer[120]; /*the line buffer for reading */
  char *eof_test; /*ptr to see when we reach the end of the file */
  char *foo; /*the pointer to the current string token */
  FILE *in_file; /* the input file */
  FILE *out_file; /*the output file */
  char *out_prefix = NULL; /*the prefix of the output file */
  int have_prefix = 0;
  char out_name[500]; /*the output name */
  char out_format[1000];
  char *in_name = NULL; /*the name of the input file */
  unsigned int n_out = 0; /*keeps track of which output file is being written*/
  int done;
  char current_flag;
  int nFrames;
  int nZeroes;
  double count;

  int startFrame = 1;
  int endFrame;
  int span = 0;
  int currentCount = 0;
  int haveEnd = 0;

  struct linked_xyz *current_frame;
  struct linked_xyz *temp_frame;
  
  unsigned int n_interpolate = 0; /* number of frames to interpolate */
  double dx, dy, dz; /* temp variables for interpolating distances */
  

  char pov_dir[500]; /* the pov_dir */

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

	  case 'o':
	    // -o <prefix> => the output prefix.

	    i++;
	    out_prefix = argv[i];
	    have_prefix = 1;
	    done = 1;
	    break;

	  case 'i':
	    // -i <#> => the number to interpolate

	    i++;
	    n_interpolate = atoi( argv[i] );
	    done = 1;
	    break;

	  case 'f':
	    // -f <#> => frame to render

	    i++;
	    startFrame =  atoi( argv[i] );
	    endFrame = startFrame;
	    haveEnd = 1;
	    done = 1;
	    break;

	  case 's':
	    // -s <#> => frame to start

	    i++;
	    startFrame =  atoi( argv[i] );
	    done = 1;
	    break;

	  case 'e':
	    // -e <#> => frame to end

	    i++;
	    endFrame =  atoi( argv[i] );
	    haveEnd = 1;
	    done = 1;
	    break;

	  case 'H':
 	    // -h => generate a pov-ray Header
	    
	    generate_header = 1;
	    break;

	  case 'h':
 	    // -h => draw Hydrogens
	    
	    draw_hydrogens = 1;
	    break;
	
	  case 'b':
	    // -b => draw bonds

	    draw_bonds = 1;
	    break;

	  case 'a':
	    // -a => draw the atoms

	    draw_atoms = 1;
	    break;

	  case 'r':
	    // -r => regenerate bonds

	    regenerateBonds = 1;
	    break;

	  case 'p':
	    // -r => draw periodic box

	    draw_box = 1;
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
      
      if( in_name != NULL ){
	fprintf( stderr, 
		 "Error at \"%s\", program does not currently support\n"
		 "more than one input file.\n"
		 "\n",
		 argv[i]);
	usage();
      }
      
      in_name = argv[i];
    }
  }
  
  if(in_name == NULL){
    usage();
  }
  
  
  
  in_file = fopen(in_name, "r");
  if(in_file == NULL){
    printf("Cannot open file: %s\n", in_name);
    exit(8);
  }
  

  
  if(access(POV_DIR, F_OK)){
    /*create the pov directory*/
    mkdir(POV_DIR, dir_mode);
  }
  strcpy(pov_dir, POV_DIR); strcat(pov_dir, "/");
  
  
  // initialize atom type parser
  
  initializeParser();
  initBondList();
  
  // count the number of frames 

  printf( "Counting the number of frames..." );
  fflush(stdout);

  nFrames = frameCount( in_name );

  printf( "done.\n"
	  "%d frames found\n",
	  nFrames);
  fflush(stdout);

  // create the output string

  nZeroes = 1;
  count = (double)( nFrames * (n_interpolate+1) );
  while( count >= 10.0 ){
    count /= 10.0;
    nZeroes++;
  }

  if(!have_prefix){
    out_prefix = strtok(in_name, ".");
  }

  sprintf( out_format, "%s%s%%0%dd.pov", pov_dir, out_prefix, nZeroes );

  // start reading the first frame

  eof_test = fgets(read_buffer, sizeof(read_buffer), in_file);
  
  current_frame = (struct linked_xyz *)malloc(sizeof(struct linked_xyz));
  current_frame->next = NULL;
  

  if( haveEnd ) span = endFrame - startFrame;
  done = 0;
  if( span < 0 ) done = 1;
  while( (eof_test != NULL) && !done ){
    
    (void)sscanf(read_buffer, "%d", &n_atoms);
    current_frame->r = 
      (struct coords *)calloc(n_atoms, sizeof(struct coords));

    /*read and toss the comment line */

    eof_test = fgets(read_buffer, sizeof(read_buffer), in_file);
    if(eof_test == NULL){
      printf("error in reading file\n");
      exit(8);
    }
    
    // unless we need to get the box size
    
    if( draw_box ){
      foo = strtok(read_buffer, " ,;\t");
      if(foo == NULL){
	printf("error in reading file\n");
	exit(8);
      }
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading file\n");
	exit(8);
      }
      current_frame->boxX = atof( foo );

      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading file\n");
	exit(8);
      }
      current_frame->boxY = atof( foo );

      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading file\n");
	exit(8);
      }
      current_frame->boxZ = atof( foo );
    }

    for( i=0; i < n_atoms; i++){
      
      eof_test = fgets(read_buffer, sizeof(read_buffer), in_file);
      if(eof_test == NULL){
	printf("error in reading file\n");
	exit(8);
      }

      foo = strtok(read_buffer, " ,;\t");
      if(foo == NULL){
	printf("error in reading file\n");
	exit(8);
      }
      (void)strcpy(current_frame->r[i].name, foo); /*copy the atom name */

      /* next we grab the positions */
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading file\n");
	exit(8);
      }
      (void)sscanf(foo, "%lf",&current_frame->r[i].x);
      if(current_frame->r[i].x > big_x) big_x = current_frame->r[i].x;
      if(current_frame->r[i].x < small_x) small_x = current_frame->r[i].x;
            
  
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading file\n");
	exit(8);
      }
      (void)sscanf(foo, "%lf", &current_frame->r[i].y);
      if(current_frame->r[i].y > big_y) big_y = current_frame->r[i].y;
      if(current_frame->r[i].y < small_y) small_y = current_frame->r[i].y;


      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading file\n");
	exit(8);
      }
      (void)sscanf(foo, "%lf", &current_frame->r[i].z);
      if(current_frame->r[i].z > big_z) big_z = current_frame->r[i].z;
      if(current_frame->r[i].z < small_z) small_z = current_frame->r[i].z;

    }
    currentCount++;
    
    
    if( currentCount >= startFrame ){
      if(n_interpolate && current_frame->next != NULL){
	
	temp_frame = current_frame->next;
	
	for(i = 0; i < n_interpolate; i++){
	  
	  /* open the new output file */
	  
	  sprintf(out_name, out_format, currentCount );
	  out_file = fopen(out_name, "w");
	  currentCount++;
	  if(out_file == NULL){
	    printf("error opening output file: %s\n", out_name);
	    exit(8);
	  }
	  (void)fprintf(out_file,
			"// The following script was automatically generated by:\n"
			"// xyz2pov Copyright 2001 by MATTHEW A. MEINEKE\n"
			"\n"
			"#include \"pov_header.pov\"\n"
			"\n");
	  if( draw_box ){
	    dx = current_frame->boxX - temp_frame->boxX;
	    dy = current_frame->boxY - temp_frame->boxY;
	    dz = current_frame->boxZ - temp_frame->boxZ;
	    
	    dx /= (double)(n_interpolate + 1);
	    dy /= (double)(n_interpolate + 1);
	    dz /= (double)(n_interpolate + 1);
	    
	    fprintf( out_file,
		     "makePeriodicBox( %lf, %lf, %lf )\n"
		     "\n",
		     temp_frame->boxX + dx * (i+1),
		     temp_frame->boxZ + dz * (i+1),
		     temp_frame->boxY + dy * (i+1) );
	  }
	  
	  
	  out_coords = 
	    (struct coords *)calloc(n_atoms, sizeof(struct coords));
	  
	  for(j=0; j < n_atoms; j++){
	    dx = current_frame->r[j].x - temp_frame->r[j].x;
	    dy = current_frame->r[j].y - temp_frame->r[j].y;
	    dz = current_frame->r[j].z - temp_frame->r[j].z;
	    
	    dx /= (double)(n_interpolate + 1);
	    dy /= (double)(n_interpolate + 1);
	    dz /= (double)(n_interpolate + 1);
	    
	    strcpy(out_coords[j].name, temp_frame->r[j].name);
	    out_coords[j].x = temp_frame->r[j].x + dx * (i+1);
	    out_coords[j].y = temp_frame->r[j].y + dy * (i+1);
	    out_coords[j].z = temp_frame->r[j].z + dz * (i+1);
	  }
	  
	  pov_write(out_file, out_coords, n_atoms, draw_hydrogens, draw_bonds, 
		    draw_atoms);
	  free(out_coords);
	  (void)fclose(out_file);
	}
      }
      
      /* open the new output file */
      
      sprintf(out_name, out_format, currentCount );
      out_file = fopen(out_name, "w");
      if(out_file == NULL){
	printf("error opening output file: %s\n", out_name);
	exit(8);
      }
      (void)fprintf(out_file,
		    "// The following script was automatically generated by:\n"
		    "// xyz2pov Copyright 2001 by MATTHEW A. MEINEKE\n"
		    "\n"
		    "#include \"pov_header.pov\"\n"
		    "\n");
      
      if( draw_box ){
	
	fprintf( out_file,
		 "makePeriodicBox( %lf, %lf, %lf )\n"
		 "\n",
		 current_frame->boxX,
		 current_frame->boxZ,
		 current_frame->boxY );
      }
      
      
      
      out_coords = 
	(struct coords *)calloc(n_atoms, sizeof(struct coords));
      
      for(i = 0; i < n_atoms; i++){
	strcpy(out_coords[i].name, current_frame->r[i].name);
	out_coords[i].x = current_frame->r[i].x;
	out_coords[i].y = current_frame->r[i].y;
	out_coords[i].z = current_frame->r[i].z;
      }
      pov_write(out_file, out_coords, n_atoms, draw_hydrogens, draw_bonds, 
		draw_atoms);
      free(out_coords);
      
      (void)fclose(out_file);
    }
    
    /*free up memory */

    temp_frame = current_frame->next;
    current_frame->next = NULL;

    if(temp_frame != NULL){

      free(temp_frame->r);
      free(temp_frame);
    }

    /* make a new frame */

    temp_frame = (struct linked_xyz *)malloc(sizeof(struct linked_xyz));
    temp_frame->next = current_frame;
    current_frame = temp_frame;
    
    eof_test = fgets(read_buffer, sizeof(read_buffer), in_file);

    if( haveEnd ){
      if( currentCount >= (endFrame + n_interpolate * span) ) done = 1;
    }
  }

  (void)fclose(in_file);


  if(generate_header){
    
    dx = big_x - small_x;
    dy = big_y - small_y;
    dz = big_z - small_z;

    rsqr = dx * dx + dy * dy + dz * dz;
    diagonal = sqrt(rsqr);
    diagonal *= 0.5;
    
    // calculate the center 

    dx = big_x + small_x;
    dy = big_y + small_y;
    dz = big_z + small_z;

    dx /= 2.0;
    dy /= 2.0;
    dz /= 2.0;
        

    /*note the y and z axis is exchanged for the different coordinate
      system in pov-ray*/
    
   
    out_file = fopen("pov_header.pov", "w");

    fprintf(out_file,
	    "// The following script was automatically generated by:\n"
	    "// xyz2pov Copyright 2001 by MATTHEW A. MEINEKE\n"
	    "\n"
	    "\n"
	    "background { rgb <1.0, 1.0, 1.0> }\n"
	    "\n"
	    "\n"
	    );
    
    fprintf(out_file,
	    "//******************************************************\n"
	    "// Declare the resolution, camera, and light sources.\n"
	    "//******************************************************\n"
	    "\n"
	    "// NOTE: if you plan to render at a different resoltion,\n"
	    "// be sure to update the following two lines to maintain\n"
	    "// the correct aspect ratio.\n"
	    "\n"
	    "#declare Width = 640.0;\n"
	    "#declare Height = 480.0;\n"
	    "#declare Ratio = Width / Height;\n"
	    "\n"
	    "#declare ATOM_SPHERE_FACTOR = 0.2;\n"
	    "#declare BOND_RADIUS = 0.1;\n"
	    "\n"
	    "// declare camera, light, and system variables\n"
	    "\n"
	    "#declare sysCenterX = %lf;\n"
	    "#declare sysCenterY = %lf;\n"
	    "#declare sysCenterZ = %lf;\n"
	    "\n"
	    "#declare zoom = %lf;\n"
	    "\n",
	    dx, dz, dy,
	    diagonal );
    
    fprintf(out_file,
	    "#declare cameraLookX = sysCenterX;\n"
	    "#declare cameraLookY = sysCenterY;\n"
	    "#declare cameraLookZ = sysCenterZ;\n"
	    "\n"
	    "#declare cameraX = cameraLookX;\n"
	    "#declare cameraY = cameraLookY;\n"
	    "#declare cameraZ = cameraLookZ - zoom;\n"
	    "\n"
	    "#declare lightAx = cameraX;\n"
	    "#declare lightAy = cameraY;\n"
	    "#declare lightAz = cameraZ;\n"
	    "\n"
	    "#declare lightBx = cameraX - zoom;\n"
	    "#declare lightBy = cameraY + zoom;\n"
	    "#declare lightBz = cameraZ;\n"
	    "\n"
	    "#declare boxCenterX = cameraLookX;\n"
	    "#declare boxCenterY = cameraLookY;\n"
	    "#declare boxCenterZ = cameraLookZ;\n"
	    "\n"
	    "// declare the cameras and the light sources\n"
	    "\n"
	    "camera{\n"
	    "  location < cameraX, cameraY, cameraZ>\n"
	    "  right < Ratio , 0, 0>\n"
	    "  look_at < cameraLookX, cameraLookY, cameraLookZ >\n" 
	    "}\n"
	    "\n"
	    "light_source{\n"
	    "  < lightAx, lightAy, lightAz >\n"
	    "  rgb < 1.0, 1.0, 1.0 > }\n"
	    "\n"
	    "light_source{\n"
	    "  < lightBx, lightBy, lightBz >\n"
	    "  rgb < 1.0, 1.0, 1.0 > }\n"
	    "\n"
	    "\n"
	    "//************************************************************\n"
	    "// Set whether or not to rotate the system.\n"
	    "//\n"
	    "// To Rotate, set ROTATE to 1.0 (true),\n"
	    "// Then set the Euler Angles PHI, THETA, and PSI in degrees.\n"
	    "//************************************************************\n"
	    "\n"
	    "#declare ROTATE = 0.0;\n"
	    "#declare PHI = 0.0;\n"
	    "#declare THETA = 0.0;\n"
	    "#declare PSI = 0.0;\n"
	    "\n"
	    "#if(ROTATE)\n"
	    "  #declare phi_r = radians(PHI);\n"
	    "  #declare theta_r = radians(THETA);\n"
	    "  #declare psi_r = radians(PSI);\n"
	    "\n"
	    "  #declare A11 = cos(phi_r) * cos(psi_r) - sin(phi_r) * cos(theta_r) * sin(psi_r);\n"
	    "  #declare A12 = sin(phi_r) * cos(psi_r) + cos(phi_r) * cos(theta_r) * sin(psi_r);\n"
	    "  #declare A13 = sin(theta_r) * sin(psi_r);\n"
	    "\n"
	    "  #declare A21 = -cos(phi_r) * sin(psi_r) - sin(phi_r) * cos(theta_r) * cos(psi_r);\n"
	    "  #declare A22 = -sin(phi_r) * sin(psi_r) + cos(phi_r) * cos(theta_r) * cos(psi_r);\n"
	    "  #declare A23 = sin(theta_r) * cos(psi_r);\n"
	    "\n"
	    "  #declare A31 = sin(phi_r) * sin(theta_r);\n"
	    "  #declare A32 = -cos(phi_r) * sin(theta_r);\n"
	    "  #declare A33 = cos(theta_r);\n"
	    "\n"
	    "#end\n"
	    "\n"
	    "\n"
	    "//************************************************************\n"
	    "// declare the periodic box macro\n"
	    "//************************************************************\n"
	    "\n"
	    "#macro makePeriodicBox( lengthX, lengthY, lengthZ )\n"
	    "\n"
	    "  #local addX = lengthX / 2.0;\n"
	    "  #local addY = lengthY / 2.0;\n"
	    "  #local addZ = lengthZ / 2.0;\n"
	    "\n"
	    "  #local colorR = 0.90;\n"
	    "  #local colorG = 0.91;\n"
	    "  #local colorB = 0.98;\n"
	    "\n"
	    "  #local pipeWidth = 0.4;\n"
	    "\n"
	    "  // 1\n"
	    "  cylinder{\n"
	    "    < boxCenterX-addX, boxCenterY-addY, boxCenterZ-addZ >,\n"
	    "    < boxCenterX-addX, boxCenterY+addY, boxCenterZ-addZ >,\n"
	    "    pipeWidth\n"
	    "    texture{\n"
	    "      pigment{ rgb < colorR, colorG, colorB > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "\n"
	    "  // 2\n"
	    "  cylinder{\n"
	    "    < boxCenterX-addX, boxCenterY-addY, boxCenterZ+addZ >,\n"
	    "    < boxCenterX-addX, boxCenterY+addY, boxCenterZ+addZ >,\n"
	    "    pipeWidth\n"
	    "    texture{\n"
	    "      pigment{ rgb < colorR, colorG, colorB > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "\n"
	    "  // 3\n"
	    "  cylinder{\n"
	    "    < boxCenterX+addX, boxCenterY-addY, boxCenterZ-addZ >,\n"
	    "    < boxCenterX+addX, boxCenterY+addY, boxCenterZ-addZ >,\n"
	    "    pipeWidth\n"
	    "    texture{\n"
	    "      pigment{ rgb < colorR, colorG, colorB > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "\n"
	    "  // 4\n"
	    "  cylinder{\n"
	    "    < boxCenterX+addX, boxCenterY-addY, boxCenterZ+addZ >,\n"
	    "    < boxCenterX+addX, boxCenterY+addY, boxCenterZ+addZ >,\n"
	    "    pipeWidth\n"
	    "    texture{\n"
	    "      pigment{ rgb < colorR, colorG, colorB > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "\n"
	    "  // 5\n"
	    "  cylinder{\n"
	    "    < boxCenterX-addX, boxCenterY-addY, boxCenterZ-addZ >,\n"
	    "    < boxCenterX-addX, boxCenterY-addY, boxCenterZ+addZ >,\n"
	    "    pipeWidth\n"
	    "    texture{\n"
	    "      pigment{ rgb < colorR, colorG, colorB > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "\n"
	    "  // 6\n"
	    "  cylinder{\n"
	    "    < boxCenterX-addX, boxCenterY-addY, boxCenterZ+addZ >,\n"
	    "    < boxCenterX+addX, boxCenterY-addY, boxCenterZ+addZ >,\n"
	    "    pipeWidth\n"
	    "    texture{\n"
	    "      pigment{ rgb < colorR, colorG, colorB > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "\n"
	    "  // 7\n"
	    "  cylinder{\n"
	    "    < boxCenterX+addX, boxCenterY-addY, boxCenterZ+addZ >,\n"
	    "    < boxCenterX+addX, boxCenterY-addY, boxCenterZ-addZ >,\n"
	    "    pipeWidth\n"
	    "    texture{\n"
	    "      pigment{ rgb < colorR, colorG, colorB > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "\n"
	    "  // 8\n"
	    "  cylinder{\n"
	    "    < boxCenterX+addX, boxCenterY-addY, boxCenterZ-addZ >,\n"
	    "    < boxCenterX-addX, boxCenterY-addY, boxCenterZ-addZ >,\n"
	    "    pipeWidth\n"
	    "    texture{\n"
	    "      pigment{ rgb < colorR, colorG, colorB > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "\n"
	    "  // 9\n"
	    "  cylinder{\n"
	    "    < boxCenterX-addX, boxCenterY+addY, boxCenterZ-addZ >,\n"
	    "    < boxCenterX-addX, boxCenterY+addY, boxCenterZ+addZ >,\n"
	    "    pipeWidth\n"
	    "    texture{\n"
	    "      pigment{ rgb < colorR, colorG, colorB > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "\n"
	    "  // 10\n"
	    "  cylinder{\n"
	    "    < boxCenterX-addX, boxCenterY+addY, boxCenterZ+addZ >,\n"
	    "    < boxCenterX+addX, boxCenterY+addY, boxCenterZ+addZ >,\n"
	    "    pipeWidth\n"
	    "    texture{\n"
	    "      pigment{ rgb < colorR, colorG, colorB > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "\n"
	    "  // 11\n"
	    "  cylinder{\n"
	    "    < boxCenterX+addX, boxCenterY+addY, boxCenterZ+addZ >,\n"
	    "    < boxCenterX+addX, boxCenterY+addY, boxCenterZ-addZ >,\n"
	    "    pipeWidth\n"
	    "    texture{\n"
	    "      pigment{ rgb < colorR, colorG, colorB > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "\n"
	    "  // 12\n"
	    "  cylinder{\n"
	    "    < boxCenterX+addX, boxCenterY+addY, boxCenterZ-addZ >,\n"
	    "    < boxCenterX-addX, boxCenterY+addY, boxCenterZ-addZ >,\n"
	    "    pipeWidth\n"
	    "    texture{\n"
	    "      pigment{ rgb < colorR, colorG, colorB > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "\n"
	    "#end\n"
	    "\n"
	    "\n");




    make_header_macros(out_file);

    fclose(out_file);
  }

  return 0;
  
}



/***************************************************************************
 * prints out the usage for the command line arguments, then exits.
 ***************************************************************************/

void usage(){
  (void)fprintf(stderr, 
		"The proper usage is: %s [options] <xyz_file>\n"
		"\n"
		"Options:\n"
		"\n"
		"   -o <prefix>   the output file prefix\n"
		"   -H            generate a pov-ray header file\n"
		"   -i <#>        number of frames to interpolate\n"
		"   -h            draw hydrogens\n"
		"   -b            draw bonds\n"
		"   -a            draw atoms\n"
		"   -p            draw periodic box\n"
		"   -r            regenerate bond\n"
		"   -f <#>        render frame <#> only\n"
		"   -s <#>        begin at frame <#> (inclusive)\n"
		"   -e <#>        end at frame <#> (inclusive)\n"
		"\n",
		program_name);
  exit(8);
}
