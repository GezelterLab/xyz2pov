#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "atom_parser.h"
#include "pov_writer.h"
#include "file_stuff.h"
#define POV_DIR "./pov"

struct linked_xyz{
  struct coords *r;
  struct linked_xyz *next;
};


struct linked_coords{
  struct coords r;
  struct linked_coords *next;
};


char *program_name; /*the name of the program */
int draw_bonds = 1; /* boolean to draw bonds or not */
int draw_hydrogens = 0; /*boolean to draw hydrogens */
int draw_atoms = 0; /*boolean to draw atoms */

void usage(void);

int main(argc, argv)
     int argc;
     char *argv[];
{
  
  char name_foo[10];
  char *token;
  struct coords *out_coords;
  struct linked_coords *lc_head;
  struct linked_coords *lc_temp;
  struct linked_coords *lc_temp2;
  int i; /* loop counters */
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
  FILE *in_file; /* the input file */
  FILE *out_file; /*the output file */
  char *out_prefix = NULL; /*the prefix of the output file */
  char out_name[500]; /*the output name */
  char temp_name[500];
  char *in_name = NULL; /*the name of the input file */
  unsigned int n_out = 0; /*keeps track of which output file is being written*/

  struct linked_xyz *current_frame;
  
  double dx, dy, dz; /* temp variables for interpolating distances */

  char pov_dir[500]; /* the pov_dir */

  program_name = argv[0]; /*save the program name in case we need it*/
  
  for( i = 0; i < argc; i++){
    
    if(argv[i][0] =='-'){
      
      /* argv[i][1] is the actual option character */
      
      switch(argv[i][1]){
	
	/* -H => generate header */
	
      case 'H':
	generate_header = 1;
	break;
	
	
	/* -f <name> => the xyz input file
	 *     [i+1] actually starts the name
	 */

      case 'f':
	in_name = argv[i+1];
	break;

	/* -o <name> => the orientation file prefix
	 *     [i+1] actually starts the name
	 */
	
      case 'o':
	out_prefix = argv[i+1];
	break;
	
      default:
	(void)fprintf(stderr, "Bad option %s\n", argv[i]);
	usage();
      }
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
  
  if(out_prefix == NULL){
    out_prefix = strtok(in_name, ".");
  }
  
  if(access(POV_DIR, F_OK)){
    /*create the pov directory*/
    mkdir(POV_DIR, dir_mode);
  }
  strcpy(pov_dir, POV_DIR); strcat(pov_dir, "/");

  
  // initialize atom type parser

  initializeParser();

  // start reading the frame

  lc_head = (struct linked_coords *)malloc(sizeof(struct linked_coords));
  lc_head->next = NULL;
  
  eof_test = fgets(read_buffer, sizeof(read_buffer), in_file);
  
  current_frame = (struct linked_xyz *)malloc(sizeof(struct linked_xyz));
  current_frame->next = NULL;
  

  while(eof_test != NULL){
    
    token = strtok(read_buffer, "\t ");
    token = strtok(NULL, "\t ");
    
    (void)sscanf(token, "%d", &n_atoms);
    
    token = strtok(NULL, "\t ");
    
    (void)sscanf(token, "%s", &name_foo);
    
    i=0;
    while( isalpha( name_foo[i] ) ){
      
      lc_head->r.name[i] = name_foo[i];
      i++;
    }
    lc_head->r.name[i] = '\0';
    
    token = strtok(NULL, "\t ");
    token = strtok(NULL, "\t ");
    
    token = strtok(NULL, "\t ");
    (void)sscanf(token, "%lf", &lc_head->r.x);

    token = strtok(NULL, "\t ");
    (void)sscanf(token, "%lf", &lc_head->r.y);

    token = strtok(NULL, "\t ");
    (void)sscanf(token, "%lf", &lc_head->r.z);
    
    lc_temp = (struct linked_coords *)malloc(sizeof(struct linked_coords));
    
    lc_temp->next = lc_head;
    lc_head = lc_temp;

    eof_test = fgets(read_buffer, sizeof(read_buffer), in_file);
  }

  current_frame->r = 
    (struct coords *)calloc(n_atoms, sizeof(struct coords));

  lc_temp = lc_head->next;
  
  for( i=0; i < n_atoms; i++ ){
    
    strcpy(current_frame->r[i].name, lc_temp->r.name);

    current_frame->r[i].x = lc_temp->r.x;
    if(current_frame->r[i].x > big_x) big_x = current_frame->r[i].x;
    if(current_frame->r[i].x < small_x) small_x = current_frame->r[i].x;

    current_frame->r[i].y = lc_temp->r.y;
    if(current_frame->r[i].y > big_y) big_y = current_frame->r[i].y;
    if(current_frame->r[i].y < small_y) small_y = current_frame->r[i].y;

    current_frame->r[i].z = lc_temp->r.z;
    if(current_frame->r[i].z > big_z) big_z = current_frame->r[i].z;
    if(current_frame->r[i].z < small_z) small_z = current_frame->r[i].z;
    
    lc_temp = lc_temp->next;
  }

  
  /*clean up the memory */
  
  lc_temp = lc_head->next;
  
  while(lc_temp != NULL){
    lc_temp2 = lc_temp->next;
    free(lc_temp);
    lc_temp = lc_temp2;
  }

  
  /* open the new output file */
  
  strcpy(out_name, pov_dir);
  make_filename(temp_name, out_prefix, n_out);
  strcat(out_name, temp_name); strcat(out_name, ".pov");
  out_file = fopen(out_name, "w");
  n_out++;
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
  
  (void)fclose(in_file);


  if(generate_header){
    
    dx = big_x - small_x;
    dy = big_y - small_y;
    dz = big_z - small_z;
    
    rsqr = dx * dx + dy * dy + dz * dz;
    diagonal = sqrt(rsqr);
    diagonal *= 0.5;
    
    dx /= 2.0;
    dy /= 2.0;
    dz /= 2.0;
    
    dx += small_x;
    dy += small_y;
    dz += small_z;
    

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
	    "#declare zoom = %lf;\n"
	    "\n"
	    "#declare ATOM_SPHERE_FACTOR = 0.2;\n"
	    "#declare BOND_RADIUS = 0.1;\n"
	    "\n"
	    "camera{\n"
	    "  location < %lf, %lf, %lf - zoom>\n"
	    "  right < Ratio , 0, 0>\n"
	    "  look_at < %lf, %lf, %lf >\n" 
	    "}\n"
	    "\n",
	    diagonal,
	    dx, dz, dy,
	    dx, dz, dy);

    fprintf(out_file,
	    "light_source{\n"
	    "  < %lf, %lf, %lf - zoom >\n"
	    "  rgb < 1.0, 1.0, 1.0 > }\n",
	    dx, dz, dy);
	
    fprintf(out_file,
	    "light_source{\n"
	    "  < %lf - zoom , %lf + zoom, %lf - zoom >\n"
	    "  rgb < 1.0, 1.0, 1.0 > }\n"
	    "\n"
	    "\n",
	    dx, dz, dy);

    fprintf(out_file,
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
		"The proper usage is: %s [options] -f <xyz_file>\n\n"
		"Options:\n"
		"   -o <name>  the output file prefix\n"
		"   -H         generate a pov-ray header file\n",
		program_name);
  exit(8);
}
