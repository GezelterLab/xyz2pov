#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pov_writer.h"
#include "atom_parser.h"

struct bond{
  int i;
  int j;
};

struct linked_bond_list{
  struct bond the_bond;
  struct linked_bond_list *next;
};

void make_bonds(struct coords *, int);

struct linked_bond_list *bl_head;
  
void clean_bonds(void);

void initBondList(void){
  bl_head = NULL;
}

void pov_write(FILE *out_file, struct coords *the_coords, int n_atoms,
	       int d_hydrogens, int d_bonds, int d_atoms, int d_vectors){

  int i,j; /*loop counters */
  int skip_atom, skip_bond, test1, test2; /*booleans */
  int gb_atom, gbdp_atom;
  double dx, dy, dz; /* used in making the bonds */
  
  struct linked_bond_list *current_bond; /*keeps track of the linked list*/

  for(i = 0; i < n_atoms; i++){
    update_types(the_coords[i].name);
  }  
  
  if(d_atoms){
    
    fprintf(out_file,
	    "//************************************************************\n"
	    "// The list of atoms\n"
	    "//************************************************************\n"
	    "\n"
	    "\n");
    
    for(i = 0; i < n_atoms; i++){
      
      skip_atom = 0;     
      
      if(!d_hydrogens){
	skip_atom = !strcmp("H", the_coords[i].name);
      }
      
      if(!skip_atom){	       
        
        gb_atom = !strcmp("GB", the_coords[i].name);
        gbdp_atom = !strcmp("GBDP", the_coords[i].name);
        
        if (gb_atom) {
          fprintf(out_file, 
                  "make_%s_ellipse( %lf, %lf, %lf, %lf, %lf, %lf, %lf)\n",
                  the_coords[i].name,
                  the_coords[i].x,
                  the_coords[i].z,
                  the_coords[i].y,
                  the_coords[i].charge,
                  the_coords[i].ux,
                  the_coords[i].uz,
                  the_coords[i].uy);
        } else {
          if (gbdp_atom) {
            fprintf(out_file, 
                    "make_%s_shaded_ellipse( %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf)\n",
                    the_coords[i].name,
                    the_coords[i].x,
                    the_coords[i].z,
                    the_coords[i].y,
                    the_coords[i].charge,
                    the_coords[i].ux,
                    the_coords[i].uz,
                    the_coords[i].uy,
                    the_coords[i].vx,
                    the_coords[i].vz,
                    the_coords[i].vy);
          } else {
            fprintf(out_file, 
                    "make_%s_atom( %lf, %lf, %lf )\n",
                    the_coords[i].name,
                    the_coords[i].x,
                    the_coords[i].z,
                    the_coords[i].y);
          }
        }
      }
    }
    
    fprintf(out_file,
	    "\n"
	    "\n");
  }

      
  if (d_vectors) {

    fprintf(out_file,
	    "//************************************************************\n"
	    "// The list of vectors\n"
	    "//************************************************************\n"
	    "\n"
	    "\n");
    
    for(i = 0; i < n_atoms; i++){
      
      if (the_coords[i].hasVector) {
        fprintf(out_file,
                "make_%s_vector(%lf, %lf, %lf, %lf, %lf, %lf)\n",
                the_coords[i].name,
                the_coords[i].x,
                the_coords[i].z,
                the_coords[i].y,
                the_coords[i].ux,
                the_coords[i].uz,
                the_coords[i].uy);
      }
    }
    
    fprintf(out_file,
            "\n"
            "\n");
  }
      
  if(d_bonds){
    
    fprintf(out_file,
	    "//************************************************************\n"
	    "// The list of bonds\n"
	    "//************************************************************\n"
	    "\n"
	    "\n");
    
    if( bl_head == NULL ) make_bonds(the_coords, n_atoms);
    
    current_bond = bl_head->next;
    
    while(current_bond != NULL){
      
      skip_bond = 0;
      
      i = current_bond->the_bond.i;
      j = current_bond->the_bond.j;
      
      if(!d_hydrogens){

	test1 = !strcmp("H", the_coords[i].name);
	test2 = !strcmp("H", the_coords[j].name);

	skip_bond = (test1 || test2);
      }

      if(!skip_bond){
	
	dx = (the_coords[j].x - the_coords[i].x) / 2.0;
	dy = (the_coords[j].y - the_coords[i].y) / 2.0;
	dz = (the_coords[j].z - the_coords[i].z) / 2.0;
	
	fprintf(out_file,
	       "make_%s_bond( %lf, %lf, %lf, %lf, %lf, %lf )\n",
	       the_coords[i].name,
	       the_coords[i].x,
	       the_coords[i].z,
	       the_coords[i].y,
	       (the_coords[i].x + dx),
	       (the_coords[i].z + dz),
	       (the_coords[i].y + dy));
	
	fprintf(out_file,
	       "make_%s_bond( %lf, %lf, %lf, %lf, %lf, %lf )\n",
	       the_coords[j].name,
	       the_coords[j].x,
	       the_coords[j].z,
	       the_coords[j].y,
	       (the_coords[j].x - dx),
	       (the_coords[j].z - dz),
	       (the_coords[j].y - dy));

	fprintf(out_file, "\n");
      }
      
      current_bond = current_bond->next;
    }
  
    if( regenerateBonds )clean_bonds();
  }
}


void make_bonds(struct coords *the_coords, int n_atoms){
  
  int i, j; /*counters */
  struct linked_bond_list *temp_bond; /*bond place holder */
  
  const double bond_fudge = 1.12; // a fudge factor
  struct atom type_i, type_j; /* holds the atom types */

  int test; /* booleans */
  double dx, dy, dz, dr2, dcv, dcv2; // used to determine bond existence 

  
  bl_head = (struct linked_bond_list *)malloc(sizeof(struct linked_bond_list));
  bl_head->next = NULL;

  for(i = 0; i < (n_atoms - 1); i++){
    
    for(j = (i+1); j < n_atoms; j++){
      
      dx = the_coords[j].x - the_coords[i].x;
      dy = the_coords[j].y - the_coords[i].y;
      dz = the_coords[j].z - the_coords[i].z;
      
      dr2 = dx * dx + dy * dy + dz * dz;
      
      test = !findAtomType(the_coords[i].name, &type_i);
      if(test){
	fprintf(stderr, "Atom Type %s, not found!\n",
		the_coords[i].name);
	exit(8);
      }

      test = !findAtomType(the_coords[j].name, &type_j);
      if(test){
	fprintf(stderr, "Atom Type %s, not found!\n",
		the_coords[j].name);
	exit(8);
      }
      
      
      dcv = bond_fudge * (type_i.covalentRadii + type_j.covalentRadii);
      dcv2 = dcv * dcv;

      if(dr2 <= dcv2){
	
	temp_bond = 
	  (struct linked_bond_list *)malloc(sizeof(struct linked_bond_list));
	
	bl_head->the_bond.i = i;
	bl_head->the_bond.j = j;
	
	temp_bond->next = bl_head;
	bl_head = temp_bond;
      }
    }
  }
}
	

void clean_bonds(){
  struct linked_bond_list *current_bond;
  struct linked_bond_list *next_bond; /* place holders */
  

    current_bond = bl_head->next;

    while(current_bond != NULL){

      next_bond = current_bond->next;
      free(current_bond);
      current_bond = next_bond;
    }
    
    bl_head->next = NULL;
    free( bl_head );
    bl_head = NULL;
}


void make_header_macros(FILE *out_file){

  struct linked_atom *type_list;  // list of all atom types
  struct linked_atom *current_type; // current atom type

  char *name;
  double red, green, blue;
  double radius;
  
  type_list = get_type_list();
  current_type = type_list->next;
  
  while(current_type != NULL){
    
    name = current_type->myAtom.name;
    radius = current_type->myAtom.vanDerWallRadii;
    red = ((double)current_type->myAtom.red) / 255.0;
    green =  ((double)current_type->myAtom.green) / 255.0;
    blue =  ((double)current_type->myAtom.blue) / 255.0;

    
    
    fprintf(out_file,
	    "//****************************************************\n"
	    "// DEFINE %s MACROS\n" 
	    "//****************************************************\n"
	    "\n"
	    "#macro make_%s_bond "
	    "(end_1x, end_1y, end_1z, end_2x, end_2y, end_2z)\n"
	    "\n"
	    "  #local x1 = end_1x;\n"
	    "  #local y1 = end_1y;\n"
	    "  #local z1 = end_1z;\n"
	    "  #local x2 = end_2x;\n"
	    "  #local y2 = end_2y;\n"
	    "  #local z2 = end_2z;\n"
	    "\n"
	    "  #if(ROTATE)\n"
 	    "    #local x1_new = rotatePointX + A11 * (x1-rotatePointX) + A12 * (y1-rotatePointY) + A13 * (z1-rotatePointZ);\n"
	    "    #local y1_new = rotatePointY + A21 * (x1-rotatePointX) + A22 * (y1-rotatePointY) + A23 * (z1-rotatePointZ);\n"
	    "    #local z1_new = rotatePointZ + A31 * (x1-rotatePointX) + A32 * (y1-rotatePointY) + A33 * (z1-rotatePointZ);\n"
	    "\n"
 	    "    #local x2_new = rotatePointX + A11 * (x2-rotatePointX) + A12 * (y2-rotatePointY) + A13 * (z2-rotatePointZ);\n"
	    "    #local y2_new = rotatePointY + A21 * (x2-rotatePointX) + A22 * (y2-rotatePointY) + A23 * (z2-rotatePointZ);\n"
	    "    #local z2_new = rotatePointZ + A31 * (x2-rotatePointX) + A32 * (y2-rotatePointY) + A33 * (z2-rotatePointZ);\n"
	    "\n"
	    "  #else\n"
 	    "    #local x1_new = x1;"
	    "    #local y1_new = y1;"
	    "    #local z1_new = z1;"
	    "\n"
 	    "    #local x2_new = x2;"
	    "    #local y2_new = y2;"
	    "    #local z2_new = z2;"
	    "\n"
	    "  #end\n"
	    "\n"
	    "  cylinder{\n"
	    "    < x1_new, y1_new, z1_new >,\n"
	    "    < x2_new, y2_new, z2_new >,\n"
	    "    BOND_RADIUS\n"
	    "    texture{\n"
	    "      pigment{ rgb < %lf, %lf, %lf > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "#end\n"
	    "#macro make_%s_atom "
	    "(center_x, center_y, center_z)\n"
	    	    "\n"
	    "  #local x1 = center_x;\n"
	    "  #local y1 = center_y;\n"
	    "  #local z1 = center_z;\n"
	    "\n"
	    "  #if(ROTATE)\n"
	    "\n"
 	    "    #local x1_new = rotatePointX + A11 * (x1-rotatePointX) + A12 * (y1-rotatePointY) + A13 * (z1-rotatePointZ);\n"
	    "    #local y1_new = rotatePointY + A21 * (x1-rotatePointX) + A22 * (y1-rotatePointY) + A23 * (z1-rotatePointZ);\n"
	    "    #local z1_new = rotatePointZ + A31 * (x1-rotatePointX) + A32 * (y1-rotatePointY) + A33 * (z1-rotatePointZ);\n"
	    "\n"
	    "  #else\n"
	    "\n"
 	    "    #local x1_new = x1;"
	    "    #local y1_new = y1;"
	    "    #local z1_new = z1;"
	    "\n"
	    "  #end\n"
	    "\n"
	    "  sphere{\n"
	    "    < x1_new, y1_new, z1_new >,\n"
	    "    ATOM_SPHERE_FACTOR * %lf\n"
	    "    texture{\n"
	    "      pigment{ rgb < %lf, %lf, %lf > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "#end\n"
	    "#macro make_%s_ellipse "
	    "(center_x, center_y, center_z, ecc, u_x, u_y, u_z)\n"
	    	    "\n"
	    "  #local x1 = center_x;\n"
	    "  #local y1 = center_y;\n"
	    "  #local z1 = center_z;\n"
	    "  #local x2 = u_x;\n"
	    "  #local y2 = u_y;\n"
	    "  #local z2 = u_z;\n"
	    "\n"
	    "  #if(ROTATE)\n"
	    "\n"
 	    "    #local x1_new = rotatePointX + A11 * (x1-rotatePointX) + A12 * (y1-rotatePointY) + A13 * (z1-rotatePointZ);\n"
	    "    #local y1_new = rotatePointY + A21 * (x1-rotatePointX) + A22 * (y1-rotatePointY) + A23 * (z1-rotatePointZ);\n"
	    "    #local z1_new = rotatePointZ + A31 * (x1-rotatePointX) + A32 * (y1-rotatePointY) + A33 * (z1-rotatePointZ);\n"
	    "\n"
 	    "    #local x2_new = rotatePointX + A11 * (x2-rotatePointX) + A12 * (y2-rotatePointY) + A13 * (z2-rotatePointZ);\n"
	    "    #local y2_new = rotatePointY + A21 * (x2-rotatePointX) + A22 * (y2-rotatePointY) + A23 * (z2-rotatePointZ);\n"
	    "    #local z2_new = rotatePointZ + A31 * (x2-rotatePointX) + A32 * (y2-rotatePointY) + A33 * (z2-rotatePointZ);\n"
	    "\n"
	    "  #else\n"
	    "\n"
 	    "    #local x1_new = x1;"
	    "    #local y1_new = y1;"
	    "    #local z1_new = z1;"
	    "\n"
 	    "    #local x2_new = x2;"
	    "    #local y2_new = y2;"
	    "    #local z2_new = z2;"
	    "\n"
	    "  #end\n"
	    "\n"
            "    #local myUlen = sqrt(x2_new*x2_new + y2_new*y2_new + z2_new*z2_new);\n"
            "    #local uux = x2_new / myUlen;\n"
            "    #local uuy = y2_new / myUlen;\n"
            "    #local uuz = z2_new / myUlen;\n"
            "    #local myTheta = -degrees(acos(uuz));\n"
            "    #local myPsi = -degrees(atan(uux/uuy));\n"
            "    #local myScale = ATOM_SPHERE_FACTOR * %lf;\n"
            "\n"
	    "  sphere{\n"
	    "    < 0, 0, 0 >, 1\n"
	    "    texture{\n"
	    "      pigment{\n"
	    "        average\n"
	    "        pigment_map{\n"
	    "          [1.0 grad1]\n"
	    "          [1.0 grad2]\n"
	    "          [1.0 grad3]\n"
	    "          [5.0 gradz]\n"
	    "        }\n"
	    "      }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "    scale<myScale,myScale,ecc*myScale>\n"
            "    rotate<myTheta,0,myPsi>\n"
            "    translate< x1_new, y1_new, z1_new>\n"
	    "  }\n"
	    "#end\n"
	    "#macro make_%s_shaded_ellipse "
	    "(center_x, center_y, center_z, ecc, u_x, u_y, u_z, v_x, v_y, v_z)\n"
	    	    "\n"
	    "  #local x1 = center_x;\n"
	    "  #local y1 = center_y;\n"
	    "  #local z1 = center_z;\n"
	    "  #local x2 = u_x;\n"
	    "  #local y2 = u_y;\n"
	    "  #local z2 = u_z;\n"
	    "  #local x3 = v_x;\n"
	    "  #local y3 = v_y;\n"
	    "  #local z3 = v_z;\n"
	    "\n"
	    "  #if(ROTATE)\n"
	    "\n"
 	    "    #local x1_new = rotatePointX + A11 * (x1-rotatePointX) + A12 * (y1-rotatePointY) + A13 * (z1-rotatePointZ);\n"
	    "    #local y1_new = rotatePointY + A21 * (x1-rotatePointX) + A22 * (y1-rotatePointY) + A23 * (z1-rotatePointZ);\n"
	    "    #local z1_new = rotatePointZ + A31 * (x1-rotatePointX) + A32 * (y1-rotatePointY) + A33 * (z1-rotatePointZ);\n"
	    "\n"
 	    "    #local x2_new = rotatePointX + A11 * (x2-rotatePointX) + A12 * (y2-rotatePointY) + A13 * (z2-rotatePointZ);\n"
	    "    #local y2_new = rotatePointY + A21 * (x2-rotatePointX) + A22 * (y2-rotatePointY) + A23 * (z2-rotatePointZ);\n"
	    "    #local z2_new = rotatePointZ + A31 * (x2-rotatePointX) + A32 * (y2-rotatePointY) + A33 * (z2-rotatePointZ);\n"
	    "\n"
 	    "    #local x3_new = rotatePointX + A11 * (x3-rotatePointX) + A12 * (y3-rotatePointY) + A13 * (z3-rotatePointZ);\n"
	    "    #local y3_new = rotatePointY + A21 * (x3-rotatePointX) + A22 * (y3-rotatePointY) + A23 * (z3-rotatePointZ);\n"
	    "    #local z3_new = rotatePointZ + A31 * (x3-rotatePointX) + A32 * (y3-rotatePointY) + A33 * (z3-rotatePointZ);\n"
	    "\n"
	    "  #else\n"
	    "\n"
 	    "    #local x1_new = x1;"
	    "    #local y1_new = y1;"
	    "    #local z1_new = z1;"
	    "\n"
 	    "    #local x2_new = x2;"
	    "    #local y2_new = y2;"
	    "    #local z2_new = z2;"
	    "\n"
 	    "    #local x3_new = x3;"
	    "    #local y3_new = y3;"
	    "    #local z3_new = z3;"
	    "\n"
	    "  #end\n"
	    "\n"
            "    #local myUlen = sqrt(x2_new*x2_new + y2_new*y2_new + z2_new*z2_new);\n"
            "    #local uux = x2_new / myUlen;\n"
            "    #local uuy = y2_new / myUlen;\n"
            "    #local uuz = z2_new / myUlen;\n"
            "    #local myVlen = sqrt(x3_new*x3_new + y3_new*y3_new + z3_new*z3_new);\n"
            "    #local vvx = x3_new / myVlen;\n"
            "    #local vvy = y3_new / myVlen;\n"
            "    #local vvz = z3_new / myVlen;\n"
            "\n"
            "    #local myTheta = degrees(acos(uuz));\n"
            "    #local myPsi = -degrees(atan(uux/uuy));\n"
            "    #local myPhi = degrees(acos(vvz));\n"
            "    #local myScale = ATOM_SPHERE_FACTOR * %lf;\n"
            "\n"
	    "  sphere{\n"
	    "    < 0, 0, 0 >, 1\n"
	    "    texture{\n"
	    "      pigment{\n"
	    "        average\n"
	    "        pigment_map{\n"
	    "          [1.0 grad1]\n"
	    "          [1.0 grad2]\n"
	    "          [1.0 grad3]\n"
	    "          [5.0 gradz]\n"
	    "        }\n"
	    "      }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "    scale<myScale,myScale,ecc*myScale>\n"
            "    rotate<myTheta,myPhi,myPsi>\n"
            "    translate< x1_new, y1_new, z1_new>\n"
	    "  }\n"
	    "#end\n"
	    "#macro make_%s_vector "
	    "(center_x, center_y, center_z, ux, uy, uz)\n"
	    "\n"
            "  #local vx = VECTOR_SCALE * ux;\n"
            "  #local vy = VECTOR_SCALE * uy;\n"
            "  #local vz = VECTOR_SCALE * uz;\n"
	    "  #local x1 = center_x - 0.5 * vx;\n"
	    "  #local y1 = center_y - 0.5 * vy;\n"
	    "  #local z1 = center_z - 0.5 * vz;\n"
	    "  #local x2 = center_x + 0.5 * vx;\n"
	    "  #local y2 = center_y + 0.5 * vy;\n"
	    "  #local z2 = center_z + 0.5 * vz;\n"
            "  #local v2 = vx*vx + vy*vy + vz*vz;\n"
            "  #local vl  = sqrt(v2);\n"
	    "  #local x3 = x1 + vx * (1.0 - CONE_FRACTION);\n"
	    "  #local y3 = y1 + vy * (1.0 - CONE_FRACTION);\n"
	    "  #local z3 = z1 + vz * (1.0 - CONE_FRACTION);\n"
	    "\n"
	    "  #if(ROTATE)\n"
 	    "    #local x1_new = rotatePointX + A11 * (x1-rotatePointX) + A12 * (y1-rotatePointY) + A13 * (z1-rotatePointZ);\n"
	    "    #local y1_new = rotatePointY + A21 * (x1-rotatePointX) + A22 * (y1-rotatePointY) + A23 * (z1-rotatePointZ);\n"
	    "    #local z1_new = rotatePointZ + A31 * (x1-rotatePointX) + A32 * (y1-rotatePointY) + A33 * (z1-rotatePointZ);\n"
	    "\n"
 	    "    #local x2_new = rotatePointX + A11 * (x2-rotatePointX) + A12 * (y2-rotatePointY) + A13 * (z2-rotatePointZ);\n"
	    "    #local y2_new = rotatePointY + A21 * (x2-rotatePointX) + A22 * (y2-rotatePointY) + A23 * (z2-rotatePointZ);\n"
	    "    #local z2_new = rotatePointZ + A31 * (x2-rotatePointX) + A32 * (y2-rotatePointY) + A33 * (z2-rotatePointZ);\n"
	    "\n"
 	    "    #local x3_new = rotatePointX + A11 * (x3-rotatePointX) + A12 * (y3-rotatePointY) + A13 * (z3-rotatePointZ);\n"
	    "    #local y3_new = rotatePointY + A21 * (x3-rotatePointX) + A22 * (y3-rotatePointY) + A23 * (z3-rotatePointZ);\n"
	    "    #local z3_new = rotatePointZ + A31 * (x3-rotatePointX) + A32 * (y3-rotatePointY) + A33 * (z3-rotatePointZ);\n"
	    "\n"
	    "  #else\n"
 	    "    #local x1_new = x1;"
	    "    #local y1_new = y1;"
	    "    #local z1_new = z1;"
	    "\n"
 	    "    #local x2_new = x2;"
	    "    #local y2_new = y2;"
	    "    #local z2_new = z2;"
	    "\n"
 	    "    #local x3_new = x3;"
	    "    #local y3_new = y3;"
	    "    #local z3_new = z3;"
	    "\n"
	    "  #end\n"
	    "\n"
	    "  cylinder{\n"
	    "    < x1_new, y1_new, z1_new >,\n"
	    "    < x3_new, y3_new, z3_new >,\n"
	    "    STICK_RADIUS\n"
	    "    texture{\n"
	    "      pigment{ rgb < %lf, %lf, %lf > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "  cone{\n"
	    "    < x2_new, y2_new, z2_new >, 0.0\n"
	    "    < x3_new, y3_new, z3_new >, CONE_RADIUS\n"
	    "    texture{\n"
	    "      pigment{ rgb < %lf, %lf, %lf > }\n"
	    "      finish{\n"
	    "        ambient .2\n"
	    "        diffuse .6\n"
	    "        specular 1\n"
	    "        roughness .001\n"
	    "        metallic\n"
	    "      }\n"
	    "    }\n"
	    "  }\n"
	    "#end\n"
	    "\n"
	    "\n",
	    name,
	    name,
	    red, green, blue,
	    name,
	    radius,
	    red, green, blue,
	    name,
	    radius,
	    name,
	    radius,
	    name,
	    red, green, blue,
	    red, green, blue);
    
    current_type = current_type->next;
  }
}
