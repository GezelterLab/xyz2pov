#ifndef _POV_WRITER_
#define _POV_WRITER_

struct coords{
  double x;
  double y;
  double z;
  int hasCharge;
  double charge;
  int hasVector;
  int hasSecondVector;
  double ux;
  double uy;
  double uz;
  double vx;
  double vy;
  double vz;
  char name[30];
};

extern void pov_write(FILE *out_file, struct coords *the_coords, int n_atoms, 
		      int d_hydrogens, int d_bonds, int d_atoms, 
                      int d_vectors);

extern void make_header_macros(FILE *out_file);

extern int regenerateBonds;

extern void initBondList( void );

#endif
