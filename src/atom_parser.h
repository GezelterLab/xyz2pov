#ifndef _ATOM_PARSER_
#define _ATOM_PARSER_

struct atom {
  char name[30];
  double vanDerWallRadii;
  double covalentRadii;
  int red;
  int green;
  int blue;
};

struct linked_atom {
  struct atom myAtom;
  struct linked_atom *next;
};


extern void initializeParser(void);

extern void update_types(char *);

extern struct linked_atom *get_type_list(void);

extern void clean_type_list(void);

extern int findAtomType(char *, struct atom *);
     
extern int get_n_inUse(void);


#endif
