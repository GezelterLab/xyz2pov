#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atom_parser.h"


struct linked_atom *head_main;
struct linked_atom *head_inUse;
int n_inUse;

int findAtomType_internal(char *, struct atom *);

void initializeParser(){
  char *token;
  char *delim = "\t ,;";
  struct linked_atom *temp_linked;
  const int buffer_size = 300;
  char lineBuffer[buffer_size]; 
  
  char *in_filename;
  char *in_file_env = "_xyz2pov_AtomTypes_";
  FILE *in_file;

  head_main = (struct linked_atom *)malloc(sizeof(struct linked_atom));
  head_main->next = NULL;
  
  head_inUse = (struct linked_atom *)malloc(sizeof(struct linked_atom));
  head_inUse->next = NULL;

  n_inUse = 0;

  in_filename = getenv(in_file_env);

  in_file = fopen(in_filename, "r");

  if(in_file == NULL){

    printf("Error reading AtomTypes file, is the _xyz2pov_AtomTypes_\n"
	   "enviorment variable set?\n");
    exit(8);
  }
  
  while(fgets(lineBuffer, sizeof(lineBuffer), in_file) != NULL){
    
    if(lineBuffer[0] == '#' || lineBuffer[0] == '\n'){
      continue;
    }

    token = strtok(lineBuffer, delim);
    sscanf(token, "%s", head_main->myAtom.name);
    strtok(NULL, delim);
    strtok(NULL, delim);
    strtok(NULL, delim);
    token = strtok(NULL, delim);
    sscanf(token, "%lf", &head_main->myAtom.vanDerWallRadii);
    token = strtok(NULL, delim);
    sscanf(token, "%lf", &head_main->myAtom.covalentRadii);
    token = strtok(NULL, delim);
    sscanf(token, "%d", &head_main->myAtom.red);
    token = strtok(NULL, delim);
    sscanf(token, "%d", &head_main->myAtom.green);
    token = strtok(NULL, delim);
    sscanf(token, "%d", &head_main->myAtom.blue);
    
    temp_linked = (struct linked_atom *)malloc(sizeof(struct linked_atom));
    temp_linked->next = head_main;
    head_main = temp_linked;
  }

  fclose(in_file);
  return;
}
  

int findAtomType_internal(char *typeKey, struct atom *dummy_plug){

  struct linked_atom *link;

  link = head_main->next;

  while(link != NULL){

    if(!strcmp(link->myAtom.name, typeKey)){
      strcpy(dummy_plug->name, link->myAtom.name);
      dummy_plug->vanDerWallRadii = link->myAtom.vanDerWallRadii;
      dummy_plug->covalentRadii = link->myAtom.covalentRadii;
      dummy_plug->red = link->myAtom.red;
      dummy_plug->green = link->myAtom.green;
      dummy_plug->blue = link->myAtom.blue;
      
      return 1;
    }
    link = link->next;
  }
  return 0;
}


void update_types(char *new_key){
  
  int found = 0;
  
  struct linked_atom *link;
  
  link = head_inUse->next;
  
  while(link != NULL){
    
    if(!strcmp(link->myAtom.name, new_key)){
      found = 1;
    }
  
    link = link->next;
  }

  if(!found){
    found = findAtomType_internal(new_key, &head_inUse->myAtom);
    if(!found){
      fprintf(stderr, "Atom Type %s, not found!\n",
	      new_key);
      exit(8);
    }
    
    link = (struct linked_atom *)malloc(sizeof(struct linked_atom));

    link->next = head_inUse;
    head_inUse = link;

    n_inUse++;
  }

  
}

int get_n_inUse(void){
  return n_inUse;
}


int findAtomType(char *typeKey, struct atom *dummy_plug){

  struct linked_atom *link;

  link = head_inUse->next;

  while(link != NULL){

    if(!strcmp(link->myAtom.name, typeKey)){
      strcpy(dummy_plug->name, link->myAtom.name);
      dummy_plug->vanDerWallRadii = link->myAtom.vanDerWallRadii;
      dummy_plug->covalentRadii = link->myAtom.covalentRadii;
      dummy_plug->red = link->myAtom.red;
      dummy_plug->green = link->myAtom.green;
      dummy_plug->blue = link->myAtom.blue;
      
      return 1;
    }
    link = link->next;
  }
  return 0;
}


struct linked_atom *get_type_list(){
  
  return head_inUse;
}

void clean_type_list(){

  struct linked_atom *current_atom;
  struct linked_atom *next_atom; /* place holders */
  
  
  current_atom = head_inUse->next;
  
  while(current_atom != NULL){
    
    next_atom = current_atom->next;
    free(current_atom);
    current_atom = next_atom;
  }
  
  head_inUse->next = NULL;
}
