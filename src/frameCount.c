#define _FILE_OFFSET_BITS 64

#include <stdlib.h>
#include <stdio.h>

#include "frameCount.h"


int frameCount( char* in_name ){
  
  FILE *in_file;
  int nFrames = 0;
  int i, j;
  int lineNum = 0;
  char readBuffer[2000];
  
  in_file = fopen( in_name, "r" );
  if( in_file == NULL ){
    printf( "Error opening \"%s\"\n", in_name );
    exit(8);
  }
  
  fgets( readBuffer, sizeof( readBuffer ), in_file );
  lineNum++;
  if( feof( in_file ) ){
    printf( "File ended unexpectedly at line %d\n", lineNum );
    exit(8);
  }
  
  while( !feof( in_file ) ){
    
    i = atoi(readBuffer);

    fgets( readBuffer, sizeof( readBuffer ), in_file );
    lineNum++;    
    if( feof( in_file ) ){
      printf( "File ended unexpectedly at line %d, in frame%d\n",
	      lineNum, nFrames+1 );
      exit(8);
    }

    
    for(j=0; j<i; j++){
      
      fgets( readBuffer, sizeof( readBuffer ), in_file );
      lineNum++;    
      if( feof( in_file ) ){
	printf( "File ended unexpectedly at line %d,"
		" with atom %d, in frame %d\n", lineNum, j, nFrames+1 );
	exit(8);
      }
      
    }
    
    nFrames++;
    
    fgets( readBuffer, sizeof( readBuffer ), in_file );
    lineNum++;
  }
  
  fclose( in_file );
  
  return nFrames;
}
