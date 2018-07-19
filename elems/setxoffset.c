/* setxoffset.c - Shift x particle distribution */

/* 2008.03.27 make setxoffset.c by T. Miyajima (Cornell Univ. and KEK)
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "elem.h"

//#define OFFSET_DEBUG /* Print numpar, avg, offset */

void setxoffset_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par ;
  int i, len ;
  int j_xy;
  double offset;

  /* j_xy = 0 : setxoffset */
  /* j_xy = 1 : setyoffset */
  j_xy = 0;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(set,offset)\n", gptgetname(init) ) ;

  const char *name  = gptgetargstring(init,1) ;
  offset   = gptgetargdouble(init,2) ;


  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;


/* Debug: print numpar, avg, std */
#ifdef OFFSET_DEBUG
  printf("numpar=%d, offset=%e\n", len, offset);
#endif


  /* Scale distribution */
  for( i=0 ; i<len ; i++ ) {
    par[i].Wr[j_xy] = par[i].Wr[j_xy] + offset;
  }

}
