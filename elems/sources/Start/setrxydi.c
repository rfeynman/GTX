/* setrxydi.c - Set radial distribution in x coordinate */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "elem.h"

void setrxydist_init(gptinit *init)
{
  gptparset *set ;
  double *dist ;
  gptinitpar *par ;
  const char *name ;
  int i, len ;

  if( gptgetargnum(init)<2 )
    gpterror( "Syntax: %s(set,DIST)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  /* Set distribution */
  dist = gptgetdistribution(init,len,2,gptdist_2DR,set) ;
  for( i=0 ; i<len ; i++ )
  {
    if( dist[i]<0 ) gpterror( "%s: A negative radius is not allowed\n", gptgetname(init) ) ;
    par[i].Wr[0] = dist[i] ;
  }
  gptfree( dist ) ;
}
