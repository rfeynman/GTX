/* myrmax.c - Kill particle when R is too large */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

struct myrmax_info
{
  double myrmax2 ;
  double zlen ;
} ;

static int myrmax_sim(gptpar *par,double t,struct myrmax_info *info) ;

static  FILE * fp = 0;

void myrmax_init(gptinit *init)
{
  struct myrmax_info *info ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=1 && numarg!=2 )
    gpterror( "Syntax: %s(ECS,R,[L])\n", gptgetname(init) ) ;

  info = (struct myrmax_info *)gptmalloc( sizeof(struct myrmax_info) ) ;

  double myrmax2 = gptgetargdouble(init,1) ;
  if( myrmax2 > 0.) myrmax2 *= myrmax2;
  info->myrmax2 = myrmax2;

  if( numarg==2 )
    info->zlen = gptgetargdouble(init,2)/2 ;
  else
    info->zlen = DBL_MAX ;

  if(!fp)
  {
	fp = fopen("lostAt.txt", "w");
  	if(!fp) gpterror( " cant open lostAt.txt\n");
  }
  gptaddEBelement( init, myrmax_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}

static int myrmax_sim(gptpar *par,double t,struct myrmax_info *info)
{
  if( fabs(Z) > info->zlen ) return( 0 ) ;
  if( X*X+Y*Y < info->myrmax2 ) return( 0 ) ;

  if( par->alive )
  {
  	double E = 0.51104 * (  par->G -1);
  	fprintf(fp, "%p  %e %e %e   %e %e %e   %e %e  %d\n", par,
		  par->Wr[0],
		  par->Wr[1],
		  par->Wr[2],
		  par->GBr[0],
		  par->GBr[1],
		  par->GBr[2], t, E, par->ID);
  }
  gptremoveparticle(par) ;
  return( 1 ) ;
}
