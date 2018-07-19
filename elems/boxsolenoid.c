/* boxsolenoid.c - Solenoid  with thin lense fringe fields*/

/* This file is part of GPT and belongs to:
 * Jorg Kewisch, Brookhaven National Laboratory
 * boxsolenoid(ECS, Bs, L, d)
 *	Bs = solenoid field (tesla)
 *	L  = Length of the box
 *	d  = Length of the fringe field
 *	The total length of the field is L+d
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct boxsolenoid_info
{
  double Bs ;
  double L ;
  double d ;
  double lmd2 ;
  double lmp2 ;
  double br;
} ;

static int boxsolenoid_sim(gptpar *par,double t,struct boxsolenoid_info *info) ;


void boxsolenoid_init(gptinit *init)
{
  struct boxsolenoid_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(ECS,Bs, L, d)\n", gptgetname(init) ) ;

  info = (struct boxsolenoid_info *)gptmalloc( sizeof(struct boxsolenoid_info) ) ;

  info->Bs = gptgetargdouble(init,1) ;
  info->L  = gptgetargdouble(init,2) ;
  info->d  = gptgetargdouble(init,3) ;
  info->lmd2 = (info->L - info->d) / 2.;	// from center to max z for longitudinal field
  info->lmp2 = (info->L + info->d) / 2.;	// from center to max z for radial field
  info->br = info->Bs/(2.*info->d);


  gptaddEBelement( init, boxsolenoid_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int boxsolenoid_sim(gptpar *par,double t,struct boxsolenoid_info *info)
{
  double zabs = fabs(Z);
  if( zabs > info->lmp2 ) return( 0 ) ;
  if( zabs < info->lmd2 )
  {
	  BX=0.;
	  BY=0.;
	  BZ=info->Bs;
  }
  else
  {
  	BX = info->br * X;
  	BY = info->br * Y;
  	BZ = info->Bs *  (info->lmp2-zabs)/info->d;
	if(Z < 0.)
	{
		BX = -BX;
		BY = -BY;
	}
  }

  return( 1 ) ;
}
