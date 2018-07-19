/* quad.c - Quadrupole lens */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct myquad_info
{
  double G ;	// gradient
  double zlen ; //  length
  double rlen ; //   radius
} ;

static int myquad_sim(gptpar *par, double t, struct myquad_info *info) ;


void myquad_init(gptinit *init)
{
  struct myquad_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(ECS,L,G,R)\n", gptgetname(init) ) ;

  info = (struct myquad_info *)gptmalloc( sizeof(struct myquad_info) ) ;

  info->zlen = gptgetargdouble(init,1)/2. ;
  info->G    = gptgetargdouble(init,2) ;
  info->rlen = gptgetargdouble(init,3)/2. ;

  gptaddEBelement( init, myquad_sim, gptfree, GPTELEM_GLOBAL, info ) ;


  	struct axis  * fromaxis = init->paxis;
  	gpttransform fromTransform = fromaxis->a;
	

  	gpttransform ecsTransform = init->e;

  	// calc offset of dipole in the toaxis coordinate system
  	double o[3];
  	for(int i=0; i<3; i++)
  	{
  	 	o[i] = fromTransform.o[i];
  		for(int j=0; j<3; j++)
	  		o[i] += fromTransform.m[i][j] * ecsTransform.o[j];
  	}
  	printf("o = %e %e %e\n", o[0], o[1], o[2]);


	// cacl matrix  of dipole in the toaxis coordinate system

  	double m[3][3];
  	for(int i=0; i<3; i++)
  	for(int j=0; j<3; j++)
  	{
		m[i][j] =0.;
  		for(int k=0; k<3; k++)
		{
			m[i][j] += fromTransform.m[i][k] * ecsTransform.m[k][j];
		}
  	}

  	printf("d = %e %e %e\n", m[0][2], m[1][2], m[2][2]);


	double r1 = -info->rlen * m[0][0] - info->zlen * m[0][2];
	double z1 = -info->rlen * m[2][0] - info->zlen * m[2][2];
	double r2 =  info->rlen * m[0][0] - info->zlen * m[0][2];
	double z2 =  info->rlen * m[2][0] - info->zlen * m[2][2];
	double r3 =  info->rlen * m[0][0] + info->zlen * m[0][2];
	double z3 =  info->rlen * m[2][0] + info->zlen * m[2][2];
	double r4 = -info->rlen * m[0][0] + info->zlen * m[0][2];
	double z4 = -info->rlen * m[2][0] + info->zlen * m[2][2];

	FILE * fp = fopen("_.agr","w");
  	fprintf(fp, "%e %e\n",  o[2] +  z1, o[0] + r1 );
  	fprintf(fp, "%e %e\n",  o[2] +  z2, o[0] + r2 );
  	fprintf(fp, "%e %e\n",  o[2] +  z3, o[0] + r3 );
  	fprintf(fp, "%e %e\n",  o[2] +  z4, o[0] + r4 );
  	fprintf(fp, "%e %e\n",  o[2] +  z1, o[0] + r1 );
  	fprintf(fp, "&\n");
	fclose(fp);





}


static int myquad_sim(gptpar *par, double t, struct myquad_info *info)
{
  if( fabs(Z) > info->zlen ) return( 0 ) ;

  BX = info->G * Y ;
  BY = info->G * X ;
  BY = 0.001;

  return( 1 ) ;
}
