/* quad.c - Quadrupole lens */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#define XTIMEPLOT
#include <stdio.h>
#include <math.h>
#include "elem.h"

struct kasten_info
{
  double G ;	// gradient
  double zlen ; //  length
  double rlen ; //   radius
  double time ; //  
  int samples;
} ;

// static int kasten_sim(gptpar *par, double t, struct kasten_info *info) ;

static FILE * fkast = 0;
static FILE * fline = 0;
static FILE * fcent = 0;
#ifdef TIMEPLOT
static int kasten_sim(gptpar *par, double t, struct kasten_info *info);
static void kasten_exit(struct kasten_info *info) ;
#endif

void kasten_init(gptinit *init)
{
  struct kasten_info *info = (struct kasten_info *)gptmalloc( sizeof(struct kasten_info) ) ;

  gptbuildECS( init ) ;

  int  numc = gptgetargnum(init);
  if( numc !=2 && numc != 3)
    gpterror( "Syntax: %s(ECS,L,R)\n", gptgetname(init) ) ;


  info->zlen = gptgetargdouble(init,1)/2. ;
  info->rlen = gptgetargdouble(init,2)/2. ;
  info->time=0.;
  info->samples=0;
  const char * tag  = (numc == 3) ? gptgetargstring(init,3) : "";



  if(!fkast)
  {
	fkast = fopen("kasten.agr", "w");
  	if(!fkast)
  	{
		perror("kasten.agr");
		gpterror( "can't open kasten.agr\n");
	}
  }

  if(!fline)
  {
	fline = fopen("line.agr", "w");
  	if(!fline)
  	{
		perror("line.agr");
		gpterror( "can't open line.agr\n");
	}
  	fprintf(fline, "0 0\n");
  }

  if(!fcent)
  {
	fcent = fopen("cent.agr", "w");
  	if(!fcent)
  	{
		perror("cent.agr");
		gpterror( "can't open cent.agr\n");
	}
  	fprintf(fcent, "0 0\n");
  }

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
//   	printf("o = %e %e %e\n", o[0], o[1], o[2]);


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

//   	printf("d = %e %e %e\n", m[0][2], m[1][2], m[2][2]);

	double zz = ecsTransform.o[2];
  	fprintf(fline, "%e %e\n",  zz - info->zlen, 0.);
  	fprintf(fline, "%e %e\n",  zz - info->zlen, info->rlen);
  	fprintf(fline, "%e %e\n",  zz + info->zlen, info->rlen);
  	fprintf(fline, "%e %e\n",  zz + info->zlen, 0.);

	double r1 = -info->rlen * m[0][0] - info->zlen * m[0][2];
	double z1 = -info->rlen * m[2][0] - info->zlen * m[2][2];
	double r2 =  info->rlen * m[0][0] - info->zlen * m[0][2];
	double z2 =  info->rlen * m[2][0] - info->zlen * m[2][2];
	double r3 =  info->rlen * m[0][0] + info->zlen * m[0][2];
	double z3 =  info->rlen * m[2][0] + info->zlen * m[2][2];
	double r4 = -info->rlen * m[0][0] + info->zlen * m[0][2];
	double z4 = -info->rlen * m[2][0] + info->zlen * m[2][2];

	fprintf(fcent, "@    s1 comment \"%s\"\n", tag);
  	fprintf(fcent, "%e %e\n",  o[2], o[0] );


	fprintf(fkast, "@    s1 comment \"%s\"\n", tag);
  	fprintf(fkast, "%e %e\n",  o[2] +  z1, o[0] + r1 );
  	fprintf(fkast, "%e %e\n",  o[2] +  z2, o[0] + r2 );
  	fprintf(fkast, "%e %e\n",  o[2] +  z3, o[0] + r3 );
  	fprintf(fkast, "%e %e\n",  o[2] +  z4, o[0] + r4 );
  	fprintf(fkast, "%e %e\n",  o[2] +  z1, o[0] + r1 );
  	fprintf(fkast, "&\n");


#ifdef TIMEPLOT
	gptaddEBelement( init, kasten_sim, kasten_exit, GPTELEM_GLOBAL, info ) ;
#endif

}

#ifdef TIMEPLOT
static int kasten_sim(gptpar *par, double t, struct kasten_info *info)
{
// 	if( fabs(Z) < info->zlen)
	{
		info->time += t;
		info->samples++;
		printf("kasten_sim  %e %e %e\n", info->zlen, Z, t);
	}
	return 0;
}
static void kasten_exit(struct kasten_info *info) 
{
	printf("kasten_exit %e %d\n", info->time,info->samples);
}
#endif
