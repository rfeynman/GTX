/* map3D3Sym_B.c - 3D rectangular field map for the magnetic field */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 * modified for a magnet that is symmetric in all 3 planes. Only 1/8 of the field table is given, Jorg Kewisch, BNL,2015
 */

#include <algorithm>
using namespace std ;

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>

#include "elem.h"


struct catleg_info
{
	char name[20];
	double x1, z1;	// points of arc on the pole faces
	double c1, s1;	// rotation from the magnet coordinate system to the pole face system
  	double B ;	      /* Specified B field */
  	double dl ;         /* Offset for edge (difference between real and e2fective length) */
  	double b1, b2 ;     /* Parameters of Enge fringe field */ 
  	double length;

  	struct axis *toaxis ;
//   	FILE * fdiag;
// 	void print()
// 	{
// 		printf("%s\n", name);
// 		printf("x1     = %20.12e\n", x1);
// 		printf("z1     = %20.12e\n", z1);
// 		printf("c1     = %20.12e\n", c1);
// 		printf("s1     = %20.12e\n", s1);
// 		printf("B      = %20.12e\n", B );
// 		printf("dl     = %20.12e\n", dl);
// 		printf("b1     = %20.12e\n", b1);
// 		printf("b2     = %20.12e\n", b2);
// 		printf("length = %20.12e\n", length);
// 	}
} ;

// static void printTransform(const char * name, gpttransform trm)
// {
//   printf("axis %s:\n", name);
//   for(int i=0; i<3; i++)
//   {
// 	  printf("%f  %f  %f      %f\n", trm.m[i][0], trm.m[i][1], trm.m[i][2], trm.o[i]);
//   }
// }

static int catleg_sim(gptpar *par, double t, struct catleg_info *info) ;
static void catleg_exit(struct catleg_info *info) ;
static int numcat=0;

// dogleg("wcs",  0.013057, 0.000000, 38.938675,      9.848078e-01,  0.000000e+00,  -1.736482e-01,     0.000000e+00,  1.000000e+00,  0.000000e+00, "20_deg_dipole_map.bin",   -dogfield, 1, "ccs1");
// catleg("wcs",  0.013057, 0.000000, 38.938675,      9.848078e-01,  0.000000e+00,  -1.736482e-01,     0.000000e+00,  1.000000e+00,  0.000000e+00, "ccs1", 0.3, w20,  dogfield, 20.4, w20/4, 0., 0.1, 0., 1);



//	catleg(ECS, file)
void catleg_init(gptinit *init)
{
  	gptbuildECS( init ) ;
  	int argc = gptgetargnum(init) ;
  	if( argc!=10)
	{
		printf("num args %d\n", argc);
    		gpterror( "Syntax: %s(ECS, toaxis, arclength, angle, Bfield, phi_in, phi_out, dl, b1, b2, flipx)\n", gptgetname(init) ) ;
	}

/* Get axisnames */
  	const char  *toaxisname ;
  	struct axis *fromaxis, *toaxis ;
  	fromaxis= init->paxis;
  	toaxisname  =gptgetargstring(init,1) ;
  	if( (  toaxis=getaxis(  toaxisname))==NULL ) gpterror( "Specified axis \"%s\" not found\n", toaxisname ) ;

/* Fill coe2ficients for local use. */
  	double B      = gptgetargdouble(init,2) ;   /* Specified B field */
  	double length = gptgetargdouble(init,3) ;	// arc cength
  	double angle  = gptgetargdouble(init,4) ;	// delection angle
  	if(fabs(angle) < 1e-6) gpterror( "angle is less than 1e-6 rad, does not make sense\n");


/* Get angles of pole-faces */
	double phi_in  = gptgetargdouble(init,5) ;	/* in- and output angles in bend plane */
    	double phi_out = gptgetargdouble(init,6) ;
    	double dl      = gptgetargdouble(init,7) ;     /* Offset for edge (difference between real and e2fective length) */
    	double b1      = gptgetargdouble(init,8) ;     /* Parameters of Enge fringe field */ 
  	if( b1<0 ) gpterror( "%s: b1 must be nonnegative\n", gptgetname(init) ) ;
    	double b2      = gptgetargdouble(init,9) ;
  	int flipx      = gptgetargdouble(init,10) ;;	// if true flip the magnet around the longitudinal axis by 180 degree
	if(flipx) { B = -B; angle = -angle; phi_in = -phi_in; phi_out = -phi_out; }
	double rho = length/angle;	// deflection radius


  	struct catleg_info *info1= (struct catleg_info *)gptmalloc( sizeof(struct catleg_info) ) ;
  	struct catleg_info *info2= (struct catleg_info *)gptmalloc( sizeof(struct catleg_info) ) ;
  	char diagFile[200];
//   	sprintf(diagFile, "cat1Diag%d", numcat);
//   	info1->fdiag = fopen(diagFile, "w");
//   	sprintf(diagFile, "cat2Diag%d", numcat++);
//   	info2->fdiag = fopen(diagFile, "w");


	double ca = cos(angle/2.);
	double sa = sin(angle/2.);



	// info1 for Z <= 0
	strcpy(info1->name, "info1");
  	info1->B  = B;	      /* Specified B field */
  	info1->dl = dl;         /* Offset for edge (difference between real and e2fective length) */
  	info1->b1 = b1 ;     /* Parameters of Enge fringe field */ 
  	info1->b2 = b2 ;     /* Parameters of Enge fringe field */ 
	info1->length = length*0.75;

	info1->x1 = rho* (1. - ca);
	info1->z1 = -rho* sa;
	info1->c1 = cos(M_PI - angle/2. + phi_in);
	info1->s1 = sin(M_PI - angle/2. + phi_in);
	info1->toaxis = toaxis;


	// info2 for Z> 0
	strcpy(info2->name, "info2");
  	info2->B  = B;	      /* Specified B field */
  	info2->dl = dl;         /* Offset for edge (difference between real and e2fective length) */
  	info2->b1 = b1 ;     /* Parameters of Enge fringe field */ 
  	info2->b2 = b2 ;     /* Parameters of Enge fringe field */ 
	info2->length = length*0.75;

	info2->x1 = rho* (1. - ca);
	info2->z1 = -info1->z1;
	info2->c1 = cos(angle/2. - phi_out);
	info2->s1 = sin(angle/2. - phi_out);
	info2->toaxis = 0;

// 	printTransform("init", init->e);
// 	printTransform(init->paxis->name, init->paxis->a);

  	gptaddEBelement( init, catleg_sim, catleg_exit, GPTELEM_GLOBAL, info1 ) ;
// 	info1->print();



	// the init structure has the has the ECS, i.e the axis for the element and the transformation for the element relative to 
	// the axis has the transform from the wcs to itself
	// strategy: calculate element coordinates in the wcs system, then transform to the to-axis

  	gpttransform fromTransform = fromaxis->a;		// description of the from-axis in the wcs system: first offset, then rotation
  	gpttransform toTransform = toaxis->a;			// description of the to-axis in the wcs system: first offset, then rotation
  	gpttransform ecsTransform = init->e;			// description of the magnet in the from-axis system
//   	printTransform("element", ecsTransform);


  
  	// calc offset of dipole in the toaxis coordinate system
  	double m[3];
  	for(int i=0; i<3; i++)
  	{
  	 	m[i] = fromTransform.o[i] - toTransform.o[i];	// 
  		for(int j=0; j<3; j++)
	  		m[i] += fromTransform.m[i][j] * ecsTransform.o[j];
  	}
//   	printf("m = %e %e %e\n", m[0], m[1], m[2]);

  	double p[3];
  	for(int i=0; i<3; i++)
  	{
		p[i] =0.;
  		for(int j=0; j<3; j++)
		{
	  		p[i] += toTransform.m[j][i] * m[j];
// 	  		printf("%d %d  %e %e %e \n", i, j, p[i] , toTransform.m[j][i] , m[j]);
		}
  	}
//   	printf("p = %e %e %e\n", p[0], p[1], p[2]);
  	for(int i=0; i<3; i++)
  	{
		init->e.o[i] =p[i];
// 		printf("i = %d %e\n", i, init->e.o[i]);
  	}
//   	printf("m = %e %e %e\n", init->e.o[0], init->e.o[1], init->e.o[2]);

	// cacl matrix  of dipole in the toaxis coordinate system
  	double mm[3][3];
  	for(int i=0; i<3; i++)
  	for(int j=0; j<3; j++)
  	{
		mm[i][j] =0.;
  		for(int k=0; k<3; k++)
		{
			mm[i][j] += fromTransform.m[i][k] * ecsTransform.m[k][j];
		}
  	}


  	for(int i=0; i<3; i++)
  	for(int j=0; j<3; j++)
  	{
		init->e.m[i][j] =0.;
  		for(int k=0; k<3; k++)
		{
			init->e.m[i][j] += toTransform.m[k][i] * mm[k][j];
		}
  	}

//   	printTransform("2nd", init->e);



  	init->paxis = info1->toaxis;	// the axis where to put this element
// 	printTransform("init", init->e);
// 	printTransform(init->paxis->name, init->paxis->a);
  	gptaddEBelement( init, catleg_sim, catleg_exit, GPTELEM_GLOBAL, info2 ) ;
// 	info2->print();

//   	printf("done\n");
}


static int catleg_sim(gptpar *par, double t, struct catleg_info *info)
{
	if( fabs(Z) > info->length) return 0;
	if(info->toaxis  && Z > 0) par->newaxis = info->toaxis ;

	double z2 = Z - info->z1;
	double x2 = X - info->x1;

	double z3 =  z2* info->c1 + x2*  info->s1 - info->dl;

	double f = info->b1*z3 + info->b2*(z3*z3-Y*Y);
	double h = Y*(info->b1 + 2. * info->b2 * z3);

	double ef  = exp(f);
	double sh  = sin(h);
	double ch  = cos(h);
	double efch = ef*ch;
	double efsh = ef*sh;
	double nom = info->B / (1. + 2.* efch + ef*ef);

	double bz =   nom * efsh;
	double by =  -nom * (1.+efch);

	BY =   by;
	BX = - info->s1 * bz;
	BZ =   info->c1 * bz;


// 	double winkel = atan2(par->GBr[0],par->GBr[2])*180./M_PI;
//   	fprintf(info->fdiag, "%15.7e   %15.7e %15.7e %15.7e    %15.7e %15.7e %15.7e %s\n", t, par->r[0], par->r[1], par->r[2], BX, BY, BZ,  info->name);
//   	fprintf(info->fdiag, "%d  %15.7e   %15.7e %15.7e %15.7e       %15.7e %15.7e %15.7e     %15.7e %15.7e %15.7e   %15.7e    %15.7e %15.7e %15.7e %s\n",
// 			par->ID, t, X, Y, Z, par->Wr[0], par->Wr[1], par->Wr[2], par->GBr[0], par->GBr[1], par->GBr[2], winkel,  BX, BY, BZ,  info->name);
// 	
// 	EX = EY = EZ = 0.;
  	return 1;
}

static void catleg_exit(struct catleg_info *info)
{
}


