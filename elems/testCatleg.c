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

#define gptmalloc malloc

struct catleg_info
{
	double x1, z1;	// points of arc on the pole faces
	double c1, s1;	// rotation from the magnet coordinate system to the pole face system
  	double B ;	      /* Specified B field */
  	double dl ;         /* Offset for edge (difference between real and e2fective length) */
  	double b1, b2 ;     /* Parameters of Enge fringe field */ 
  	double length;

  	FILE * fdiag;
	void print()
	{
		printf("%s\n", p);
		printf("x1     = %20.12e\n", x1);
		printf("z1     = %20.12e\n", z1);
		printf("c1     = %20.12e\n", c1);
		printf("s1     = %20.12e\n", s1);
		printf("B      = %20.12e\n", B );
		printf("dl     = %20.12e\n", dl);
		printf("b1     = %20.12e\n", b1);
		printf("b2     = %20.12e\n", b2);
		printf("length = %20.12e\n", length);
	}
} ;

  	struct catleg_info *info1= (struct catleg_info *)gptmalloc( sizeof(struct catleg_info) ) ;
  	struct catleg_info *info2= (struct catleg_info *)gptmalloc( sizeof(struct catleg_info) ) ;

// void printTransform(const char * name, gpttransform trm)
// {
//   printf("axis %s:\n", name);
//   for(int i=0; i<3; i++)
//   {
// 	  printf("%f  %f  %f      %f\n", trm.m[i][0], trm.m[i][1], trm.m[i][2], trm.o[i]);
//   }
// }

static int numcat=0;

// dogleg("wcs",  0.013057, 0.000000, 38.938675,      9.848078e-01,  0.000000e+00,  -1.736482e-01,     0.000000e+00,  1.000000e+00,  0.000000e+00, "20_deg_dipole_map.bin",   -dogfield, 1, "ccs1");
// catleg("wcs",  0.013057, 0.000000, 38.938675,      9.848078e-01,  0.000000e+00,  -1.736482e-01,     0.000000e+00,  1.000000e+00,  0.000000e+00, "ccs1", 0.3, w20,  dogfield, 20.4, w20/4, 0., 0.1, 0., 1);



//	catleg(ECS, file)
void catleg_init(double length, double angle, double B, double phi_in, double phi_out, double dl, double b1, double b2, int flipx ) 
{


/* Fill coe2ficients for local use. */
	double rho = length/angle;	// deflection radius


/* Get angles of pole-faces */
  	if( b1<0 ) printf( "car: b1 must be nonnegative\n" ) ;
	if(flipx) { B = -B; angle = -angle; phi_in = -phi_in; phi_out = -phi_out; }




  	char diagFile[200];
  	sprintf(diagFile, "catDiag%d", numcat++);
  	info1->fdiag = fopen(diagFile, "w");

  	info1->B = B;	      /* Specified B field */
  	info1->dl = dl;         /* Offset for edge (difference between real and e2fective length) */
  	info1->b1 = b1 ;     /* Parameters of Enge fringe field */ 
  	info1->b2 = b2 ;     /* Parameters of Enge fringe field */ 
	info1->length = length;

	double ca = cos(angle/2.);
	double sa = sin(angle/2.);
	info1->x1 = rho* (1. - ca);
	info1->z1 = rho* sa;
	info1->c1 = cos(M_PI - angle/2. + phi_in);
	info1->s1 = sin(M_PI - angle/2. + phi_in);



	// we add the second half of the magnet to the toaxis CCS
  	memcpy(info2, info1, sizeof(struct catleg_info) ) ;

	info2->z1 = -info1->z1;
	info2->c1 = cos(angle/2. - phi_out);
	info2->s1 = sin(angle/2. - phi_out);

	info1->print("info1");
	info2->print("info2");

//   	printf("done\n");
}

double BX; double BY; double BZ;
static int catleg_sim(double X, double Y, double Z)
{

  	struct catleg_info *info=  Z > 0 ? info2 : info1;
	if( fabs(Z) > info->length) return 0;

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
	double nom = 1. + 2.* efch + ef*ef;

	double bz =   info->B * efsh      / nom;
	double by =  -info->B * (1.+efch) / nom;

	BY = by;
	BX =   info->s1 * bz;
	BZ =   info->c1 * bz;


  printf(" %20.12f %20.12f %20.12f |  %20.12e %20.12e %20.12e\n",  X, Y, Z, BX, BY, BZ);

// 	info = info2;
// 
// 	z2 = Z - info->z1;
// 	x2 = X - info->x1;
// 
// 	z3 =  z2* info->c1 + x2*  info->s1 - info->dl;
// 
// 	f = info->b1*z3 + info->b2*(z3*z3-Y*Y);
// 	h = Y*(info->b1 + 2. * info->b2 * z3);
// 
// 	ef  = exp(f);
// 	sh  = sin(h);
// 	ch  = cos(h);
// 	efch = ef*ch;
// 	efsh = ef*sh;
// 	nom = 1. + 2.* efch + ef*ef;
// 
// 	bz =   info->B * efsh      / nom;
// 	by =  -info->B * (1.+efch) / nom;
// 
// 	BY = by;
// 	BX =   info->s1 * bz;
// 	BZ =   info->c1 * bz;
// 
// 
//   printf("z==0   %20.12f %20.12f %20.12f |  %20.12e %20.12e %20.12e\n",  X, Y, Z, BX, BY, BZ);
// 
  	return 1;
}

static void catleg_exit(struct catleg_info *info)
{
}




main(int argc, char ** argv)
{

	int flip;
	double delta, delta0, deltap, deltam;;
	printf("flip, delta? ");
	scanf("%d%lf", &flip, &delta);
	delta *= 0.001;
	delta0 =0.;
	deltap =  delta;
	deltam = -delta;
	delta = 2. * delta;

	double dogfield =    -0.00926  ;
	double dogdl=0.;
	double dogb1 = 35; 
	double dogb2 = 0;
	double w20 = 20.*M_PI/180.;

	catleg_init( 0.3, w20,  dogfield, w20/4, w20/4, dogdl, dogb1, dogb2, 1);
	catleg_sim( 0., 0., 0.);
	double bx0 = BX; double by0 = BY; double bz0 = BZ;

	printf("init\n");
	double X; double Y; double Z;
	while(1)
	{
		printf("X,Y,Z? ");

		scanf("%lf%lf%lf", &X, &Y, &Z);
		printf("null  ");
		catleg_sim( X*0.001+delta0, Y*0.001+delta0, Z*0.001+delta0);
		double bx000 = BX; double by000 = BY; double bz000 = BZ;

		printf("m00   ");
		catleg_sim( X*0.001+deltam, Y*0.001+delta0, Z*0.001+delta0);
		double bxm00 = BX; double bym00 = BY; double bzm00 = BZ;
		printf("p00   ");
		catleg_sim( X*0.001+deltap, Y*0.001+delta0, Z*0.001+delta0);
		double bxp00 = BX; double byp00 = BY; double bzp00 = BZ;

		printf("0m0   ");
		catleg_sim( X*0.001+delta0, Y*0.001+deltam, Z*0.001+delta0);
		double bx0m0 = BX; double by0m0 = BY; double bz0m0 = BZ;
		printf("0p0   ");
		catleg_sim( X*0.001+delta0, Y*0.001+deltap, Z*0.001+delta0);
		double bx0p0 = BX; double by0p0 = BY; double bz0p0 = BZ;

		printf("00m   ");
		catleg_sim( X*0.001+delta0, Y*0.001+delta0, Z*0.001+deltam);
		double bx00m = BX; double by00m = BY; double bz00m = BZ;
		printf("00p   ");
		catleg_sim( X*0.001+delta0, Y*0.001+delta0, Z*0.001+deltap);
		double bx00p = BX; double by00p = BY; double bz00p = BZ;


		double dbxdx = ( bxp00 - bxm00)/ delta;
		double dbydx = ( byp00 - bym00)/ delta;
		double dbzdx = ( bzp00 - bzm00)/ delta;
		double dbxdy = ( bx0p0 - bx0m0)/ delta;
		double dbydy = ( by0p0 - by0m0)/ delta;
		double dbzdy = ( bz0p0 - bz0m0)/ delta;
		double dbxdz = ( bx00p - bx00m)/ delta;
		double dbydz = ( by00p - by00m)/ delta;
		double dbzdz = ( bz00p - bz00m)/ delta;

		double grad  =  delta / by0 *( dbxdx + dbydy + dbzdz);
		double curly =  delta / by0 *( dbxdz - dbzdx);
		double curlz =  delta / by0 *( dbydx - dbxdy); 
		double curlx =  delta / by0 *( dbzdy - dbydz); 

		printf("grad %e curly %e curlz %e curlx %e\n", grad, curly, curlz, curlx);
	}

}
