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
	char name[20];
	double x1, z1;	// points of arc on the pole faces
	double c1, s1;	// rotation from the magnet coordinate system to the pole face system
  	double B ;	      /* Specified B field */
  	double dl ;         /* Offset for edge (difference between real and e2fective length) */
  	double b1, b2 ;     /* Parameters of Enge fringe field */ 
  	double length;

	void print()
	{
		printf("%s\n", name);
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



//	catleg(ECS, file)
void catleg_init(double length, double angle, double B, double phi_in, double phi_out, double dl, double b1, double b2, int flipx ) 
{
  	if( b1<0 ) printf( "car: b1 must be nonnegative\n" ) ;
	if(flipx) { B = -B; angle = -angle; phi_in = -phi_in; phi_out = -phi_out; }
	double rho = length/angle;	// deflection radius


	// info2 for Z> 0
  	info2->B = B;	      /* Specified B field */
  	info2->dl = dl;         /* Offset for edge (difference between real and e2fective length) */
  	info2->b1 = b1 ;     /* Parameters of Enge fringe field */ 
  	info2->b2 = b2 ;     /* Parameters of Enge fringe field */ 
	info2->length = length;
	strcpy(info2->name, "info2");

	double ca = cos(angle/2.);
	double sa = sin(angle/2.);
	info2->x1 = rho* (1. - ca);
	info2->z1 = rho* sa;
	info2->c1 = cos(angle/2. - phi_out);
	info2->s1 = sin(angle/2. - phi_out);
	info2->print();



	// we add the second half of the magnet to the toaxis CCS
  	memcpy(info1, info2, sizeof(struct catleg_info) ) ;

	info1->name[4]='1';
	info1->z1 = -info1->z1;
	info1->c1 = cos(M_PI - angle/2. + phi_in);
	info1->s1 = sin(M_PI - angle/2. + phi_in);


	info1->print();
	info2->print();

//   	printf("done\n");
}

static int catleg_sim(double X, double Y, double Z, double & BX, double & BY, double & BZ)
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

	double bx0, by0, bz0, bx000, by000, bz000, bxm00, bym00, bzm00, bxp00, byp00, bzp00, bx0m0, by0m0, bz0m0, bx0p0, by0p0, bz0p0, bx00m, by00m, bz00m, bx00p, by00p, bz00p;

	catleg_init( 0.3, w20,  dogfield, w20/4, w20/4, dogdl, dogb1, dogb2, flip);
	catleg_sim( 0., 0., 0., bx0, by0, bz0);
	printf("init\n");

	double X; double Y; double Z;
	while(1)
	{
		printf("X,Y,Z? ");
		scanf("%lf%lf%lf", &X, &Y, &Z);

		printf("null  "); catleg_sim( X*0.001+delta0, Y*0.001+delta0, Z*0.001+delta0, bx000, by000, bz000);

		printf("m00   "); catleg_sim( X*0.001+deltam, Y*0.001+delta0, Z*0.001+delta0, bxm00, bym00, bzm00);
		printf("p00   "); catleg_sim( X*0.001+deltap, Y*0.001+delta0, Z*0.001+delta0, bxp00, byp00, bzp00);

		printf("0m0   "); catleg_sim( X*0.001+delta0, Y*0.001+deltam, Z*0.001+delta0, bx0m0, by0m0, bz0m0);
		printf("0p0   "); catleg_sim( X*0.001+delta0, Y*0.001+deltap, Z*0.001+delta0, bx0p0, by0p0, bz0p0);

		printf("00m   "); catleg_sim( X*0.001+delta0, Y*0.001+delta0, Z*0.001+deltam, bx00m, by00m, bz00m);
		printf("00p   "); catleg_sim( X*0.001+delta0, Y*0.001+delta0, Z*0.001+deltap, bx00p, by00p, bz00p);


		double grad  =  1. / by0 *( ( bxp00 - bxm00) + ( by0p0 - by0m0) + ( bz00p - bz00m));
		double curly =  1. / by0 *( ( bx00p - bx00m) - ( bzp00 - bzm00));
		double curlz =  1. / by0 *( ( byp00 - bym00) - ( bx0p0 - bx0m0)); 
		double curlx =  1. / by0 *( ( bz0p0 - bz0m0) - ( by00p - by00m)); 
		printf("grad %e curly %e curlz %e curlx %e\n", grad, curly, curlz, curlx);
	}

}
