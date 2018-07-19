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


struct dogleg_info
{

  /* Array size information */
  float  xdelta, ydelta, zdelta, xmin;
  int numx, numy, numz ;

  /* Field information */
  float *fx, *fy, *fz;

  int flipx;	// if true flip the magnet around the longitudinal axis by 180 degree
  int * refCount;

  struct axis  *fromaxis, *toaxis ;

} infos;



//	dogleg(ECS, file)
void dogleg_init(int flipx, double ffac)
{
  struct dogleg_info *info, *info2 ;

  /* Commandline parameters */
  int numc ;

  /* Grid characteristics */
  int points , fp;
  ssize_t nb, nbs;
  int iheader[8];
  float * fheader;

  const char *filename = "/home/jorg/LowEnergyCooling/Dima/realSolenoidField/20_deg_dipole/20_deg_dipole_map.bin";




  info =  &infos;

  info->flipx = flipx ;
//   printf("flip %d\n", info->flipx);
  info->fromaxis = 0;
  
  info->toaxis=0;


  fheader = (float *) iheader;
  fp = open(filename, O_RDONLY);
  if( fp < 0 )
  {
	perror(filename);
	exit(-1);
  }
  nb = read(fp, iheader, 32);
  if(32 != nb) printf( "first read only %ld bytes, must be 32\n", nb);
  if(10071951 != iheader[0]) printf( "wrong magic number in dogleg file\n");
  info->numx = iheader[2] ; 
  info->numy = iheader[3] ; 
  info->numz = iheader[4] ; 
  info->xdelta = fheader[5] ; 
  info->ydelta = fheader[6] ; 
  info->zdelta = fheader[7] ;
  info->xmin = (info->numx-1)/2 * info->xdelta;

  points = iheader[2]*iheader[3]*iheader[4];
  nbs = 3 * points  * sizeof(float);  /* number of bytes for all 3 field components */
  info->fx = ((float *)malloc( nbs )) ;
  info->refCount = (int *) (info->fx+1);
  *(info->refCount) = 1;
  nb = read(fp, info->fx, nbs );
  if(nbs != nb) printf( "second  read only %ld bytes, must be %ld\n", nb, nbs);
  for(int i=0; i< 3 * points; i++) info->fx[i] *= ffac;
  info->fy = info->fx + points;
  info->fz = info->fy + points;
  close(fp);


}


static int dogleg_sim(double X, double Y, double Z, double & BX, double & BY, double & BZ)
{
  int i,j,k ;
  double fu,fv,fw, gu,gv,gw ;
  int N000, N001, N010, N011, N100, N101, N110, N111;
  double F000, F001, F010, F011, F100, F101, F110, F111;
  int numx, numy, numz;
//   printf("x y z %f %f %f\n", X,Y, Z);
  struct dogleg_info *info, *info2 ;
  info =  &infos;

  /* calc field in the positive quadrant and use symmetry */

  double XX= X;	// need to copy so not screw up gpt
  double YY= Y;
  double ZZ= Z;

  if(info->flipx) { XX =  -XX; YY = - YY; }


  double XXX = XX;
  XXX  += info->xmin;
  double YYY = fabs(YY);
  double ZZZ = fabs(ZZ);
  /* Calculate master index and offset fractions */
  i= XXX/info->xdelta;
  j= YYY/info->ydelta;
  k= ZZZ/info->zdelta;

  numx = info->numx;
  numy = info->numy;
  numz = info->numz;

//   printf("y %12.4e %12.4e %12.4e %4d %4d %4d \n", X, Y, Z, i, j, k);

  if(  ( i  <   0    ) || ( i >= numx-1 ) || ( j >= numy-1 ) || ( k >= numz-1 ) ) return( 0 ) ;

  fu= (XXX-i*info->xdelta)/info->xdelta ;
  fv= (YYY-j*info->ydelta)/info->ydelta ;
  fw= (ZZZ-k*info->zdelta)/info->zdelta ;
  gu = 1.-fu ;
  gv = 1.-fv ;
  gw = 1.-fw ;

  /* Calculate boundary offsets */
  N000 = (i*numy+j)*numz+k;
  N001 = N000+1;
  N010 = N000+numz;
  N011 = N010+1;
  N100 = N000+numy*numz;
  N101 = N100+1;
  N110 = N100+numz;
  N111 = N110+1;

  F000 = gu*gv*gw ;
  F001 = gu*gv*fw ;
  F010 = gu*fv*gw ;
  F011 = gu*fv*fw ;
  F100 = fu*gv*gw ;
  F101 = fu*gv*fw ;
  F110 = fu*fv*gw ;
  F111 = fu*fv*fw ;

  BX = F000*info->fx[N000] + F001*info->fx[N001] + F010*info->fx[N010] + F011*info->fx[N011] + F100*info->fx[N100] + F101*info->fx[N101] + F110*info->fx[N110] + F111*info->fx[N111] ;
  BY = F000*info->fy[N000] + F001*info->fy[N001] + F010*info->fy[N010] + F011*info->fy[N011] + F100*info->fy[N100] + F101*info->fy[N101] + F110*info->fy[N110] + F111*info->fy[N111] ;
  BZ = F000*info->fz[N000] + F001*info->fz[N001] + F010*info->fz[N010] + F011*info->fz[N011] + F100*info->fz[N100] + F101*info->fz[N101] + F110*info->fz[N110] + F111*info->fz[N111] ;


//   if(X < 0.)     //   flipping x , we ain't goona flip x
//   {
//           //  By constant, dBy/dy constant, Bz constant, dBz/dz constant  => dBx/dx constant => Bx changes
//           BX = -BX ;
//   }
  if(YY < 0.)       //  flipping y
  {
          //  By constant, dBy/dy changes, Bz changes, dBz/dz changes  => dBx/dx changes => Bx changes
          BX = -BX ;
          BZ = -BZ ;
  }
  if(ZZ < 0.)      //  flipping z
  {
          //  By constant, dBy/dy constant, Bz changes, dBz/dz constant  => dBx/dx constant => Bx constant
          BZ = -BZ;
  }
//   printf("x %12.4f %12.4e %12.4e %4d %4d %4d %12.4e\n", X, Y, Z, i, j, k, BY);



  if(info->flipx)
  {
	  BX = -BX;
	  BY = -BY;
  }


/*      flip to the other axis  */
// 	printf("axs %s -> %s  %e\n", par->paxis->name, info->toaxis->name, X);


  printf(" dog %20.12f %20.12f %20.12f |  %20.12e %20.12e %20.12e\n",  X, Y, Z, BX, BY, BZ);
  return( 1 ) ;
}

static void dogleg_exit(struct dogleg_info *info)
{
  	(*(info->refCount))--;
	if(*(info->refCount) == 0 ) free(info->fx);		// is used 2 times add flag 
  	free(info) ;
}


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

	printf("flip	%d\n", flipx);
	printf("angle	%15.7e\n", angle);
	printf("rho  	%15.7e\n", rho  );




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

//   	printf("done\n");
}

static int catleg_sim(double X, double Y, double Z, double & BX, double & BY, double & BZ)
{

  	struct catleg_info *info=  Z > 0 ? info2 : info1;
	if( fabs(Z) > info->length) return 0;
// 	info->print();

	printf("z1,x1 %15.4f %15.4f \n", info->z1,info->x1);
	double z2 = Z - info->z1;
	double x2 = X - info->x1;
// 	printf("z2,x2 %15.4f %15.4f \n", z2,x2);

	double z3 =  z2* info->c1 + x2*  info->s1 - info->dl;
// 	printf("z3,x3 %15.4f %15.4f \n", z2,0.);

	double f = info->b1*z3 + info->b2*(z3*z3-Y*Y);
	double h = Y*(info->b1 + 2. * info->b2 * z3);
// 	printf("b1,b2 %15.4f %15.4f \n",info->b1, info->b2);

	double ef  = exp(f);
	double sh  = sin(h);
	double ch  = cos(h);
	double efch = ef*ch;
	double efsh = ef*sh;
	double nom = 1. + 2.* efch + ef*ef;
// 	printf("ef,nom %15.7e %15.7e\n", ef,  nom);
// 	printf("sh,ch %15.4f %15.4f \n", sh,ch);

	double bz =   info->B * efsh      / nom;
	double by =  -info->B * (1.+efch) / nom;

// 	printf("bz,by %15.7e %15.7e\n", bz,by);
	BY = by;
	BX =  -info->s1 * bz;
	BZ =   info->c1 * bz;
// 	printf("BZ,BY,BX %15.7e %15.7e %15.7e \n", BZ,BY,BX);


  	printf(" cat %20.12f %20.12f %20.12f |  %20.12e %20.12e %20.12e\n",  X, Y, Z, BX, BY, BZ);

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

	double dogfield =    3.6900549496e-05  ;
	double catfield =    -9.269261732096e-03  ;
	double catdl=0.;
	double catb1 = 35; 
	double catb2 = 0;
	double w20 = 20.*M_PI/180.;

	double bx0cat, by0cat, bz0cat;
	double bx0dog, by0dog, bz0dog;
	double bx000, by000, bz000, bxm00, bym00, bzm00, bxp00, byp00, bzp00, bx0m0, by0m0, bz0m0, bx0p0, by0p0, bz0p0, bx00m, by00m, bz00m, bx00p, by00p, bz00p;
	double grad, curly, curlz, curlx;

	catleg_init( 0.3, w20,  catfield, w20/4., w20/4., catdl, catb1, catb2, flip);
	catleg_sim( 0., 0., 0., bx0cat, by0cat, bz0cat);
	dogleg_init(flip, dogfield);
	dogleg_sim( 0., 0., 0., bx0dog, by0dog, bz0dog);
	printf("init\n");

	double X; double Y; double Z;
	while(1)
	{
		printf("X,Y,Z? ");
		scanf("%lf%lf%lf", &X, &Y, &Z);
		X *= 1e-3;
		Y *= 1e-3;
		Z *= 1e-3;

		catleg_sim( X+delta0, Y+delta0, Z+delta0, bx000, by000, bz000);

// 		catleg_sim( X+deltam, Y+delta0, Z+delta0, bxm00, bym00, bzm00);
// 		catleg_sim( X+deltap, Y+delta0, Z+delta0, bxp00, byp00, bzp00);
// 
// 		catleg_sim( X+delta0, Y+deltam, Z+delta0, bx0m0, by0m0, bz0m0);
// 		catleg_sim( X+delta0, Y+deltap, Z+delta0, bx0p0, by0p0, bz0p0);
// 
// 		catleg_sim( X+delta0, Y+delta0, Z+deltam, bx00m, by00m, bz00m);
// 		catleg_sim( X+delta0, Y+delta0, Z+deltap, bx00p, by00p, bz00p);
// 
// 
// 		grad  =   ( ( bxp00 - bxm00) + ( by0p0 - by0m0) + ( bz00p - bz00m)  )/ by0cat;
// 		curly =   ( ( bx00p - bx00m) - ( bzp00 - bzm00)                     )/ by0cat;
// 		curlz =   ( ( byp00 - bym00) - ( bx0p0 - bx0m0)                     )/ by0cat; 
// 		curlx =   ( ( bz0p0 - bz0m0) - ( by00p - by00m)                     )/ by0cat; 
// 		printf("grad %e curly %e curlz %e curlx %e\n", grad, curly, curlz, curlx);



		dogleg_sim( X+delta0, Y+delta0, Z+delta0, bx000, by000, bz000);

// 		dogleg_sim( X+deltam, Y+delta0, Z+delta0, bxm00, bym00, bzm00);
// 		dogleg_sim( X+deltap, Y+delta0, Z+delta0, bxp00, byp00, bzp00);
// 
// 		dogleg_sim( X+delta0, Y+deltam, Z+delta0, bx0m0, by0m0, bz0m0);
// 		dogleg_sim( X+delta0, Y+deltap, Z+delta0, bx0p0, by0p0, bz0p0);
// 
// 		dogleg_sim( X+delta0, Y+delta0, Z+deltam, bx00m, by00m, bz00m);
// 		dogleg_sim( X+delta0, Y+delta0, Z+deltap, bx00p, by00p, bz00p);
// 
// 
// 		grad  =   ( ( bxp00 - bxm00) + ( by0p0 - by0m0) + ( bz00p - bz00m)  )/ by0dog;
// 		curly =   ( ( bx00p - bx00m) - ( bzp00 - bzm00)                     )/ by0dog;
// 		curlz =   ( ( byp00 - bym00) - ( bx0p0 - bx0m0)                     )/ by0dog; 
// 		curlx =   ( ( bz0p0 - bz0m0) - ( by00p - by00m)                     )/ by0dog; 
// 		printf("grad %e curly %e curlz %e curlx %e\n", grad, curly, curlz, curlx);
	}

}
