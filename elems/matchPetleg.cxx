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
#include <pthread.h>


#define gptmalloc malloc

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "Solver.hxx"
#include "tools.hxx"
#include "ObjectiveFunction.hxx"

#ifdef WIN32
#include <crtdbg.h>
#endif

pthread_mutex_t rtxMutex= PTHREAD_MUTEX_INITIALIZER;

inline double square(double x) { return x*x;}

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
  printf("num %d %d %d\n", info->numx, info->numy, info->numz); 

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


//   printf(" dog %20.12f %20.12f %20.12f |  %20.12e %20.12e %20.12e\n",  X, Y, Z, BX, BY, BZ);
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

// 	printf("flip	%d\n", flipx);
// 	printf("angle	%15.7e\n", angle);
// 	printf("rho  	%15.7e\n", rho  );




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
// 	info2->print();



	// we add the second half of the magnet to the toaxis CCS
  	memcpy(info1, info2, sizeof(struct catleg_info) ) ;

	info1->name[4]='1';
	info1->z1 = -info1->z1;
	info1->c1 = cos(M_PI - angle/2. + phi_in);
	info1->s1 = sin(M_PI - angle/2. + phi_in);
// 	info1->print();

//   	printf("done\n");
}

static int catleg_sim(double X, double Y, double Z, double & BX, double & BY, double & BZ)
{

  	struct catleg_info *info=  Z > 0 ? info2 : info1;
	if( fabs(Z) > info->length) return 0;
// 	info->print();

// 	printf("z1,x1 %15.4f %15.4f \n", info->z1,info->x1);
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


//   	printf(" cat %20.12f %20.12f %20.12f |  %20.12e %20.12e %20.12e\n",  X, Y, Z, BX, BY, BZ);

  	return 1;
}

static void catleg_exit(struct catleg_info *info)
{
}


class M2 : public UnconstrainedObjectiveFunction
{
	double catfield;
	double catdl;
	double catb1; 
	double catb2;
	double dcatfield;
	double dcatdl;
	double dcatb1; 
	double dcatb2;
	double bestgoal;
	int nrun;
  	int numx;
  	int numy;
  	int numz;
  	double xdelta;
  	double ydelta;
  	double zdelta;
  	double xmin;
public:
	M2(double catfield, double catdl, double catb1, double catb2);
    	double eval(Vector v, int *nerror, int resource);
};

M2::M2( double catfield, double catdl, double catb1, double catb2)
{
	this->catfield =    catfield  ;
	this->catdl    =    catdl;
	this->catb1    =    catb1;
	this->catb2    =    catb2;
	dcatfield = 1e-6;
	dcatdl = 1e-2;
	dcatb1 = 5.; 
	dcatb2 =10.;
	xOptimal.setSize(4); 	// sets the number of variables
	xStart.setSize(4);
	bestgoal = 1e33;
	nrun = 0;
  	numx = infos.numx;
  	numy = infos.numy;
  	numz = infos.numz;
  	xdelta = infos.xdelta;
  	ydelta = infos.ydelta;
  	zdelta = infos.zdelta;
  	xmin = - (numx-1)/2 * xdelta;
}

int main(int argc, char ** argv)
{




    	double rhoStart=1e-0, rhoEnd=1e-5;
    	int niter=100000;

	double dogfield =    3.6900549496e-05  ;
	dogleg_init(0, dogfield);

	double catfield =    -9.249261732096e-03  ;
	double catdl =  0.;
	double catb1 = 35.; 
	double catb2 =  0.;
    	M2 *mof = new M2(catfield, catdl, catb1, catb2);
    	ObjectiveFunction *of = mof;
	CONDOR(rhoStart, rhoEnd, niter, of);
}


double M2::eval(Vector v, int *nerror, int resource)
{
	if(nerror) *nerror=0;
    	double *p=v;
    	double *start=xStart;

	double bx0cat, by0cat, bz0cat;
	double bx0dog, by0dog, bz0dog;
	double grad, curly, curlz, curlx;

	const double w20 = 20.*M_PI/180.;
	double tcatfield = catfield + dcatfield*p[0];
	double tcatdl    = catdl + dcatdl*p[1];
	double tcatb1    = catb1 + dcatb1*p[2];
	double tcatb2    = catb2 + dcatb2*p[3];
	catleg_init( 0.3, w20,  tcatfield, w20/4., w20/4., tcatdl, tcatb1, tcatb2, 0);
// 	printf("init\n");

	double sumx =0.;
	double sumy =0.;
	double sumz =0.;
	int n=0;
	for(int i=73; i< 127; i++)
	for(int j=0; j< numy; j++)
	for(int k=0; k< numz; k++)
	{
		double X; double Y; double Z;
		X = xmin + i*xdelta;
		Y = j*ydelta;
		Z = k*zdelta;
		catleg_sim( X, Y, Z, bx0cat, by0cat, bz0cat);
		dogleg_sim( X, Y, Z, bx0dog, by0dog, bz0dog);
		sumx += square(bx0cat - bx0dog);
		sumy += square(by0cat - by0dog);
		sumz += square(bz0cat - bz0dog);
		n++;
	}
	double sum = sumx + sumy + sumz;
	printf("%8d   sum = %15.7e %15.7e %15.7e %15.7e : %15.7e\n", nrun++, tcatfield, tcatdl, tcatb1, tcatb2, sum/n);
	return sum;

}
