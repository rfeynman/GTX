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


double BX, BY, BZ;

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

  FILE * fdiag;
} infos;

// void printTransform(const char * name, gpttransform trm)
// {
//   printf("axis %s:\n", name);
//   for(int i=0; i<3; i++)
//   {
// 	  printf("%f  %f  %f      %f\n", trm.m[i][0], trm.m[i][1], trm.m[i][2], trm.o[i]);
//   }
// }

// static int dogleg_sim(double X, double Y, double Z);
// static void dogleg_exit(struct dogleg_info *info) ;
static int numDog=0;



//	dogleg(ECS, file)
void dogleg_init(int flipx)
{
  struct dogleg_info *info, *info2 ;

  /* Commandline parameters */
  int numc ;
  double ffac ;

  /* Grid characteristics */
  int points , fp;
  ssize_t nb, nbs;
  int iheader[8];
  float * fheader;

  const char *filename = "/home/jorg/LowEnergyCooling/GPT/uturn/20_deg_dipole_map.bin";



  /* Initialize GDF */
//   printf("numc %d\n", numc);
//   filename = gptgetargstring(init,1) ;
//   printf("file %s\n", filename);
  ffac = 1. ;
//   printf("factor %e\n", ffac);

  info =  &infos;

  info->flipx = flipx ;
//   printf("flip %d\n", info->flipx);
  info->fromaxis = 0;
  
  info->toaxis=0;

  char diagFile[200];
  sprintf(diagFile, "dogDiag%d", numDog++);
  info->fdiag = fopen(diagFile, "w");

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


static int dogleg_sim(double X, double Y, double Z)
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


  printf(" %12.4f %12.4f %12.4f |  %12.4f %12.4f %12.4f\n",  X, Y, Z, BX, BY, BZ);
  return( 1 ) ;
}

static void dogleg_exit(struct dogleg_info *info)
{
  	(*(info->refCount))--;
	if(*(info->refCount) == 0 ) free(info->fx);		// is used 2 times add flag 
  	free(info) ;
}




main(int argc, char ** argv)
{

	int flip;
	double delta, delta0, deltap, deltam;;
	scanf("%d%lf", &flip, &delta);
	delta *= 0.001;
	delta0 =0.;
	deltap =  delta;
	deltam = -delta;
	delta = 2. * delta;

	dogleg_init(flip);
	printf("init\n");
	double X; double Y; double Z;
	while(1)
	{

		scanf("%lf%lf%lf", &X, &Y, &Z);
		printf("acll\n");
		printf("null  ");
		dogleg_sim( X*0.001+delta0, Y*0.001+delta0, Z*0.001+delta0);
		double bx000 = BX; double by000 = BY; double bz000 = BZ;

		printf("m00   ");
		dogleg_sim( X*0.001+deltam, Y*0.001+delta0, Z*0.001+delta0);
		double bxm00 = BX; double bym00 = BY; double bzm00 = BZ;
		printf("p00   ");
		dogleg_sim( X*0.001+deltap, Y*0.001+delta0, Z*0.001+delta0);
		double bxp00 = BX; double byp00 = BY; double bzp00 = BZ;

		printf("0m0   ");
		dogleg_sim( X*0.001+delta0, Y*0.001+deltam, Z*0.001+delta0);
		double bx0m0 = BX; double by0m0 = BY; double bz0m0 = BZ;
		printf("0p0   ");
		dogleg_sim( X*0.001+delta0, Y*0.001+deltap, Z*0.001+delta0);
		double bx0p0 = BX; double by0p0 = BY; double bz0p0 = BZ;

		printf("00m   ");
		dogleg_sim( X*0.001+delta0, Y*0.001+delta0, Z*0.001+deltam);
		double bx00m = BX; double by00m = BY; double bz00m = BZ;
		printf("00p   ");
		dogleg_sim( X*0.001+delta0, Y*0.001+delta0, Z*0.001+deltap);
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

		double grad  =  dbxdx + dbydy + dbzdz;
		double curly =  dbxdz - dbzdx;
		double curlz =  dbydx - dbxdy; 
		double curlx =  dbzdy - dbydz; 

		printf("curl %e %e %e %e\n", grad, curly, curlz, curlx);
	}

}
