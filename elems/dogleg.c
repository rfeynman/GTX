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


struct dogleg_info
{

	char name[20];	// info name
  /* Array size information */
  float  xdelta, ydelta, zdelta, xmin;
  int numx, numy, numz ;

  /* Field information */
  float *fx, *fy, *fz;

  int flipx;	// if true flip the magnet around the longitudinal axis by 180 degree
  int * refCount;

  struct axis  *fromaxis, *toaxis ;

  FILE * fdiag;
} ;

// void printTransform(const char * name, gpttransform trm)
// {
//   printf("axis %s:\n", name);
//   for(int i=0; i<3; i++)
//   {
// 	  printf("%f  %f  %f      %f\n", trm.m[i][0], trm.m[i][1], trm.m[i][2], trm.o[i]);
//   }
// }

static int dogleg_sim(gptpar *par, double t, struct dogleg_info *info) ;
static void dogleg_exit(struct dogleg_info *info) ;
static int numDog=0;



//	dogleg(ECS, file)
void dogleg_init(gptinit *init)
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

  const char *filename ;

  const char  *toaxisname ;

  gptbuildECS( init ) ;

  /* Commandline parsing */
  numc=gptgetargnum(init) ;
  if( numc!=3 && numc!=4 )
    gpterror( "Syntax: %s(ECS,mapfile, scalefactor, flipx, [toAxis])\n", gptgetname(init) ) ;

// dogleg("wcs",  0.013057, 0.000000, 38.938675,      9.848078e-01,  0.000000e+00,  -1.736482e-01,     0.000000e+00,  1.000000e+00,  0.000000e+00, "20_deg_dipole_map.bin",   -dogfield, 1, "ccs1");

  /* Initialize GDF */
//   printf("numc %d\n", numc);
  filename = gptgetargstring(init,1) ;
//   printf("file %s\n", filename);
  ffac = gptgetargdouble(init,2) ;
//   printf("factor %e\n", ffac);

  info = (struct dogleg_info *)gptmalloc( sizeof(struct dogleg_info) ) ;
  strcpy(info->name, "info1");

  info->flipx = gptgetargint(init,3) ;
//   printf("flip %d\n", info->flipx);
  info->fromaxis = init->paxis;
  if(numc == 4)
  {
  	toaxisname  =gptgetargstring(init,4) ;
//   	printf("new axis %s\n", toaxisname);
  	if( (  info->toaxis=getaxis(  toaxisname))==NULL )
	      	gpterror( "Specified axis \"%s\" not found\n", toaxisname ) ; 
  }
  else info->toaxis=0;

  char diagFile[200];
  sprintf(diagFile, "dog1Diag%d", numDog);
  info->fdiag = fopen(diagFile, "w");

  fheader = (float *) iheader;
  fp = open(filename, O_RDONLY);
  if( fp < 0 )
  {
	perror(filename);
  	gpterror( "cant open %s \n", filename);
  }
  nb = read(fp, iheader, 32);
  if(32 != nb) gpterror( "first read only %ld bytes, must be 32\n", nb);
  if(10071951 != iheader[0]) gpterror( "wrong magic number in dogleg file\n");
  info->numx = iheader[2] ; 
  info->numy = iheader[3] ; 
  info->numz = iheader[4] ; 
  info->xdelta = fheader[5] ; 
  info->ydelta = fheader[6] ; 
  info->zdelta = fheader[7] ;
  info->xmin = (info->numx-1)/2 * info->xdelta;

  points = iheader[2]*iheader[3]*iheader[4];
  nbs = 3 * points  * sizeof(float);  /* number of bytes for all 3 field components */
  info->fx = ((float *)gptmalloc( nbs )) ;
  info->refCount = (int *) (info->fx+1);
  *(info->refCount) = 1;
  nb = read(fp, info->fx, nbs );
  if(nbs != nb) gpterror( "second  read only %ld bytes, must be %d\n", nb, nbs);
  for(int i=0; i< 3 * points; i++) info->fx[i] *= ffac;
  info->fy = info->fx + points;
  info->fz = info->fy + points;
  close(fp);
  gptaddEBelement( init, dogleg_sim, dogleg_exit, GPTELEM_GLOBAL, info ) ;


  if(numc == 4)
  {
  	info2 = (struct dogleg_info *)gptmalloc( sizeof(struct dogleg_info) ) ;
  	memcpy(info2, info, sizeof(struct dogleg_info) ) ;
  	strcpy(info->name, "info2");
  	*(info2->refCount) = 2;
  	gpttransform fromTransform = info->fromaxis->a;
//   	printTransform(info->fromaxis->name, fromTransform);

  	gpttransform ecsTransform = init->e;
//   	printTransform("element", ecsTransform);

  	gpttransform toTransform = info->toaxis->a;
//   	printTransform(info->toaxis->name, toTransform);


  	init->paxis = info->toaxis;
  
  	// calc offset of dipole in the toaxis coordinate system
  	double m[3];
  	for(int i=0; i<3; i++)
  	{
  	 	m[i] = fromTransform.o[i] - toTransform.o[i];
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
//   	printTransform(info->toaxis->name, toTransform);

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


  	sprintf(diagFile, "dog2Diag%d", numDog);
  	info2->fdiag = fopen(diagFile, "w");


  	gptaddEBelement( init, dogleg_sim, dogleg_exit, GPTELEM_GLOBAL, info2 ) ;

//   	printf("done\n");
  }
  numDog++;
}


static int dogleg_sim(gptpar *par, double t, struct dogleg_info *info)
{
  int i,j,k ;
  double fu,fv,fw, gu,gv,gw ;
  int N000, N001, N010, N011, N100, N101, N110, N111;
  double F000, F001, F010, F011, F100, F101, F110, F111;
  int numx, numy, numz;
//   printf("x z %f %f\n", X,Z);

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



//   if(XX < 0.)     //   flipping x , we ain't goona flip x
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
	if(info->toaxis && par->paxis == info->fromaxis && Z > 0)
	{
// 		printf("axs %p -> %p \n", par->paxis, info->toaxis);
// 		printf("axs %s -> %s \n", par->paxis->name, info->toaxis->name);
		par->newaxis = info->toaxis ;
// 		printf(" flipped %e %d\n", par->n, par->ID);
	}


	double winkel = atan2(par->GBr[0],par->GBr[2])*180./M_PI;
  	fprintf(info->fdiag, "%d  %15.7e   %15.7e %15.7e %15.7e       %15.7e %15.7e %15.7e     %15.7e %15.7e %15.7e   %15.7e    %15.7e %15.7e %15.7e %s\n",
			par->ID, t, X, Y, Z, par->Wr[0], par->Wr[1], par->Wr[2], par->GBr[0], par->GBr[1], par->GBr[2], winkel,  BX, BY, BZ,  info->name);
  return( 1 ) ;
}

static void dogleg_exit(struct dogleg_info *info)
{
  	(*(info->refCount))--;
	if(*(info->refCount) == 0 ) free(info->fx);		// is used 2 times add flag 
  	free(info) ;
}


