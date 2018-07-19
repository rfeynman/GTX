/* map3D3Sym_B.c - 3D rectangular field map for the magnetic field */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 * modified for a magnet that is symmetric in all 3 planes. Only 1/8 of the field table is given, Jorg Kewisch, BNL,2015
 */

// #include <algorithm>
// using namespace std ;

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


struct csol_info
{
  /* GDF information */
  struct gdfmem gm ;

  /* Array size information */
  float  xdelta, ydelta, zdelta, rdelta, xqmin, zqmin;
  int nsr, nsz, ndx, ndy, ndz, nqx, nqy, nqz; 
  float delta;

  /* Field information */
  float *allFields;

  float * Bsor;
  float * Bsoz;
  float * Bx[4];	// horizontal deflection
  float * By[4];
  float * Bz[4];

} ;

// void printTransform(const char * name, gpttransform trm)
// {
//   printf("axis %s:\n", name);
//   for(int i=0; i<3; i++)
//   {
// 	  printf("%f  %f  %f      %f\n", trm.m[i][0], trm.m[i][1], trm.m[i][2], trm.o[i]);
//   }
// }

static int csol_sim(gptpar *par, double t, struct csol_info *info) ;
static void csol_exit(struct csol_info *info) ;



//	csol(ECS, file)
void csol_init(gptinit *init)
{
  struct csol_info *info;

  /* Commandline parameters */
  int numc ;

  /* Grid characteristics */
  int  fp;
  ssize_t nb, nbs;
  int iheader[16];
  float * fheader;

  const char *filename ;

//   printf("**csol**\n");

  gptbuildECS( init ) ;

  /* Commandline parsing */
  numc=gptgetargnum(init) ;
  if( numc!=6 )
    gpterror( "Syntax: %s(ECS,mapfile, solfactor, hdipfactor, vdipfactor, quadfactor, skewfactor)\n", gptgetname(init) ) ;


  /* Initialize GDF */
//   printf("numc %d\n", numc);
  filename = gptgetargstring(init,1) ;
//   printf("file %s\n", filename);
 	double solfac  = gptgetargdouble(init,2) ;
 	double hdipfac = gptgetargdouble(init,3) ;
 	double vdipfac = gptgetargdouble(init,4) ;
 	double quadfac = gptgetargdouble(init,5) ;
 	double skewfac = gptgetargdouble(init,6) ;
//   printf("factor %e\n", ffac);



  fheader = (float *) iheader;
  fp = open(filename, O_RDONLY);
  if( fp < 0 )
  {
	perror(filename);
  	gpterror( "cant open %s \n", filename);
  }
  nb = read(fp, iheader, 64);
  if(64 != nb) gpterror( "%s: first read only %d bytes, must be 32\n", gptgetname(init), nb);
  if(10071951 != iheader[0]) gpterror( "%s: wrong magic number in csol file\n", gptgetname(init));

  	info = (struct csol_info *) gptmalloc( sizeof(struct csol_info ) ) ;
	info->nsr    = iheader[2] ; 
	info->nsz    = iheader[3] ; 
	info->ndx    = iheader[4] ; 
	info->ndy    = iheader[5] ; 
	info->ndz    = iheader[6] ; 
	info->nqx    = iheader[7] ; 
	info->nqy    = iheader[8] ; 
	info->nqz    = iheader[9] ; 
	info->xdelta = fheader[11];
	info->ydelta = fheader[12];
	info->zdelta = fheader[13];
	info->rdelta = fheader[14];

// 	for(int i=0;i<10; i++)
// 		printf("if %d = %d\n", i,iheader[i]);
// 	for(int i=11;i<15; i++)
// 		printf("df %d = %e\n", i,fheader[i]);
	info->xqmin = (info->nqx-1)/2 * info->xdelta;
	info->zqmin = (info->nqz-1)/2 * info->zdelta;

  // number of data points per field table	
  int spoints = iheader[2]*iheader[3];
  int dpoints = iheader[4]*iheader[5]*iheader[6];
  int qpoints = iheader[7]*iheader[8]*iheader[9];
//   printf(" points %d %d %d \n", spoints, dpoints, qpoints);


  // We have 1 solenoid with Br and Bz2 (1x2=2)
  // We have 2 solenoid with Bx and By and Bz (2x3=6)
  // We have 2 solenoid with Bx and By and Bz (2x3=6)
  nbs = (2*spoints +6*dpoints +6*qpoints) * sizeof(float);  /* number of bytes for all 3 field components */


  info->allFields = ((float *)gptmalloc( nbs )) ;
  nb = read(fp, info->allFields, nbs );
  close(fp);
  if(nbs != nb) gpterror( "second  read only %d bytes, must be %d\n", nb, nbs);

  info->Bsor = info->allFields;
  info->Bsoz = info->Bsor   + spoints;

  info->Bx[0] = info->Bsoz  + spoints;
  info->By[0] = info->Bx[0] + dpoints;
  info->Bz[0] = info->By[0] + dpoints;

  info->Bx[1] = info->Bz[0] + dpoints;
  info->By[1] = info->Bx[1] + dpoints;
  info->Bz[1] = info->By[1] + dpoints;

  info->Bx[2] = info->Bz[1] + dpoints;
  info->By[2] = info->Bx[2] + qpoints;
  info->Bz[2] = info->By[2] + qpoints;

  info->Bx[3] = info->Bz[2] + qpoints;
  info->By[3] = info->Bx[3] + qpoints;
  info->Bz[3] = info->By[3] + qpoints;


//   info->Bsor = info->allFields;
//   info->Bsoz = info->Bsor   + spoints;

//   printf( "info->Bsor[]= %e\n", info->Bsor[0]);
//   printf( "info->Bsoz[]= %e\n", info->Bsoz[0]);
//   for(int f=0; f<4; f++)
//   {
//   	printf( "info->Bx[%d]= %e\n", f,info->Bx[f][0]);
//   	printf( "info->By[%d]= %e\n", f,info->By[f][0]);
//   	printf( "info->Bz[%d]= %e\n", f,info->Bz[f][0]);
//   }


  for(int i=0; i< 2*spoints; i++) info->Bsor[i] *= solfac;
  for(int i=0; i< 3*dpoints; i++) info->Bx[0][i] *= hdipfac;
  for(int i=0; i< 3*dpoints; i++) info->Bx[1][i] *= vdipfac;
  for(int i=0; i< 3*qpoints; i++) info->Bx[2][i] *= quadfac;
  for(int i=0; i< 3*qpoints; i++) info->Bx[3][i] *= skewfac;



  gptaddEBelement( init, csol_sim, csol_exit, GPTELEM_GLOBAL, info ) ;


}


static int csol_sim(gptpar *par, double t, struct csol_info *info)
{
  int i,j,k ;
  double fu,fv,fw, gu,gv,gw ;
  int N000, N001, N010, N011, N100, N101, N110, N111;
  double F000, F001, F010, F011, F100, F101, F110, F111;

  //	inside the solenoid?
  double ZZ = fabs(Z);
  k= ZZ/info->zdelta;
  if( k >= info->nsz-1 )
  {
//   	printf("rzxy %e %e %e %e : %e %e %e %e\n", t, X, Y, Z, 0., 0., 0., 0.);
	return( 0 ) ;
  }

  /* calc field for the solenoid  */
  double R = sqrt(X*X+Y*Y);

  i= R/info->rdelta;
  if( i >= info->nsr-1 )
  {
//   	printf("rzxy %e %e %e %e : %e %e %e %e\n", t, X, Y, Z, 0., 0., 0., 0.);
	return( 0 ) ;
  }

  fu= (R -i*info->rdelta)/info->rdelta ; gu = 1.-fu ;
  fv= (ZZ-k*info->zdelta)/info->zdelta ; gv = 1.-fv ;

  /* Calculate boundary offsets */
  N000 = i*info->nsz+k;
  N001 = N000+1;
  N010 = N000+info->nsz;
  N011 = N010+1;

  F000 = gu*gv ;
  F001 = gu*fv ;
  F010 = fu*gv ;
  F011 = fu*fv ;


  double BR;
  BR = F000*info->Bsor[N000] + F001*info->Bsor[N001] + F010*info->Bsor[N010] + F011*info->Bsor[N011];
  BZ = F000*info->Bsoz[N000] + F001*info->Bsoz[N001] + F010*info->Bsoz[N010] + F011*info->Bsoz[N011];
  if( Z > 0.) BR = - BR;
  if(R > 1e-10)
  {
  	BX = BR * X/R;
  	BY = BR * Y/R;
  }
  else
  {
  	BX = 0.;
  	BY = 0.;
  }
//   printf("rzxy %e %e %e %e : r %e z %e x %e %e\n", t, X, Y, Z, BR, BZ, BX, BY);


// ---------------   end solenoid field ----------------------------------------------------


//  calc field for the correctors
//  all corrcetors have the same Z range
  	ZZ  = Z + info->zqmin;
  	k= ZZ/info->zdelta;
  	if( k < 0 || k >= info->ndz-1 )
  {
//   	printf("rzxy %e %e %e %e : %e %e %e %e\n", t, X, Y, Z, 0., 0., 0., 0.);
	return( 0 ) ;
  }
  for(int f=0; f<4; f++)
  {
  	double XX, YY;
	int numx, numy, numz;
	if( f < 2)
	{
  		XX = fabs(X);
  		YY = fabs(Y);
		numx = info->ndx;
		numy = info->ndy;
		numz = info->ndz;
	}
	else
	{
  		XX  = X + info->xqmin;
  		YY  = Y + info->xqmin;
		numx = info->nqx;
		numy = info->nqy;
		numz = info->nqz;
	}
  	/* Calculate master index and offset fractions */
  	i= XX/info->xdelta;
  	j= YY/info->ydelta;
  	k= ZZ/info->zdelta;


//   printf("y %12.4e %12.4e %12.4e %4d %4d %4d \n", X, Y, Z, i, j, k);

  	if( i  <    0   ) continue;
  	if( j  <    0   ) continue;
  	if( i >= numx-1 ) continue;
  	if( j >= numy-1 ) continue;

  	fu= (XX-i*info->xdelta)/info->xdelta ; gu = 1.-fu ;
  	fv= (YY-j*info->ydelta)/info->ydelta ; gv = 1.-fv ;
  	fw= (ZZ-k*info->zdelta)/info->zdelta ; gw = 1.-fw ;

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

	float * fx = info->Bx[f];
	float * fy = info->By[f];
	float * fz = info->Bz[f];
  	float BBX = F000*fx[N000] + F001*fx[N001] + F010*fx[N010] + F011*fx[N011] + F100*fx[N100] + F101*fx[N101] + F110*fx[N110] + F111*fx[N111] ;
  	float BBY = F000*fy[N000] + F001*fy[N001] + F010*fy[N010] + F011*fy[N011] + F100*fy[N100] + F101*fy[N101] + F110*fy[N110] + F111*fy[N111] ;
  	float BBZ = F000*fz[N000] + F001*fz[N001] + F010*fz[N010] + F011*fz[N011] + F100*fz[N100] + F101*fz[N101] + F110*fz[N110] + F111*fz[N111] ;

	if(f == 0)	// horizontal deflectin dipole, field mainly in y direction
	{
  		if(X < 0.)     //   flipping x 
  		{
          		//  By constant, dBy/dy constant, Bz constant, dBz/dz constant  => dBx/dx constant => Bx changes
          		BBX = -BBX ;
  		}
  		if(Y < 0.)       //  flipping y
  		{
          		//  By constant, dBy/dy changes, Bz changes, dBz/dz changes  => dBx/dx changes => Bx changes
          		BBX = -BBX ;
          		BBZ = -BBZ ;
  		}
	}

	if(f == 1)	// vertical dipole, field in x direction
	{
  		if(Y < 0.)     //   flipping y 
  		{
          		//  Bx constant, dBx/dx constant, Bz constant, dBz/dz constant  => dBy/dy constant => By changes
          		BBY = -BBY ;
  		}
  		if(X < 0.)       //  flipping x
  		{
          		//  Bx constant, dBx/dx changes, Bz changes, dBz/dz changes  => dBy/dy changes => By changes
          		BBY = -BBY ;
          		BBZ = -BBZ ;
  		}
	}
// 	if(f > 1)	// do nothing: quadrupoles, full field map defined

	BX += BBX;
	BY += BBY;
	BZ += BBZ;
// 	if(f == 0)	// horizontal deflectin dipole, field mainly in y direction
//   		printf("x %12.4e %12.4e %12.4e %4d %4d %4d %12.4e %12.4e %12.4e\n", X, Y, Z, i, j, k, BX, BY, BZ);
//   		printf("x %12.4e %12.4e %12.4e %4d %4d %4d %12.4e %12.4e %12.4e\n", X, Y, Z, i, j, k, BBX, BBY, BBZ);


   }




  return( 1 ) ;
}

static void csol_exit(struct csol_info *info)
{
	free(info->allFields);		// is used 2 times add flag 
  	free(info) ;
}


