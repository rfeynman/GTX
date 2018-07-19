/*ionizer.c - Kill particle when R is too large */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "elem.h"


const double binwidth =0.001;
struct ionizer_info
{

  /* Coordinate and field information */
  double *z, *r2 ;
  int np ;
  double zmin;
  double zmax;
  double dz;
} ;

static int  ionizer_sim(gptpar *par,double t,struct ionizer_info *info) ;

static  FILE * fp = 0;
static  FILE * fq = 0;

inline double square(double x) { return x*x;}
void ionizer_init(gptinit *init)
{
  struct ionizer_info *info ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=3 )
    gpterror( "Syntax: %s(ECS, filename, z, r)\n", gptgetname(init) ) ;

  info = (struct ionizer_info *)gptmalloc( sizeof(struct ionizer_info) ) ;

  /* GDF information */
  struct gdfmem gm ;
  struct gdff ingdff ;
  const char * filename = gptgetargstring(init,1) ;
  gdfsrinit( filename, &ingdff ) ;
  gdfrmem( &ingdff, &gm, GDFR_ABORTONERROR | GDFR_READARRAYS ) ;
  int points;

  /* Retrieve z array */
  const char * zname = gptgetargstring(init,2) ;
  struct gdfdata * ds=gptinputgetgroup(&gm,filename,zname) ;
  double * zz=gptinputdoublearray(ds,zname,&points) ;

  /* Retrieve radius from datafile */
  const char * fzname   = gptgetargstring(init,3) ;
  double * rr = gptinputdoublearraypoints(ds,fzname,points) ;

  FILE * fa = fopen("appert.agr", "w");
  if(!fa) gpterror( " cant open appert.txt\n");
  for(int i=0; i< points; i++)
	  fprintf(fa, "%e %e\n", zz[i], rr[i]);
  fprintf(fa, "&\n");

  double zmin = zz[0];
  double zmax = zz[points-1];
  int np = int( (zmax-zmin)/binwidth);
  info->zmin = zmin;
  info->zmax = zmax;
  info->np = np;
  double dz = info->dz = (zmax-zmin)/np;

  double * r2 = info->r2 = (double *) gptmalloc( np*sizeof(double) ) ;

  double z0 = zz[0];
  double r0 = rr[0];
  r2[0] = square(r0);
  int j=1;
  for(int i=1; i< points; i++)
  {
	  double r1 = rr[i];
	  double z1 = zz[i];
	  double zt;
	  while( (zt = zmin + j*dz) <= z1)
	  {
	  	r2[j++] = square(   r0 + (r1-r0)*(zt-z0)/(z1-z0)   );
	  }
	  r0 = r1;
	  z0 = z1;
  }
  printf(" points %d  np %d j %d\n", points, np, j);

  for(int i=0; i< np; i++)
	  fprintf(fa, "%e %e\n", zmin + i*dz, sqrt(r2[i]));
  fclose(fa);

  if(!fp)
  {
	fp = fopen("lostAtM.txt", "w");
  	if(!fp) gpterror( " cant open lostAt.txt\n");
  }

  if(!fq)
  {
	fq = fopen("ionizedAtM.txt", "w");
  	if(!fq) gpterror( " cant open lostAt.txt\n");
  }
  gptaddEBelement( init, ionizer_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}

static int ionizer_sim(gptpar *par,double t,struct ionizer_info *info)
{
  if( par->alive )
  {
	int iz = int((Z - info->zmin)/info->dz+0.5);
	double R2 = square(X)+square(Y);
	if(iz < 0 ||  iz >= info->np || R2 >= info->r2[iz]   )    // electron is outside the legal area or pipe radius
	{
  		double E = 0.51104 * (  par->G -1);
		double r = sqrt( R2 );
  		fprintf(fp, "%p  %e %e %e   %e %e %e   %e %e  %e  %d\n", par,
	  		par->Wr[0],
	  		par->Wr[1],
	  		par->Wr[2],
	  		par->GBr[0],
	  		par->GBr[1],
	  		par->GBr[2], t, E, r, par->ID);
  		gptremoveparticle(par) ;
  		return( 1 ) ;
	}

// # ionization cross section of electrons on hydrogen according to
// #  martin reiser: Theory and design of charged particle beams
// 	double betaGamma2 = par->GBr[0]*par->GBr[0] + par->GBr[1]*par->GBr[1] +par->GBr[2]*par->GBr[2];
// 	if(6.5e-05 > betaGamma2) return 0; // below threshhold
// 	double beta2 = betaGamma2/(betaGamma2+1.);
// 	double ff = 6.027e-5/beta2*( 1.659e4*beta2 - 1.);
// 	double sigma = 1.301e-24/beta2*f(beta2)*(log(1.177e5*betaGamma2) - beta2);
// 	double ra = 1e18*sigma*double(random())/double(RAND_MAX);
	double ra = double(random())/double(RAND_MAX);
	if(ra > 0.99999)
	{
		double r = sqrt( R2 );
  		double E = 0.51104 * (  par->G -1);
		printf("ion %e : %e %e %e  %d\n", ra, par->Wr[0],  par->Wr[1],  par->Wr[2], par->ID);
  		fprintf(fq, "%p  %e %e %e   %e %e %e   %e %e  %e  %d\n", par,
	  		par->Wr[0],
	  		par->Wr[1],
	  		par->Wr[2],
	  		par->GBr[0],
	  		par->GBr[1],
	  		par->GBr[2], t, E, r, par->ID);
		return 0;
	}



  }
  return 0;
}
