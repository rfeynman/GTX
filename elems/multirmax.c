/*multirmax.c - Kill particle when R is too large */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"


const double binwidth =0.001;
struct multirmax_info
{

  /* Coordinate and field information */
  double *z, *r2 ;
  int np ;
  double zmin;
  double zmax;
  double dz;
} ;

static int  multirmax_sim(gptpar *par,double t,struct multirmax_info *info) ;

static  FILE * fp = 0;

inline double square(double x) { return x*x;}
void multirmax_init(gptinit *init)
{
  struct multirmax_info *info ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=3 )
    gpterror( "Syntax: %s(ECS, filename, z, r)\n", gptgetname(init) ) ;

  info = (struct multirmax_info *)gptmalloc( sizeof(struct multirmax_info) ) ;

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
  gptaddEBelement( init, multirmax_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}

static int multirmax_sim(gptpar *par,double t,struct multirmax_info *info)
{
  if( par->alive )
  {

	double R2 = square(par->Wr[0])+square(par->Wr[1]);
  	double E = 0.51104 * (  par->G -1);
	double r = sqrt( R2 );
// 	if(par->Wr[2] < 0. && E > 50. && E < 1400.)	//  according to multipac
	if(par->Wr[2] < 0. )
	{
  		fprintf(fp, "%d  %e %e %e   %e %e %e   %e %e    %e  %e  U\n", par->ID,
	  		par->Wr[0],
	  		par->Wr[1],
	  		par->Wr[2],
	  		par->GBr[0],
	  		par->GBr[1],
	  		par->GBr[2], t, tstart, E, r);
	  	par->Wr[2] = 0.;
	  	par->GBr[0] = 0.;
	  	par->GBr[1] = 0.;
	  	par->GBr[2] = 0.00125;
		par->offWr =1;
		par->offGBr =1;
		return 1;
	}
	int iz = int((par->Wr[2] - info->zmin)/info->dz+0.5);
	if(iz < 0 ||  iz >= info->np || R2 >= info->r2[iz]   )    // electron is outside the legal area or pipe radius
	{
  		fprintf(fp, "%d  %e %e %e   %e %e %e   %e %e    %e  %e \n", par->ID,
	  		par->Wr[0],
	  		par->Wr[1],
	  		par->Wr[2],
	  		par->GBr[0],
	  		par->GBr[1],
	  		par->GBr[2], t, tstart, E, r);
  		gptremoveparticle(par) ;
  		return( 1 ) ;
	}
  }
  return 0;
}
