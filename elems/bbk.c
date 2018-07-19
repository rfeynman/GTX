/* bbk3Sym_B.c - 3D rectangular field map for the magnetic field */

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
#include "elem.h"


struct bbk_info
{
  /* GDF information */
  struct gdfmem gm ;

  double * time;
  double * avgx;
  double * avgy;
  double * avgz;
  double * stdx;
  double * stdy;
  double * stdz;
  double * avgE;
  double * stdG;

  double freq, charge, t0, tmax, tstep;
  int n, numperiods;
} ;




void bbk3Sym_B_init(gptinit *init)
{
  struct bbk_info *info ;
  struct gdfdata *ds ;

  /* Commandline parameters */
  int numc ;
  const char *fxname,*fyname,*fzname ;
  const char *xname,*yname,*zname ;
  double ffac ;


  /* GDF information */
  struct gdff ingdff ;
  const char *filename ;

  gptbuildECS( init ) ;

  /* Commandline parsing */
  numc=gptgetargnum(init) ;
  if( numc!=3 )
    gpterror( "Syntax: %s(ECS, bunchCoordinates.gdf, freq, charge)\n", gptgetname(init) ) ;

  info = (struct bbk_info *)gptmalloc( sizeof(struct bbk_info) ) ;

  /* Initialize GDF */
  filename     = gptgetargstring(init,1) ;
  info->freq   = gptgetargdouble(init,2) ;
  info->charge = gptgetargdouble(init,3) ;
  gdfsrinit( filename, &ingdff ) ;
  gdfrmem( &ingdff, &info->gm, GDFR_ABORTONERROR | GDFR_READARRAYS ) ;

  // time avgx avgy avgz  stdx stdy stdz  avgE  stdG
  /* Retrieve x, y and z arrays */

  int points;
  struct gdfdata * ds = gptinputgetgroup(&info->gm,filename,xname) ;
  info->time=gptinputdoublearray(ds,"time",&points) ;
  info->avgx=gptinputdoublearraypoints(ds,"avgx",points) ;
  info->avgy=gptinputdoublearraypoints(ds,"avgy",points) ;
  info->avgz=gptinputdoublearraypoints(ds,"avgz",points) ;
  info->stdx=gptinputdoublearraypoints(ds,"stdx",points) ;
  info->stdy=gptinputdoublearraypoints(ds,"stdy",points) ;
  info->stdz=gptinputdoublearraypoints(ds,"stdz",points) ;
  info->avgE=gptinputdoublearraypoints(ds,"avgE",points) ;
  info->stdG=gptinputdoublearraypoints(ds,"stdG",points) ;


  info->n=points;
  info->tstep = info->time[1]-info->time[0];
  info->numperiods = int(info->tmax/info->t0);   
  printf("num = %d\n", info->numperiods);


  gptaddEBelement( init, bbk_sim, bbk_exit, GPTELEM_GLOBAL, info ) ;
}


static int bbk_sim(gptpar *par, double t, struct bbk_info *info)
{
	t = t - int(t / info->t0) * info->t0;
	double dpp = 1e33;
	double dzz = 1e33;
	for(int i=0; i < numperiods; i++)
	{
		double tt = t + i*t0;
		int itt = int(tt/tstep);
		double pp = pos[itt] + (pos[itt+1]-pos[itt])*(tt-time[itt])/(time[itt+1]-time[itt]);
		double dz = fabs(Z-pp);
		if( dz < dpp )
		{
			dpp = dz;
			dzz = Z - pp;
		}

		printf("%d %d %e %e %e %e\n", i, itt, time[itt], pos[itt], pos[itt+1], pp);
	}
		
	printf("dpp = %e\n", dpp);
	printf("dzz = %e\n", dzz);
}
