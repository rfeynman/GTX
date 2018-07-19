/* sc3dclas.c - Classical 3D Point-to-point space-Charge routine. */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#define ffac (1/(4*gpt_pi*gpt_eps0))

struct sc_info
{
  double zmin ;
  double zmax ;
} ;

static void sc_sim( double t, double *x, double *p, void *info) ;


void spacecharge3Dclassic_init(gptinit *init)
{
  struct sc_info *info ;
  int numarg ;

  numarg = gptgetargnum(init) ;
  if( numarg!=0 && numarg!=2 )
    gpterror( "Syntax: %s([zmin,zmax])", gptgetname(init) ) ;

  info = (struct sc_info *)gptmalloc( sizeof(struct sc_info) ) ;

  if( numarg==2 )
  {
    info->zmin=gptgetargdouble(init,1) ;
    info->zmax=gptgetargdouble(init,2) ;
  } else
  {
    info->zmin=-DBL_MAX ;
    info->zmax= DBL_MAX ;
  }

  odeaddfprfunction( ODEFNC_USR, sc_sim,info ) ;
}

static void sc_sim( double t, double *x, double *p, void *vinfo)
{
  const struct sc_info *info = (struct sc_info *)vinfo ;
  const double zmin = info->zmin ;
  const double zmax = info->zmax ;

#pragma omp parallel for schedule(static) num_threads(gptCPUs) if(gptCPUs!=1 && numpar > 16*gptCPUs)
  for(int i=0 ; i<numpar ; i++)
    if( pars[i].alive && pars[i].Wr[2] > zmin && pars[i].Wr[2] < zmax )
      for(int j=0 ; j<numpar ; j++)
  {
    if( i==j || !pars[j].alive ) continue ;
    if( pars[j].Wr[2] < zmin || pars[j].Wr[2] > zmax ) continue ;

    double rji[3] ;
    rji[0] = pars[i].Wr[0] - pars[j].Wr[0] ;
    rji[1] = pars[i].Wr[1] - pars[j].Wr[1] ;
    rji[2] = pars[i].Wr[2] - pars[j].Wr[2] ;

    double r2 = rji[0]*rji[0] + rji[1]*rji[1] + rji[2]*rji[2]  + pars[j].r2 ;
	double r3 = r2 * sqrt(r2) ;
	double mfac = ffac*pars[j].q*pars[j].n/r3 ;

    pars[i].WE[0] += mfac*rji[0] ;
    pars[i].WE[1] += mfac*rji[1] ;
    pars[i].WE[2] += mfac*rji[2] ;
  }
}
