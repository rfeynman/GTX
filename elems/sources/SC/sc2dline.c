/* sc2dline.c - Two Dimensional Space-Charge routine */

/* Author: Kees van der Geer.
 * This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 * Last modification: 12.04.99
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#define ffac 8.98755e9       /* 1 over 4 pi eps_0 */

struct sc_info
{
  double zmin ;
  double zmax ;
} ;

static void sc_sim( double t, double *x, double *p, void *info) ;


void spacecharge2Dline_init(gptinit *init)
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

    struct par *pari = &pars[i] ;
    struct par *parj = &pars[j] ;

    double Rmr[3] ;
    Rmr[0] = pari->Wr[0]-parj->Wr[0] ;
    Rmr[1] = pari->Wr[1]-parj->Wr[1] ;
    Rmr[2] = pari->Wr[2]-parj->Wr[2] ;

    double vj[3] ;
    vj[0] = parj->GBr[0]*gpt_c/parj->G ;
    vj[1] = parj->GBr[1]*gpt_c/parj->G ;
    vj[2] = parj->GBr[2]*gpt_c/parj->G ;

    double vj2 = gptVECSQR( vj ) ;
    if(vj2<=0.0)
        gpterror("Spacecharge2Dline: Particle has zero velocity\n");
      
    double vjabs = sqrt(vj2) ;
    double alphaov = gptVECINP(Rmr,vj)/vj2 ;

    double d[3] ;
    d[0] = -Rmr[0] + alphaov*vj[0] ;
    d[1] = -Rmr[1] + alphaov*vj[1] ;
    d[2] = -Rmr[2] + alphaov*vj[2] ;
    double d2 = gptVECSQR(d) ;
    if (d2 < parj->r2 ) d2 = parj->r2  ;

    double facE = 2*ffac*parj->q*parj->n/(vjabs*d2) ;
    double facB = 2e-7*parj->q*parj->n/(vjabs*d2) ;

    pari->WE[0] -= facE*d[0] ;
    pari->WE[1] -= facE*d[1] ;
    pari->WE[2] -= facE*d[2] ;

    pari->WB[0] += facB*(d[1]*vj[2] - d[2]*vj[1]) ;
    pari->WB[1] += facB*(d[2]*vj[0] - d[0]*vj[2]) ;
    pari->WB[2] += facB*(d[0]*vj[1] - d[1]*vj[0]) ;
  }
}
