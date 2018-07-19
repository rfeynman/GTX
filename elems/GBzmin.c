/* GBzmin.c: Remove particles with (signed) velocity */

#include <stdio.h>
#include <math.h>
#include "elem.h"

// GBzmin("wcs", "I", 5.0e-10, 0.0); # kill particles going backwards (gamma*beta < 0.0) after t=5e-10 s
struct GBzmin_info
{
  double tmin ;
  double GBzmin ;
} ;
static void GBzmin_sim(double tstart, double tend, double *dt, double *xstart, double *xend, void *vinfo, void *stepinfo);

void GBzmin_init(gptinit *init)
{
  struct GBzmin_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(ECS,Tmin,GBzmin)\n", gptgetname(init) ) ;

  info = (struct GBzmin_info *)gptmalloc( sizeof(struct GBzmin_info) ) ;

  info->tmin     = gptgetargdouble(init,1) ;
  info->GBzmin    = gptgetargdouble(init,2) ;

  /* Register the routine removing the appropriate particles after each succesful step */
  odeaddendfunction( ODEFNC_USR-1, (odeendfnc)GBzmin_sim, info ) ;
}


static void GBzmin_sim(double tstart, double tend, double *dt, double *xstart, double *xend, void *vinfo, void *stepinfo)
{
   struct GBzmin_info *info = (struct GBzmin_info *)vinfo ;

   double GBzmin, tmin;

   if(numpar<2) return ;

  
   tmin = info->tmin ;
   GBzmin = info->GBzmin ;

   if (tend < tmin)
     return;

   for(int i=0 ; i<numpar ; i++) {
     if( pars[i].alive ) {
       if(pars[i].GBr[2] < GBzmin) 
	 gptremoveparticle(&pars[i]) ;
     }
   }
}

