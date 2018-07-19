/* cathode.c: Remove particles with (signed) velocity */

#include <stdio.h>
#include <math.h>
#include "elem.h"

// cathode("wcs", "I", 5.0e-10, 0.0); # kill particles going backwards (gamma*beta < 0.0) after t=5e-10 s
struct cathode_info
{
  double E0 ;
} ;
static void cathode_sim(double tstart, double tend, double *dt, double *xstart, double *xend, void *vinfo, void *stepinfo);

void cathode_init(gptinit *init)
{
  struct cathode_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(ECS,E0)\n", gptgetname(init) ) ;

  info = (struct cathode_info *)gptmalloc( sizeof(struct cathode_info) ) ;

  info->E0     = gptgetargdouble(init,1) ;

  /* Register the routine removing the appropriate particles after each succesful step */
  odeaddendfunction( ODEFNC_USR-1, (odeendfnc)cathode_sim, info ) ;
}


static void cathode_sim(double tstart, double tend, double *dt, double *xstart, double *xend, void *vinfo, void *stepinfo)
{
   struct cathode_info *info = (struct cathode_info *)vinfo ;


   if(numpar<2) return ;

  
   double beta0 = sqrt(2*info->E0/0.511e6) ;


   for(int i=0 ; i<numpar ; i++) {
     if( pars[i].alive ) {
       if(pars[i].Wr[2] < -1e-7) 
       {
	       pars[i].GBr[0] = 0.00;	
	       pars[i].GBr[1] = 0.00;	
	       pars[i].GBr[2] = beta0;	// 0.4 eV
       } 
     }
   }
}

