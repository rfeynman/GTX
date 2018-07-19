
/* onAxis(ECS, "beam", n) sets the n particles with the smallest r to r=0 to used them
   to track a dispersion function
   ******   run with -j 1 option*******
   J Kewisch 2016
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "elem.h"

//#define OFFSET_DEBUG /* Print numpar, avg, offset */

struct onAxis_info
{
	const char * name;
  gptinitpar *pars ;
  int done;
  int * IDs;
  double * R2s;
  int len;
  int numParticle;
} ;

static int onAxis_sim(gptpar *par, double t, struct onAxis_info *info);
static void onAxis_exit(struct onAxis_info *info) ;

void onAxis_init(gptinit *init)
{
  struct onAxis_info *info = (struct onAxis_info *)gptmalloc( sizeof(struct onAxis_info) ) ;
  gptbuildECS( init ) ;

  printf("x\n");

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(set,numParticle)\n", gptgetname(init) ) ;

  const char *name  = gptgetargstring(init,1) ;
  info->numParticle = gptgetargdouble(init,2) ;

  printf("name %s\n", name);
  info->name = strdup(name);

//   /* Get particle set */
//   if( gpttestparset( name )==NULL )
//     gpterror( "The particle set \"%s\" does not exist\n", name ) ;
// 
//   gptparset *set = gptgetparset( name ) ;
//   printf("set %p %s %p\n", set, set->name, set->pars);
// 
//   info->pars = set->pars;
//   info->len = set->len;

//   printf("pars %p\n", info->pars);
//   printf(" x y %e %e\n", info->pars[0].Wr[0], info->pars[0].Wr[1]);
  info->done = 0;
  printf("x\n");
  info->IDs = (int *) gptmalloc(info->numParticle*sizeof(int));
  printf("x\n");
  info->R2s = (double *) gptmalloc(info->numParticle*sizeof(double));
  printf("x\n");
  for(int i=0; i< info->numParticle; i++)  
  {
	  info->IDs[i] =-1;
	  info->R2s[i]=1e33;
  }
  printf("x\n");
  gptaddEBelement( init, onAxis_sim, onAxis_exit, GPTELEM_GLOBAL, info ) ;
  printf("init\n");
}



static int onAxis_sim(gptpar *par, double t, struct onAxis_info *info)
{

	if(info->done) return 0;
  printf("sim %p %p  %d  %s\n", par, info, info->done, info->name);
  gptparset *set = par->set;  //still not right
  printf("set %p %s %p\n", set, set->name, set->pars);
  printf(" x y %e %e\n", info->pars[0].Wr[0], info->pars[0].Wr[1]);

  int len = set->len;
// 	if(Z < 0) return 0;
  /* Scale distribution */
  double GBzmax = -1e33;
  double GBzmin =  1e33;
  for(int i=0 ; i < len ; i++ )
  {
	  gptinitpar * p = (set->pars)+i;
	  double gb = p->GBr[2];
	  if(gb > GBzmax) GBzmax =gb;
	  if(gb < GBzmin) GBzmin =gb;
  }


  for(int n= 0; n< info->numParticle; n++)
  {	
	double gb1 = GBzmin + n* (GBzmax-GBzmin)/(info->numParticle-1);
	double diff =1e33;
	int nb=0;
  	for(int i=0 ; i < len ; i++ )
  	{
	  gptinitpar * p = (set->pars)+i;
	  double gb = p->GBr[2];
	  if( fabs(gb-gb1) < diff)
	  {
		  diff = abs(gb-gb1);
		  nb=i;
	  }
	}
	info->IDs[n]=nb;
  }

  printf("sum\n");
  FILE * fp = fopen("onAxis.dat","w");
  if(!fp)
  {
	  perror("onAxis.dat");
	  exit(-1);
  }
  for(int n= 0; n< info->numParticle; n++)
  {
	  int i=info->IDs[n];
	  set->pars[i].Wr[0]=0.;
	  set->pars[i].Wr[1]=0.;
	  fprintf(fp,"%d\n", set->pars[i].ID);
  }
  fclose(fp);
  info->done=1;


  return 0;

}
static void onAxis_exit(struct onAxis_info *info) 
{
}
