/* map3D_B.c - 3D rectangular field map for the magnetic field */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <algorithm>
using namespace std ;

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include "elem.h"


struct map3d_info
{
  /* GDF information */
  struct gdfmem gm ;

  /* Coordinate information, only used for sorting */
  double *x, *y, *z ;

  /* Array size information */
  double xmin, xmax, xdelta ;
  double ymin, ymax, ydelta ;
  double zmin, zmax, zdelta ;
  int numx, numy, numz ;

  /* Field information */
  double *fx, *fy, *fz;
} ;


static int map3d_sim(gptpar *par, double t, struct map3d_info *info) ;
static void map3d_exit(struct map3d_info *info) ;


/* Helper functions for sorting */
static void swapmap3D(struct map3d_info *info, int i, int j) ;
static void qsortmap3D(struct map3d_info *info, int left, int right) ;


void map3D_B_init(gptinit *init)
{
  struct map3d_info *info ;
  struct gdfdata *ds ;

  /* Commandline parameters */
  int numc ;
  const char *fxname,*fyname,*fzname ;
  const char *xname,*yname,*zname ;
  double ffac ;

  /* Grid characteristics */
  int i,j,k,points,N ;

  /* GDF information */
  struct gdff ingdff ;
  const char *filename ;

  gptbuildECS( init ) ;

  /* Commandline parsing */
  numc=gptgetargnum(init) ;
  if( numc!=8 )
    gpterror( "Syntax: %s(ECS,mapfile.gdf,x,y,z,Bx,By,Bz,Bfac)\n", gptgetname(init) ) ;

  info = (struct map3d_info *)gptmalloc( sizeof(struct map3d_info) ) ;

  /* Initialize GDF */
  filename = gptgetargstring(init,1) ;
  gdfsrinit( filename, &ingdff ) ;
  gdfrmem( &ingdff, &info->gm, GDFR_ABORTONERROR | GDFR_READARRAYS ) ;

  /* Retrieve x, y and z arrays */
  xname = gptgetargstring(init,2) ;
  yname = gptgetargstring(init,3) ;
  zname = gptgetargstring(init,4) ;

  ds=gptinputgetgroup(&info->gm,filename,xname) ;
  info->x=gptinputdoublearray(ds,xname,&points) ;
  info->y=gptinputdoublearraypoints(ds,yname,points) ;
  info->z=gptinputdoublearraypoints(ds,zname,points) ;

  /* Retrieve field from datafile */
  fxname = gptgetargstring(init,5) ;
  fyname = gptgetargstring(init,6) ;
  fzname = gptgetargstring(init,7) ;

  info->fx = gptinputdoublearraypoints(ds,fxname,points) ;
  info->fy = gptinputdoublearraypoints(ds,fyname,points) ;
  info->fz = gptinputdoublearraypoints(ds,fzname,points) ;

  ffac = gptgetargdouble(init,8) ;
  for(i=0 ; i<points ; i++) info->fx[i] *= ffac ;
  for(i=0 ; i<points ; i++) info->fy[i] *= ffac ;
  for(i=0 ; i<points ; i++) info->fz[i] *= ffac ;


  /* Sort all datapoints in correct order: x runs fastest */
  for(i=0 ; i<points ; i++) swapmap3D(info,rand()%points,rand()%points) ;
  qsortmap3D(info,0,points-1) ;

  /* Obtain file characteristics */
  info->xmin = info->x[0] ; 
  info->ymin = info->y[0] ; 
  info->zmin = info->z[0] ; 
  info->xmax = info->x[points-1] ; 
  info->ymax = info->y[points-1] ; 
  info->zmax = info->z[points-1] ;
  i=1 ;
  while(i<points && info->x[i-1]<info->x[i]) i++ ;
  info->numx = i ;
  while(i<points && info->y[i-info->numx]<info->y[i]) i+=info->numx ;
  info->numy = i/info->numx ;
  info->numz = points/info->numx/info->numy ;
  if( info->numx<2 ) gpterror( "%s: All X-coordinates are equal to %g\n", filename, info->xmin ) ;
  if( info->numy<2 ) gpterror( "%s: All Y-coordinates are equal to %g\n", filename, info->ymin ) ;
  if( info->numz<2 ) gpterror( "%s: All Z-coordinates are equal to %g\n", filename, info->zmin ) ;
  info->xdelta = (info->xmax-info->xmin)/(info->numx-1) ; 
  info->ydelta = (info->ymax-info->ymin)/(info->numy-1) ; 
  info->zdelta = (info->zmax-info->zmin)/(info->numz-1) ;

  /* Test file for correctness */
  for(k=0 ; k<info->numz ; k++) for(j=0 ; j<info->numy ; j++) for(i=0 ; i<info->numx ; i++)
  { 
    N=(k*info->numy+j)*info->numx+i ;
    if( info->x[N] < info->xmin+i*info->xdelta - info->xdelta/8 ||
        info->x[N] > info->xmin+i*info->xdelta + info->xdelta/8 ||
        info->y[N] < info->ymin+j*info->ydelta - info->ydelta/8 ||
        info->y[N] > info->ymin+j*info->ydelta + info->ydelta/8 ||
        info->z[N] < info->zmin+k*info->zdelta - info->zdelta/8 ||
        info->z[N] > info->zmin+k*info->zdelta + info->zdelta/8 )
    {
      gptwarning( "x-direction: %d samples in range [%g,%g] with delta %g\n", info->numx, info->xmin, info->xmax, info->xdelta ) ;
      gptwarning( "y-direction: %d samples in range [%g,%g] with delta %g\n", info->numy, info->ymin, info->ymax, info->ydelta ) ;
      gptwarning( "z-direction: %d samples in range [%g,%g] with delta %g\n", info->numz, info->zmin, info->zmax, info->zdelta ) ;
      gpterror( "%s: No rectangular grid with the above parameters at (%g,%g,%g)\n", filename, info->x[N], info->y[N], info->z[N] ) ;
    }
  }

  gptaddEBelement( init, map3d_sim, map3d_exit, GPTELEM_GLOBAL, info ) ;
}


static int map3d_sim(gptpar *par, double t, struct map3d_info *info)
{
  int i,j,k,N ;
  double ft,fu,fv, gt,gu,gv ;
  int N1,N2,N3,N4,N5,N6,N7,N8 ;
  double f1,f2,f3,f4,f5,f6,f7,f8 ;

  if( Z<info->zmin || Z>info->zmax ) return( 0 ) ;
  if( Y<info->ymin || Y>info->ymax ) return( 0 ) ;
  if( X<info->xmin || X>info->xmax ) return( 0 ) ;

  /* Calculate master index and offset fractions */
  i = (int)((X-info->xmin)/info->xdelta) ;
  j = (int)((Y-info->ymin)/info->ydelta) ;
  k = (int)((Z-info->zmin)/info->zdelta) ;
  N = (k*info->numy+j)*info->numx+i ;
  ft= (X-info->xmin-i*info->xdelta)/info->xdelta ;
  fu= (Y-info->ymin-j*info->ydelta)/info->ydelta ;
  fv= (Z-info->zmin-k*info->zdelta)/info->zdelta ;
  gt = 1-ft ;
  gu = 1-fu ;
  gv = 1-fv ;

  /* Calculate boundary offsets */
  N1 = N ;
  N2 = N+1 ;
  N3 = N+info->numx+1 ;
  N4 = N+info->numx ;
  N5 = N+info->numx*info->numy ;
  N6 = N+info->numx*info->numy+1 ;
  N7 = N+info->numx*(info->numy+1)+1 ;
  N8 = N+info->numx*(info->numy+1) ;

/* Consistency check in case of problems with this element
 * if(i<0 || i>=info->numx-1) gpterror( "Error in x range: %lf maps to %d\n", X, i ) ;
 * if(j<0 || j>=info->numy-1) gpterror( "Error in y range: %lf maps to %d\n", Y, j ) ;
 * if(k<0 || k>=info->numz-1) gpterror( "Error in z range: %lf maps to %d\n", Z, k ) ;
 * if(ft<0 || ft>1+16*DBL_EPSILON ) gpterror( "x fraction is %.16f at x=%f\n", ft, X ) ;
 * if(fu<0 || fu>1+16*DBL_EPSILON ) gpterror( "y fraction is %.16f at y=%f\n", fu, Y ) ;
 * if(fv<0 || fv>1+16*DBL_EPSILON ) gpterror( "z fraction is %.16f at z=%f\n", fv, Z ) ;
 */

  f1=gt*gu*gv ; f2=ft*gu*gv ; f3=ft*fu*gv ; f4=gt*fu*gv ; 
  f5=gt*gu*fv ; f6=ft*gu*fv ; f7=ft*fu*fv ; f8=gt*fu*fv ; 
  BX = f1*info->fx[N1] + f2*info->fx[N2] + f3*info->fx[N3] + f4*info->fx[N4] + 
       f5*info->fx[N5] + f6*info->fx[N6] + f7*info->fx[N7] + f8*info->fx[N8] ;
  BY = f1*info->fy[N1] + f2*info->fy[N2] + f3*info->fy[N3] + f4*info->fy[N4] +
       f5*info->fy[N5] + f6*info->fy[N6] + f7*info->fy[N7] + f8*info->fy[N8] ;
  BZ = f1*info->fz[N1] + f2*info->fz[N2] + f3*info->fz[N3] + f4*info->fz[N4] +
       f5*info->fz[N5] + f6*info->fz[N6] + f7*info->fz[N7] + f8*info->fz[N8] ;

  return( 1 ) ;
}

static void map3d_exit(struct map3d_info *info)
{
  gdfmfreechilds( &info->gm.ds ) ;
  free(info) ;
}


/* swap two points, keeping all arrays in sync */
static void swapmap3D(struct map3d_info *info, int i, int j)
{
  if( i!=j &&
      info->x[i]==info->x[j] &&
      info->y[i]==info->y[j] &&
      info->z[i]==info->z[j] )
    gpterror( "Duplicate points at (%g,%g,%g)\n", info->x[i], info->y[i], info->z[i] ) ;

  swap(info->x[i],info->x[j]) ;
  swap(info->y[i],info->y[j]) ;
  swap(info->z[i],info->z[j]) ;

  swap(info->fx[i],info->fx[j]) ;
  swap(info->fy[i],info->fy[j]) ;
  swap(info->fz[i],info->fz[j]) ;
}


/* Sort all datapoints in correct order: x runs fastest */
static void qsortmap3D(struct map3d_info *info, int left, int right)
{
  int i, last ;

  if( left>=right ) return ;

  swapmap3D(info,left,(left+right)/2) ;
  last = left ;
  for(i=left+1 ; i<=right ; i++)
    if(info->z[i]<info->z[left] ||
      (info->y[i]<info->y[left] && info->z[i]==info->z[left]) ||
      (info->x[i]<info->x[left] && info->z[i]==info->z[left] && info->y[i]==info->y[left]))
        swapmap3D(info,++last,i) ;
  swapmap3D(info,left,last) ;
  qsortmap3D(info,left,last-1) ;
  qsortmap3D(info,last+1,right) ;
}
