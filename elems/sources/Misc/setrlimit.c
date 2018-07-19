/* setrlimit.c: Sets resource limits */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#ifndef _MSC_VER
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h>
#include <errno.h>
#endif
#include "elem.h"

void setrlimit_init(gptinit *init)
{
#ifdef _MSC_VER
  gptwarning("%s: Not supported under MS-windows. Silently ignored.\n", gptgetname(init)) ;
#else
  /* Print usage line when the number of parameters is incorrect */
  int argnum = gptgetargnum(init) ;
  if( argnum!=2 && argnum!=3 )
    gpterror( "Syntax: %s(resource,softlimit,[hardlimit])\n", gptgetname(init) ) ;

  /* Read resource and limit */
  const char *description = gptgetargstring(init,1) ;

  int resource = 0 ;
       if( !strcmp(description,"CORE" ) ) resource=RLIMIT_CORE ;
  else if( !strcmp(description,"CPU"  ) ) resource=RLIMIT_CPU ;
  else if( !strcmp(description,"DATA" ) ) resource=RLIMIT_DATA ;
  else if( !strcmp(description,"FSIZE") ) resource=RLIMIT_FSIZE ;
  else gpterror( "%s: Unknown resource: %s\n", gptgetname(init),description) ;


  /* Try to set new value: get current limit and overwrite with new values */
  struct rlimit limit ;
  if(  getrlimit(resource,&limit) )
    gpterror( "%s: %s\n", gptgetname(init), strerror(errno) ) ;

  limit.rlim_cur = gptgetargdouble(init,2) ;

  if( argnum==3 )
    limit.rlim_max = gptgetargdouble(init,3) ;

  if( setrlimit(resource,&limit) )
    gpterror( "%s: %s\n", gptgetname(init), strerror(errno) ) ;
#endif
}
