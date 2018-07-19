#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <error.h>
#include <vector>



main(int argc, char ** argv)
{
	// the file contains two columns 
	// 0 = time from zero to 1/freq
	// 1 = the  positions of the  center of the bunch at that time


	std::vector<double> time, pos;

	FILE * fp = fopen("au10.time.asc", "r");
	if(!fp)
	{
		perror(argv[1]);
		exit(-1);
	}
	
	double t,p;
	while( fscanf(fp, "%le%le", &t, &p) == 2) 
	{
		time.push_back(t);
		pos.push_back(p);
	}
	fclose(fp);
	double freq=703.75e6;
	double t0=1./freq;
	printf("t0 = %e\n", t0);

	int n = time.size();
	double tmax = time[n-1];

	int numperiods = int(tmax/t0);   
	printf("num = %d\n", numperiods);

	double tstep = time[1]-time[0];

	printf("t0/tstep = %f\n", t0/tstep);
	printf("tmax = %e\n", tmax);

	double T = 4e-9;
	double Z = 3.;
	T= atof(argv[1]);
	Z= atof(argv[2]);

	t = T - int(T / t0) * t0;
	printf("t = %e\n", t);

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
