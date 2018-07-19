#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <error.h>
#include <vector>

const double clight          =    2.99792456210e8;

void mkdots(FILE * fp, int i)
{
	fprintf(fp, "@with g0\n");
	fprintf(fp, "@    s%d symbol 9\n",i);
	fprintf(fp, "@    s%d symbol size 0.060000\n",i);
	fprintf(fp, "@    s%d line type 0\n",i);
}

// this reads the gdf ascii file. readgdf reads the binary files
// readGDF <input file> [ <x|y|r|l> [time|position|ID] <value>]

int main(int argc, char ** argv)
{
	FILE * fp = fopen(argv[1], "r");
	if(!fp)
	{
		perror(argv[1]);
		exit(-1);
	}
	FILE * fo;

	char line[1000];
	char word[1000];
	int sofar;
	char pp[100];
	char what[100];
	strcpy(pp, "");

	int plot = 0;
	 if( argc == 5)
	 {
		 if( !strcmp(argv[3], "time") ) plot = 1;
		 if( !strcmp(argv[3], "position") ) plot = 2;
		 if( !strcmp(argv[3], "ID") ) plot = 4;
	 }
	 printf("%d %d\n" , argc, plot);
	 if(plot)
	 {
		fo = fopen("_gpt.agr", "w");
		if(!fo)
		{
			perror("_gpt.agr");
			exit(-1);
		}
		strcpy(what, argv[2]);
		strcpy(pp, argv[4]);
		mkdots(fo, 0);
	 }





	while(  fgets(line, 1000, fp)  )
	{
		char comnd[100];
		char timex[100];
		char posix[100];
		sscanf(line, "%s", comnd);

		if(   !strcmp("time", comnd)  )
		{
			double 	 x, y, z, G, Bx, By, Bz, rxy, m, q, nmacro, rmacro, ID, fEx, fEy, fEz, fBx, fBy, fBz;
			sscanf(line, "%s%s", comnd, timex);
			printf("%s   %s ",  comnd, timex);
			fgets(line, 1000, fp); // column header
			int nh=0;
			char * p = line;
			while( sscanf(p, "%s%n", word, &sofar) == 1 )
			{
				printf("%s ", word);
				p += sofar;
				nh++;
			}
			printf(" %d\n", nh);

			while(  fgets(line, 1000, fp)  )
			{
				 int num = sscanf(line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
						 &x, &y, &z, &G, &Bx, &By, &Bz, &rxy, &m, &q, &nmacro, &rmacro, &ID, &fEx, &fEy, &fEz, &fBx, &fBy, &fBz);
				 if(num != 19) break;
				 if( plot & 1)
				 {
					 fprintf(fo, "%e %e\n", z, G*Bz);
				 }

			}

			continue;
		}
		else
		if(   !strcmp("position", comnd)  )
		{
			enum cols { 	 x, y, z, rxy, Bx, By, Bz, G,  t, m, q, nmacro, rmacro, ID };
			std::vector<double *> particles;
// 					&x, &y, &z, &rxy, &Bx, &By, &Bz, &G,  &t, &m, &q, &nmacro, &rmacro, &ID);
			sscanf(line, "%s%s", comnd, posix);
			printf("%s   %s  ",  comnd, posix);
			int doplot = ! strcmp(posix, pp);
			fgets(line, 1000, fp); // column header
			int nh=0;
			char * p = line;
			while( sscanf(p, "%s%n", word, &sofar) == 1 )
			{
				printf("%s ", word);
				p += sofar;
				nh++;
			}
			printf(" %d\n", nh);
			if(doplot)
			{
				fprintf(stderr, "plot at %s   %s  mask %d \n",  comnd, posix, ((plot & 2)  && doplot)   );
			}
			while(  fgets(line, 1000, fp)  )
			{
				double * data = new double[14];
				 int num = sscanf(line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
						 data, data+1, data+2, data+3, data+4, data+5, data+6, data+7, data+8, data+9, data+10, data+11, data+12, data+13);
				 if(num != 14) break;

				 if( (plot & 2)  && doplot)
				 {
					 particles.push_back(data);
				 }
			}

			if( (plot & 2)  && doplot)
			{
				fprintf(stderr, "vector length %ld \n", particles.size() );
				int nd=0;
				double zavg=0.;
				double pavg =0.;
				for (int i=0; i<particles.size(); i++)
				{
					 nd++;
					 double z = particles[i][t]*703.75e6*360.;
// 					 double gbz = particles[i][G]*particles[i][Bz];
					 double gbz = particles[i][G] - 1.;    //  E/m_e
					 zavg += z;
					 pavg += gbz;
				}
				zavg /= nd;
				pavg /= nd;
	
				char cmd[200];
	        		fprintf(fo, "@with g0\n");
				sprintf(cmd, "%s, s=%s", comnd, posix);
	        		fprintf(fo, "@subtitle \"%s\"\n", cmd); 
				if( what[0] == 'x')
				{
	        			fprintf(fo, "@xaxis label \"x [m]\"\n");
	        			fprintf(fo, "@yaxis label \"x'\"\n");
				}
				if( what[0] == 'y')
				{
	        			fprintf(fo, "@xaxis label \"y [m]\"\n");
	        			fprintf(fo, "@yaxis label \"y'\"\n");
				}
				if( what[0] == 'r')
				{
	        			fprintf(fo, "@xaxis label \"r [m]\"\n");
	        			fprintf(fo, "@yaxis label \"r'\"\n");
				}
				if( what[0] == 'l')
				{
	        			fprintf(fo, "@xaxis label \"RF Phase\"\n");
	        			fprintf(fo, "@yaxis label \"Delta P/P0\"\n");
				}
	        		fprintf(fo, "@hardcopy device \"PNG\"\n");
	        		fprintf(fo, "@device \"PNG\" DPI 1200\n");
	
	
	
				for (int i=0; i<particles.size(); i++)
				{
					double ax = 0;
					double ay = 0;
					if( what[0] == 'x')
					{
					 	ax = particles[i][x];
					 	ay = particles[i][Bx]/particles[i][Bz];
					}
					if( what[0] == 'y')
					{
					 	ax = particles[i][y];
					 	ay = particles[i][By]/particles[i][By];
					}
					if( what[0] == 'r')
					{
					 	ax = particles[i][rxy];
						ay = sqrt(particles[i][x]*particles[i][Bx]/particles[i][Bz] + particles[i][y]*particles[i][By]/particles[i][Bz]);
					}
					if( what[0] == 'l')
					{
					 	ax = particles[i][t]*703.75e6*360. - zavg;	// delta z
	// 					 ax = -particles[i][t]*particles[i][Bz]*clight - zavg;
	// 					 ay = (particles[i][G]*particles[i][Bz] - pavg)/pavg;
					 	ay = (particles[i][G] - 1. - pavg)/pavg;   //  delta E over E
					}
					fprintf(fo, "%e %e\n", ax,ay);
				}
			}
			continue;
		}
		else
			printf("%s", line);


	}
}
