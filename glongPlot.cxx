
#include "grace_class.hxx"
#include "gdf.hxx"

typedef std::vector<double *> Line;



void addLine(GRACE & g, Line & line)
{
	double xmax = g.maxX();
	double ymax = g.maxY();
// 	printf("max %e %e\n", xmax, ymax);
	int c = g.addCurve();
	for(Line::iterator j=line.begin(); j != line.end(); j++)
	{
		double * d = *j;
// 		printf("l %e %e\n", d[0], d[1]);
		if(  d[0] <= xmax)
			g.addPoint(c, d[0], d[1]*ymax*0.3);
	}
}

int main(int argc, char ** argv)
{

	if(argc > 4)
	{
		fprintf(stderr, "usage: %s <gpt output file>  [t|p]  <pos> [<pos>...]\n", argv[0]);
		exit(-1);
	}
		
	GDF * gdf = new GDF(argv[1]);

	char  what = argv[2][0];

	char cwd[200];
	getcwd(cwd, 200);
	strcat(cwd,"/");
	strcat(cwd,argv[1]);

	for(int i=3; i<argc; i++)
	{

		double pos=atof(argv[i]);
		char filename[100];
		strcpy(filename, "pos");
		sprintf(filename+3, "%06.3f", pos);
		printf("file = %s\n", filename);
		GRACE phase(filename, "Longitudinal Phase", "Time [s]", "Energy [eV]");
		phase.subTitle(cwd);
		int phaseN = phase.addCurve();
		phase.mkdots(phaseN);

		for(ScreenList::iterator j=gdf->screenList->begin(); j != gdf->screenList->end(); j++)
		{
			Particles  p(*j);
			if(p.npart) 
			{
				if( fabs(p.avgPosition - pos) < 1e-2)
				{
					for(int n=0; n< p.npart; n++)
						phase.addPoint(phaseN, p.dcord[n][4], p.dcord[n][5]);
				}
			}
		}
	}




}
