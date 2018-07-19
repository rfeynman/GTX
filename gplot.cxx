
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

	int doTime=0;
	if(argv[1] && !strcmp("-t", argv[1]))
	{
		doTime =1; argv++; argc--;
	}
	if(argc < 2)
	{
		fprintf(stderr, "usage: %s [-t] <gpt output file> [<gpt output file> ....] \n", argv[0]);
		exit(-1);
	}
		
	ScreenList * sc  = new ScreenList;
	GDF * gdf[argc];
	int npp=0;
	for( int a = 1; a < argc; a++)
	{
		gdf[a] = new GDF(argv[a]);
		int np = gdf[a]->getDouble("Npart");
		if(!npp) npp=np;
		else if( np != npp)
		{
			printf("wrong NPART\n");
		}
		for(ScreenList::iterator j=gdf[a]->screenList->begin(); j != gdf[a]->screenList->end(); j++)
		{
			sc->push_back(*j);
		}
		std::sort(sc->begin(), sc->end(), Screen::screenSort);
	}

	char cwd[200];
	getcwd(cwd, 200);
	strcat(cwd,"/");
	strcat(cwd,argv[1]);

	const char * pathlength = doTime ? "Time [s]" : "Path length [m]";
	GRACE env("env.agr", "Envelopes", pathlength, "RMS size [m]");
	env.subTitle(cwd);
	int envH = env.addCurve("horizontal");
	int envV = env.addCurve("vertical");

	GRACE orb("orbit.agr", "Orbits", pathlength, "Offset [m]");
	orb.subTitle(cwd);
	int orbH = orb.addCurve("horizontal");
	int orbV = orb.addCurve("vertical");

	GRACE emitn("emitn.agr", "Normalized Emittance", pathlength, "Normalized Emittance [m]");
	emitn.subTitle(cwd);
	int emitnH = emitn.addCurve("horizontal");
	int emitnV = emitn.addCurve("vertical");

	GRACE emit4("emit4d.agr", "Normalized 4D Emittance", pathlength, "Normalized Emittance [m]");
	emit4.subTitle(cwd);
	int emit4N = emit4.addCurve();

	GRACE cool("cool.agr", "Energy Error", pathlength, "Delta E / E");
	cool.subTitle(cwd);
	int coolN = cool.addCurve();

	GRACE sigp("sigp.agr", "Energy Spread", pathlength, "Delta E / E");
	sigp.subTitle(cwd);
	int sigpN = sigp.addCurve();

	GRACE sigl("sigl.agr", "Bunch length", pathlength, "RMS size [m]");
	sigl.subTitle(cwd);
	int siglN = sigl.addCurve();

	GRACE ener("ener.agr", "Avarage Energy", pathlength, "RMS size [m]");
	ener.subTitle(cwd);
	int enerN = ener.addCurve();

	GRACE time("time_vs_path.agr", "Avarage Time", pathlength, "Time [s]");
	time.subTitle(cwd);
	int timeN = time.addCurve();

	GRACE path("path.agr", "Avarage Pathl ength", "Time [s]", pathlength);
	path.subTitle(cwd);
	int pathN = path.addCurve();




	GRACE csigp("csigp.agr", "Energy Spread", pathlength, "Energy [eV}");
	csigp.subTitle(cwd);
	int csigpN = csigp.addCurve();

	GRACE cemitn("cemitn.agr", "Normalized Emittance", pathlength, "Normalized Emittance [m]");
	cemitn.subTitle(cwd);
	int cemitnH = cemitn.addCurve("horizontal");
	int cemitnV = cemitn.addCurve("vertical");


// 	for(ScreenList::iterator j=gdf->toutList->begin(); j != gdf->toutList->end(); j++)
// 	{
// 		Particles p(*j);
// 		if(p.npart) 
// 		{
// 			p.sortz('a');
// 			p.twiss(0);
// 		}
// 	}

	for(ScreenList::iterator j=sc->begin(); j != sc->end(); j++)
	{
		Particles p(*j);
		if(p.npart) 
		{
			int npart = p.npart;
			p.sortz('a');
			p.twiss(0);
			path.addPoint(pathN, p.avgTime, p.avgPosition);
			double x = doTime ? p.avgTime : p.avgPosition;
			if( npart != npp) 
				printf(" wrong num particles at pos %f\n", x);
			env.addPoint(envH,     x, p.envx);
			env.addPoint(envV,     x, p.envy);
			orb.addPoint(orbH,     x, p.orbx);
			orb.addPoint(orbV,     x, p.orby);
			emitn.addPoint(emitnH, x, p.epsx_n);
			emitn.addPoint(emitnV, x, p.epsy_n);
			emit4.addPoint(emit4N, x, p.emit4d);
			sigp.addPoint(sigpN,   x, p.coolsigma);
			cool.addPoint(coolN,   x, p.sigma_p);
			sigl.addPoint(siglN,   x, p.sigma_l);
			ener.addPoint(enerN,   x, p.ener);
			time.addPoint(timeN,   x, p.avgTime);

			int cut = npart*0.15;
			p.twiss(0, p.dcord+cut, p.npart-2*cut);
			csigp.addPoint(csigpN,   x, p.sigma_p);
			cemitn.addPoint(cemitnH, x, p.epsx_n);
			cemitn.addPoint(cemitnV, x, p.epsy_n);
		}
	}

// 	printf("sc        parm            pos             time       type  npart\n");
// 	for(ScreenList::iterator j=gdf->allList->begin(); j != gdf->allList->end(); j++)
// 	{
// 		Screen * s = *j;
// 		printf("sc  %15.5e %15.5e %15.5e  %4d %4d\n", s->parameter, s->avgPosition, s->avgTime, s->type, s->npart);
// 
// 	}

	delete sc;
	for( int a = 1; a < argc; a++)
	{
		delete gdf[a];
	}










	FILE * fp = fopen("line.agr", "r");
	if(fp)
	{
		Line line;
		double * d = new double[2];
		double maxd = 0.;

		while( fscanf(fp, "%le%le", d, d+1)  == 2)
		{
			if(d[1] > maxd) maxd = d[1];
// 			printf("read %le %le %e\n", d[0], d[1], maxd);
			line.push_back(d);
			d = new double[2];
		}
		delete [] d;
// 		printf("maxd %e\n", maxd);


		for(Line::iterator j=line.begin(); j != line.end(); j++)
		{
			((*j)[1]) /= maxd;
		}

		printf("env\n");
		addLine(env, line);
		printf("orb\n");
		addLine(orb, line);
		printf("emitn\n");
		addLine(emitn, line);
		printf("emit4\n");
		addLine(emit4, line);
		printf("sigp\n");
		addLine(sigp, line);
		printf("cool\n");
		addLine(cool, line);
		printf("sigl\n");
		addLine(sigl, line);
		printf("ener\n");
		addLine(ener, line);
		printf("csigp\n");
		addLine(csigp, line);
		printf("csigl\n");
		addLine(cemitn, line);
		printf("time\n");
		addLine(time, line);
	}
	else
		fprintf(stderr, " line.agr not found\n");
}
