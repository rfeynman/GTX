
#include "grace_class.hxx"
#include "gdf.hxx"


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
		for(ScreenList::iterator j=gdf[a]->toutList->begin(); j != gdf[a]->toutList->end(); j++)
		{
			sc->push_back(*j);
		}
		fprintf(stderr, "%s %ld %ld\n", argv[a], sc->size(), gdf[a]->toutList->size());
	}
	std::sort(sc->begin(), sc->end(), Screen::screenSort);

	char cwd[200];
	getcwd(cwd, 200);
	strcat(cwd,"/");
	strcat(cwd,argv[1]);
	GRACE tracks("tracks.agr", "Tracks", "[m]", "[m]");
	tracks.subTitle(cwd);

	int id=12, x=0, z=2;
	ScreenList::iterator j=sc->begin(); 
	Screen * s = *j;
	int np= s->npart;
	for(int n=0; n<np; n++)
	{
		int curv = tracks.addCurve("");
	}
	for(ScreenList::iterator k = j+1; k != sc->end(); k++)
	{
		Screen * s = *k;
		for(int n=0; n < s->npart; n++)
		{
			int IDD = round(s->dcord[id][n])-1;
			if(IDD >= 0 && IDD < np)
			{
				tracks.addPoint(IDD,  s->dcord[z][n], s->dcord[x][n]);
// 				printf("%5d %10.4f %10.4f\n", IDD,  s->dcord[z][n], s->dcord[x][n]);
			}
		}
	}


	delete sc;
	for( int a = 1; a < argc; a++)
	{
		delete gdf[a];
	}










}
