
#include "gdf.hxx"


int main(int argc, char ** argv)
{

	if(argc < 4)
	{
		fprintf(stderr, "usage: %s <gdf output file> <gdf input file>  ... <gdf input file>  \n", argv[0]);
		exit(-1);
	}
		

	int ninputs = argc-2;
	GDF * gdfs[niputs];
	for(int n=0; n<ninputs; n++) gdfs[n] = new GDF(argv[n+2]);
	for(int n=1; n<ninputs; n++) gdfs[0]->addGdf(gdfs[n]);
	for(int n=0; n<ninputs; n++) delete gdfs[n];
	gdfs[0]->saveGdf();

}
