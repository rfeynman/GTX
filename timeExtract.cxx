
#include "gdf.hxx"


int main(int argc, char ** argv)
{

	if(argc != 4)
	{
		fprintf(stderr, "usage: %s <gdf input file> <gdf output file>  <time>\n", argv[0]);
		exit(-1);
	}
		
	double t=atof(argv[3]);
	GDF * gdf = new GDF(argv[1]);
	gdf->writeTimeStep(argv[2], t);


}
