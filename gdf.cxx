#include "gdf.hxx"

//                                0    1    2     3     4     5     6     7    8    9   10   11      12      13        14    15						screen
const char * columnp[] = {"position", "x", "y", "z", "rxy", "Bx", "By", "Bz", "G", "t", "m", "q", "nmacro", "rmacro", "ID", "empty"};
//                            0    1    2    3    4    5     6     7      8     9   10       11      12      13     14   						snapshot
const char * columns[] = {"time", "x", "y", "z", "G", "Bx", "By", "Bz", "rxy", "m", "q", "nmacro", "rmacro", "ID", "empty"};
//                            0    1    2    3    4    5     6     7      8     9   10       11      12      13     14    15      16     17     18     10     20	tout
const char * columnt[] = {"time", "x", "y", "z", "G", "Bx", "By", "Bz", "rxy", "m", "q", "nmacro", "rmacro", "ID", "fEx", "fEy", "fEz", "fBx", "fBy", "fBz", "empty"};



GDF::~GDF()
{
	for(GdfItemList::iterator j=gdfItemList->begin(); j != gdfItemList->end(); j++) {delete (*j); *j=0; }
	delete gdfItemList;
	for(ScreenList::iterator j=screenList->begin(); j != screenList->end(); j++)    {delete (*j); *j=0; }
	delete screenList;
	for(ScreenList::iterator j=toutList->begin(); j != toutList->end(); j++)        {delete (*j); *j=0; }
	delete toutList;
	delete allList;
	delete [] buff;
}

inline uchar GDF::readUChar(char  * &p)
{
	uchar i;
	memcpy(&i, p, 1);
	p += 1;
	return i;
}


inline int GDF::readInt(char  * &p)
{
	int i;
	memcpy(&i, p, 4);
	p += 4;
	return i;
}


inline uint GDF::readUInt(char  * &p)
{
	uint i;
	memcpy(&i, p, 4);
	p += 4;
	return i;
}


inline float GDF::readFloat(char  * &p)
{
	float f;
	memcpy(&f, p, 4);
	p += 4;
	return f;
}


inline double  GDF::readDouble(char  * &p)
{
	double  d;
	memcpy(&d, p, 8);
	p += 8;
	return d;
}



GDF::GDF(const char * filename)
{






	int fi = open(filename,O_RDONLY);
        if(fi < 0)
        {
                perror(filename);
                exit(0);
        }
        size = lseek(fi, 0, SEEK_END);
        lseek(fi, 0, SEEK_SET);
//         printf("size is %ld\n", size);
        buff = new char[size];
        ssize_t rsize = read(fi, buff, size);
//         printf("read %ld\n", rsize);
        close(fi);
	if(size != rsize)
	{
		fprintf(stderr, "read %ld bytes, sould be %ld\n", rsize, size);
		exit(-1);
	}
	char * pb = buff;
	char * pend = buff + size;

	gdfItemList = new GdfItemList;


	uint magic = readUInt(pb); 
	if(magic != GDFID)
	{
		fprintf(stderr, "wrong magic number\n");
		exit(-1);
	}

	uint cretime  = readUInt(pb);
	char creator[GDFNAMELEN]; memcpy(creator, pb, GDFNAMELEN); pb += GDFNAMELEN;
	char destin[GDFNAMELEN]; memcpy(destin, pb, GDFNAMELEN); pb += GDFNAMELEN;
	uchar gdfmaj   = readUChar(pb);
	uchar gdfmin   = readUChar(pb);
	uchar cremaj   = readUChar(pb);
	uchar cremin   = readUChar(pb);
	uchar desmaj   = readUChar(pb);
	uchar desmin   = readUChar(pb);
	uchar dummy1   = readUChar(pb);
	uchar dummy2   = readUChar(pb);


// 	printf("creator %s\n", creator);
// 	printf("destin  %s\n", destin);
// 	printf("cretime %u\n", cretime);
// 	printf("gdfmaj  %u\n", gdfmaj);
// 	printf("gdfmin  %u\n", gdfmin);
// 	printf("cremaj  %u\n", cremaj);
// 	printf("cremin  %u\n", cremin);
// 	printf("desmaj  %u\n", desmaj);
// 	printf("desmin  %u\n", desmin);
// 	printf("dummy1  %u\n", dummy1);
// 	printf("dummy2  %u\n", dummy2);

	int re = pend - pb;

	while(pb < pend)
	{
		GdfItem * thisItem = new GdfItem;
		thisItem->startPtr = pb;
		int itemLength = GDFNAMELEN+8;
		memcpy(thisItem->name, pb, GDFNAMELEN); pb += GDFNAMELEN;
		if(!strlen(thisItem->name)) strcpy(thisItem->name, "empty");
		else if (thisItem->name[0] == '@') strcpy(thisItem->name, "x0x40logo"); 
		thisItem->block_type     = readUInt(pb);
		thisItem->block_size     = readUInt(pb);

		uint data_type      = thisItem->block_type & 0xff;
		uint start_dir      = thisItem->block_type & t_dir;
		uint end_dir        = thisItem->block_type & t_edir;	// end folder  
		uint data_is_param  = thisItem->block_type & t_param;	// if 1, the data that is coming next in the file contains the scanned parameters when using MR, time when using tout, position for screen
		uint data_is_array  = thisItem->block_type & t_data;	// if 1, the data that is coming next in the file contains particle info: position, momenta, electromagnetic fields

// 		printf("item %-20s %10x %10d  %x  ", thisItem->name, thisItem->block_type, thisItem->block_size, data_type);
// 		if(start_dir) printf("start_dir  ");
// 		if(end_dir) printf("end_dir  ");
// 		if(data_is_param) printf("data_is_param  ");
// 		if(data_is_array) printf("data_is_array  ");
// 		printf("\n");


	// read a scanned parameters when using MR, time when using tout, position for screen
    		if(data_is_param)     
    		{
        		switch(data_type)
			{
				case t_s32:
				{
                			thisItem->ivalue = readInt(pb);                   //read param
					itemLength += 4;
// 					printf("data %d\n", thisItem->ivalue);
					break;
				}
				case t_flt:
				{
                			thisItem->fvalue = readFloat(pb);                   //read param
					itemLength += 4;
// 					printf("data %e\n", thisItem->dvalue);
					break;
				}
				case t_dbl:
				{
                			thisItem->dvalue = readDouble(pb);                   //read param
					itemLength += 8;
// 					printf("data %e\n", thisItem->dvalue);
					break;
				}
				case t_nul:
				{
                			// no data present
// 					printf("data empty\n");
					break;
				}
				case t_ascii:
				{
					thisItem->charvalue = pb; pb += thisItem->block_size;
					itemLength += thisItem->block_size;
// 					printf("datac %s\n", thisItem->charvalue);
					break;
				}
				default:
				{
       			         	// error, abort while loop
       				         fprintf(stderr, "Error in ''load_gdf'': unknown datatype of value\n");
					 exit(-1);
					break;
				}
			}
    		}


    // read data array with particle info: position, momenta, electromagnetic fields
    		if(data_is_array)
    		{
        		switch(data_type)
			{
				case t_flt:
				{
					thisItem->block_arrsize = thisItem->block_size / sizeof(float);
                			if (thisItem->block_size != thisItem->block_arrsize * sizeof(float) )
					{ 
       				         	fprintf(stderr, "Error in ''load_gdf'': wrong size for double array\n");
		    			}
					thisItem->farray = (float *) pb; pb += thisItem->block_size;
					itemLength += thisItem->block_size;
// 					printf("dataa %d\n", thisItem->block_size);
					break;
				}
				case t_dbl:
				{
					thisItem->block_arrsize = thisItem->block_size / sizeof(double);
                			if (thisItem->block_size != thisItem->block_arrsize * sizeof(double) )
					{ 
       				         	fprintf(stderr, "Error in ''load_gdf'': wrong size for double array\n");
		    			}
					thisItem->darray = (double *) pb; pb += thisItem->block_size;
					itemLength += thisItem->block_size;
// 					printf("dataa %d\n", thisItem->block_size);
					break;
				}

				default:
				{
                			// error, abort while loop
       				         fprintf(stderr, "Error in ''load_gdf'': unknown datatype of array\n");
				}
			}
        
		}       


		thisItem->itemLength = itemLength;
		gdfItemList->push_back(thisItem);


	}

	screenList = new ScreenList;
	toutList = new ScreenList;
	allList = new ScreenList;
	ScreenList * thisList;
	
	for(GdfItemList::iterator j=gdfItemList->begin(); j != gdfItemList->end(); j++)
	{
		GdfItem * thisItem = *j;

// 		printf("name %s\n", thisItem->name);
		const char ** cols =0;
		int ncols = 0;
		int type =0;
		int  T;	// column that containes the time
		if (! strcmp("position", thisItem->name) )
		{
			cols = columnp; ncols = 15;
			thisList= screenList;
			type =1;
			T=9;
		}
		if (! strcmp("time", thisItem->name) )
		{
			// first we assume that it is a snapshot, when we assemble the columns we may revise and say it is a tout
			cols = columns; ncols = 14;
			thisList= toutList;
			type =2;
			T=-1;
		}

// //                                0    1    2     3     4     5     6     7    8    9   10   11      12      13        14    15						screen
// const char * columnp[] = {"position", "x", "y", "z", "rxy", "Bx", "By", "Bz", "G", "t", "m", "q", "nmacro", "rmacro", "ID", "empty"};
// //                            0    1    2    3    4    5     6     7      8     9   10       11      12      13     14   						snapshot
// const char * columns[] = {"time", "x", "y", "z", "G", "Bx", "By", "Bz", "rxy", "m", "q", "nmacro", "rmacro", "ID", "empty"};
// //                            0    1    2    3    4    5     6     7      8     9   10       11      12      13     14    15      16     17     18     10     20	tout
// const char * columnt[] = {"time", "x", "y", "z", "G", "Bx", "By", "Bz", "rxy", "m", "q", "nmacro", "rmacro", "ID", "fEx", "fEy", "fEz", "fBx", "fBy", "fBz", "empty"};
		if(type)
		{
			Screen * thisScreen = new Screen;
			thisScreen->parameter = thisItem->dvalue;	// parameter
			thisScreen->avgTime   = thisItem->dvalue;	// will be overwriten later for position blocks

			thisScreen->type=type;
// 			printf("found %s at %e\n", thisItem->name, thisScreen->parameter);
			for(int i=1; i< ncols+1; i++)
			{
				j++;
				if( j == gdfItemList->end() )
				{
					printf("broken list\n");
					exit(-1);
				}
				GdfItem * thisItem = *j;
				if( type == 2 && !strcmp(cols[i], "empty") && !strcmp("fEx", thisItem->name) ) // yup, it is a tout and not a snapshot
				{
					cols = columnt; ncols = 20;
					type =3;
				}
				if( strcmp(cols[i], thisItem->name) )
				{
						fprintf(stderr, "name  mismatch position  or time = %e  %s %s\n", thisScreen->parameter, cols[i], thisItem->name);
						exit(-1);
				}
				if( !strcmp("empty", thisItem->name) ) break;

				thisScreen->dcord[i-1] = thisItem->darray;	// y
// 				printf("%3d  %s\n", i, thisItem->name); 


				if(i ==  1)	// x
				{
					thisScreen->npart = thisItem->block_arrsize;
				}
				else
				{
					if(thisScreen->npart != thisItem->block_arrsize)
					{
							fprintf(stderr, "num of particles mismatch position  or time = %e  %d %d\n", thisScreen->parameter, thisScreen->npart,   thisItem->block_arrsize);
							exit(-1);
					}
				}

				if(i == 3)	// z
				{
					thisScreen->avgPosition =0.;
					for(int i=0; i < thisItem->block_arrsize; i++)
						thisScreen->avgPosition += thisItem->darray[i];
					thisScreen->avgPosition /= thisItem->block_arrsize;
				}

				if(i == T)
				{
					thisScreen->avgTime =0.;
					for(int i=0; i < thisItem->block_arrsize; i++)
						thisScreen->avgTime += thisItem->darray[i];
					thisScreen->avgTime /= thisItem->block_arrsize;
				}
			}
			if(thisScreen->npart)
			{
				thisList->push_back(thisScreen);
				allList->push_back(thisScreen);
			}
			else
			{
				delete thisScreen;
			}
		}
	}
			

	std::sort(screenList->begin(), screenList->end(), Screen::screenSort);
	std::sort(toutList->begin(), toutList->end(), Screen::screenSort);
	std::sort(allList->begin(),    allList->end(),    Screen::allSort);

}



double GDF::getDouble(const char * name)
{
	for(GdfItemList::iterator j=gdfItemList->begin(); j != gdfItemList->end(); j++)
	{
		GdfItem * it = *j;
		if(  ! strcmp(it->name, name) ) return it->dvalue;
	}
	fprintf(stderr, "double %s not found in gdf file\n", name); 
	return 0.;
}

void GDF::addGdf(GDF * gdfs)
{
	// for simplicity this just joins the screenList, toutList and allList without deleting duplicates
	// the screen items and tout items are deleted from the gdfs->list, so that deleting all gdf objects
	// does not produce a memory leak. The pointers in the items will still point to gdfs->buff, so don't 
	// delete gdfs befor gdf.
	for(ScreenList::iterator j=gdfs->toutList->begin(); j != gdfs->toutList->end(); j++) 
	{
		toutList->push_back(*j);
		allList->push_back(*j);
		gdfs->toutList->erase(j);
	}
	for(ScreenList::iterator j=gdfs->screenList->begin(); j != gdfs->screenList->end(); j++) 
	{
		screenList->push_back(*j);
		allList->push_back(*j);
		gdfs->screenList->erase(j);
	}
	std::sort(screenList->begin(), screenList->end(), Screen::screenSort);
	std::sort(toutList->begin(), toutList->end(), Screen::screenSort);
	std::sort(allList->begin(),    allList->end(),    Screen::allSort);
}

void GDF::saveGdf()
{
	fprintf(stderr, "saveGdf() is not implemented\n");
	exit(-1);
}


void GDF::writeTimeStep(const char * filename, double time)
{
        int fo = open(filename, O_CREAT | O_WRONLY | O_TRUNC, 0644);
        if(fo < 0)
        {
                perror(filename);
                exit(0);
        }

        ssize_t wsize = write(fo, buff, 16+2*GDFNAMELEN);	// copy header from input

	double dist=1e33;
	GdfItemList::iterator   jbest = gdfItemList->end();
	for(GdfItemList::iterator j=gdfItemList->begin(); j != gdfItemList->end(); j++)
	{
		GdfItem * g = *j;
		
		if( ! strcmp(g->name , "time"))
		{
			double d =  fabs(g->dvalue- time);
			if(d < dist)
			{
				dist = d; jbest = j;
			}
		}
	}
	if(jbest == gdfItemList->end())
	{
		fprintf(stderr, "writeTimeStep: time %e not found\n", time);
		exit(-1);
	}

	GdfItem * g = *jbest;
	printf("found time at %e\n", g->dvalue);
	if( !(g->block_type & t_dir))
	{
		fprintf(stderr, "time item is not start dir\n");
		exit(-1);
	}
	for(GdfItemList::iterator j=jbest; j != gdfItemList->end(); j++)
	{
		GdfItem * g = *j;
		wsize = write(fo, g->startPtr, g->itemLength);
		if( g->block_type & t_edir)
		{
			close(fo);
			printf("done writing\n");
			return;	// end folder  
		}
	}
	printf("did not find the end_dir\n");
	exit(-1);
}


Particles::Particles(Screen * thisScreen)
{
	dcord = 0;
	scord = 0;
	npart = 0;

	coolenergy = 1.6e6;

	int G, ID;
	if(thisScreen->type == 1) { G=7; ID = 13; } // position
	else
	if(thisScreen->type == 2 || thisScreen->type == 3) { G=3; ID = 12; } // snapshot or tout
	else 
	{
		fprintf(stderr, "wrong type  in Particles::Particles at position %e\n", thisScreen->parameter);
		return;
	}

	if(thisScreen->npart <= 0)
	{
		fprintf(stderr, "no good particles in Particles::Particles at position %e\n", thisScreen->parameter);
		return;
	}

// //                                0    1    2     3     4     5     6     7    8    9   10   11      12      13        14    15
// const char * columnp[] = {"position", "x", "y", "z", "rxy", "Bx", "By", "Bz", "G", "t", "m", "q", "nmacro", "rmacro", "ID", "empty"};
// //                            0        1    2    3    4      5     6     7      8     9   10       11      12      13     14    15      16     17     18     10     20
// const char * columnt[] = {"time",     "x", "y", "z", "G",   "Bx", "By", "Bz", "rxy", "m", "q", "nmacro", "rmacro", "ID", "fEx", "fEy", "fEz", "fBx", "fBy", "fBz", "empty"};


	npart = thisScreen->npart;
	parameter = thisScreen->parameter;
	avgPosition = thisScreen->avgPosition;
	avgTime = thisScreen->avgTime;
       	dcord          = new double*[npart];
       	scord          = new double*[npart];
       	for(int i=0; i< npart; i++) 
       	{
               	scord[i] = dcord[i]       = new double[8];
	}

// 	double bz =0;
// 	for(int i=0; i<npart; i++) bz += thisScreen->dcord[6][i];
// 	bz /= npart;
// 	printf("bz = %e %e\n", avgPosition, bz);



	for(int i=0; i<npart; i++)
	{
		dcord[i][0] =  thisScreen->dcord[ 0][i];					// X
		dcord[i][1] =  thisScreen->dcord[ 4][i]/thisScreen->dcord[6][i]; 		// X'
		dcord[i][2] =  thisScreen->dcord[ 1][i];					// y
		dcord[i][3] =  thisScreen->dcord[ 5][i]/thisScreen->dcord[6][i];		// y'
		if(thisScreen->type == 1) // position
// 			dcord[i][4] =   (thisScreen->dcord[ 8][i] - avgTime)* 703.75e6*360.;	// time in RF degrees 
			dcord[i][4] =  -(thisScreen->dcord[ 8][i] - avgTime) *thisScreen->dcord[6][i]*clight;	// z = (t -<t>)* beta*clight
// 			dcord[i][4] =  -thisScreen->dcord[ 8][i]*thisScreen->dcord[6][i]*clight;	// z = t * beta*clight
// 			dcord[i][4] =  thisScreen->dcord[ 8][i]*clight;	// z = t * beta*clight
		if(thisScreen->type == 2 || thisScreen->type == 3)
			dcord[i][4] =  thisScreen->dcord[ 2][i];				// z
		dcord[i][5] = (thisScreen->dcord[ G][i] - 1.) * elmass;				// E
		dcord[i][6] =  thisScreen->dcord[ID][i];					// ID
		dcord[i][7] = 0.;							// for sortZ function
	}

}

Particles::~Particles()
{
       	if(dcord) for(int i=0; i< npart; i++) 
       	{
		delete dcord[i];
       	}
	delete  dcord;
	delete  scord;	// dcord and scord are poiting to the same data, delete those only once
	dcord = scord =0;
}


void Particles::sortz(char what)
{
	// what == 'x'	sort x
	// what == 'y'	sort y
	// what == 'r'	sort radius
	// what == 't'	sort time
	// what == 'z'	sort z
	// what == 'p'	sort momentum
	// what == 'a'	sort abs(delta p over p)

	if( what == 'x')
	{
        	for(int i = 0; i < npart; i++)
        	{
			dcord[i][7] = dcord[i][0] ;
		}
	}
	else
	if( what == 'y')
	{
        	for(int i = 0; i < npart; i++)
        	{
			dcord[i][7] = dcord[i][2] ;
		}
	}
	else
	if( what == 'r')
	{
        	for(int i = 0; i < npart; i++)
        	{
			dcord[i][7] = dcord[i][0]*dcord[i][0] + dcord[i][2]*dcord[i][2] ;
		}
	}
	else
	if( what == 't')
	{
        	for(int i = 0; i < npart; i++)
        	{
			double xp2 = dcord[i][1]*dcord[i][1];
			double yp2 = dcord[i][3]*dcord[i][3];
			dcord[i][7] = xp2 > yp2 ? xp2 : yp2;
// 			dcord[i][7] = dcord[i][1]*dcord[i][1] + dcord[i][3]*dcord[i][3] ;
		}
	}
	else
	if( what == 'z')
	{
        	for(int i = 0; i < npart; i++)
        	{
			dcord[i][7] = -dcord[i][4] ;
		}
	}
	else
	if( what == 'p')
	{
        	for(int i = 0; i < npart; i++)
        	{
			dcord[i][7] = dcord[i][5] ;
		}
	}
	else
	if( what == 'a')
	{
		double pm = 0.;
        	for(int i = 0; i < npart; i++)
			pm += dcord[i][5];
		pm /= npart;
        	for(int i = 0; i < npart; i++)
        	{
			dcord[i][7] = fabs(dcord[i][5]-pm) ;
		}
	}
	std::sort(scord, scord+npart, Particles::scordSort);
// 	for(int i=0; i<npart; i++)
// 		printf("sort %e\n", scord[i][7]);
}












// this is a backup from before mixing energy and momentum
// gonna use that to clean up
void Particles::calcSigmaMatrix(double ** pt, int np, int print)
{
	for(int j=0; j<7; j++)
	{
		for(int k=0; k<7; k++) sigmaMat[j][k] =0.;
		maxPart[j]  = -1e33;
		minPart[j]  =  1e33;
	}



        for(int i = 0; i < np; i++)
        {
// 		double savept5 = pt[i][5];
//  		double gammaL= pt[i][5]/elmass +1.; 
//  		pt[i][5] = elmass*sqrt(gammaL*gammaL-1.);

		pt[i][6] = 1.;
		for(int j=0; j<7; j++)
		{
	    		if(pt[i][j] > maxPart[j]) maxPart[j] = pt[i][j];
	    		if(pt[i][j] < minPart[j]) minPart[j] = pt[i][j];

			for(int k=j; k<7; k++) sigmaMat[j][k] += pt[i][j] * pt[i][k];
		}
// 		pt[i][5]=savept5;
	}


	for(int j=0; j<7; j++)
	{
		absPart[j] = fmax(fabs(maxPart[j]), fabs(minPart[j]));
	}


	for(int j=0; j<7; j++)
	{
		for(int k=j; k<7; k++)
		{
			sigmaMat[j][k] /= np;
			if( j != k) sigmaMat[k][j] = sigmaMat[j][k];
		}
	}

	for(int j=0; j<7; j++)
	{
		for(int k=j; k<7; k++)
		{
			sigmaMat2[k][j] = sigmaMat2[j][k] = sigmaMat[j][k] - sigmaMat[j][6]*sigmaMat[k][6];
		}
	}



	if(print)
	{
		printf("sigmaMat = \n");
		for(int j=0; j<7; j++)
		{
			for(int k=0; k<7; k++) printf("%12.4e ", sigmaMat[j][k]);
			printf("\n");
		}
		printf("\n");


		printf("sigmaMat2 = \n");
		for(int j=0; j<7; j++)
		{
			for(int k=0; k<7; k++) printf("%12.4e ", sigmaMat2[j][k]);
			printf("\n");
		}
		printf("\n");
	}
}


void Particles::twiss(int print)
{
// 	printf("************* twiss **************");
// 	setresults(1e33);
	twiss(print, dcord, npart);
// 	printresults("_new");
// 	printf("************* twiss ------********");
}
// this is a backup from before mixing energy and momentum
// gonna use that to clean up
void Particles::twiss(int print, double ** pt, int np)
{
	double freq = 1e33;
	calcSigmaMatrix(pt,  np, print);
	xmax = maxPart[0];
	ymax = maxPart[2];
	zmax = maxPart[4]; // length in meters;
	xmin = minPart[0];
	ymin = minPart[2];
	zmin = minPart[4]; // length in meters;
	amax = absPart[0]; if(absPart[2] > absPart[0] ) amax = absPart[2];
	dl_max = fabs(maxPart[4] - sigmaMat[4][6]);
	dl_min = fabs(minPart[4] - sigmaMat[4][6]);
	dp_max = fabs(maxPart[5] - sigmaMat[5][6])/sigmaMat[5][6];
	dp_min = fabs(minPart[5] - sigmaMat[5][6])/sigmaMat[5][6];
	tmax = dl_max/(freq*1e6 *360.); // length in seconds
	tmin = dl_min/(freq*1e6 *360.); // length in seconds

	orbx = sigmaMat[0][6];
	orbxp= sigmaMat[1][6];
	orby = sigmaMat[2][6];
	orbyp= sigmaMat[3][6];
	phas = sigmaMat[4][6];
	ener = sigmaMat[5][6];


	double gamma= sigmaMat[6][5]/ elmass+1.;
	gammaE = gamma;
	double betagamma = sqrt(gamma*gamma-1.);
	double betaE = betagamma/gamma;
	double gg = sqrt(square(sigmaMat[6][5]/elmass)+1.);
	double Eavg = (gg-1.)*elmass;
	ener1=Eavg;

        disl   = sigmaMat2[4][4] == 0 ? 0 : sigmaMat2[0][4]/sigmaMat2[4][4];
        dislp  = sigmaMat2[4][4] == 0 ? 0 : sigmaMat2[1][4]/sigmaMat2[4][4];
        dispx  = sigmaMat2[5][5] == 0 ? 0 : sigmaMat2[0][5]/sigmaMat2[5][5];
        dispxp = sigmaMat2[5][5] == 0 ? 0 : sigmaMat2[1][5]/sigmaMat2[5][5];
        dispy  = sigmaMat2[5][5] == 0 ? 0 : sigmaMat2[2][5]/sigmaMat2[5][5];
        dispyp = sigmaMat2[5][5] == 0 ? 0 : sigmaMat2[3][5]/sigmaMat2[5][5];

	sigma_p=sqrt(sigmaMat2[5][5]) / sigmaMat[6][5];
	sigma_l=sqrt(sigmaMat2[4][4]); // length in meters
	sigma_t=sqrt(sigmaMat2[4][4])/(freq*1e6 *360.); // length in seconds
	blen = sqrt(sigmaMat2[4][4]);
	coolsigma=sqrt(sigmaMat[5][5] - 2.* sigmaMat[6][5]*coolenergy +  coolenergy*coolenergy)/coolenergy;

	temperature = sqrt(gamma * gamma * ( sigmaMat2[1][1] + sigmaMat2[3][3])  + sigmaMat2[5][5]/(sigmaMat[6][5]*sigmaMat[6][5]) );
	temperature = sqrt(gamma * gamma * ( sigmaMat2[1][1] + sigmaMat2[3][3])  + sigmaMat2[5][5]/(sigmaMat[6][5]*sigmaMat[6][5]) + 5e-6/40. + 25e-8);   // add RHIC proton temperature
	temperxp = gamma * sqrt(sigmaMat2[1][1]);
	temperyp = gamma * sqrt(sigmaMat2[3][3]);

	double tempMat[7][7];
	memcpy(tempMat, sigmaMat2, 7*7*sizeof(double));
	double emit6d = dminv( (double *) tempMat, 7);

	magnetization = betagamma*(sigmaMat2[0][3]-sigmaMat2[1][2]);

	emit4d = betagamma*sqrt(sqrt(
	      sigmaMat2[0][0]*sigmaMat2[1][1]*sigmaMat2[2][2]*sigmaMat2[3][3]
	   -  sigmaMat2[0][0]*sigmaMat2[1][1]*sigmaMat2[2][3]*sigmaMat2[2][3]
	   -  sigmaMat2[0][0]*sigmaMat2[1][2]*sigmaMat2[1][2]*sigmaMat2[3][3]
	   +2*sigmaMat2[0][0]*sigmaMat2[1][2]*sigmaMat2[1][3]*sigmaMat2[2][3]
	   -  sigmaMat2[0][0]*sigmaMat2[1][3]*sigmaMat2[1][3]*sigmaMat2[2][2]
	   -  sigmaMat2[0][1]*sigmaMat2[0][1]*sigmaMat2[2][2]*sigmaMat2[3][3]
	   +  sigmaMat2[0][1]*sigmaMat2[0][1]*sigmaMat2[2][3]*sigmaMat2[2][3]
	   +2*sigmaMat2[0][1]*sigmaMat2[1][2]*sigmaMat2[0][2]*sigmaMat2[3][3]
	   -2*sigmaMat2[0][1]*sigmaMat2[1][2]*sigmaMat2[0][3]*sigmaMat2[2][3]
	   -2*sigmaMat2[0][1]*sigmaMat2[1][3]*sigmaMat2[0][2]*sigmaMat2[2][3]
	   +2*sigmaMat2[0][1]*sigmaMat2[1][3]*sigmaMat2[0][3]*sigmaMat2[2][2]
	   -  sigmaMat2[1][1]*sigmaMat2[0][2]*sigmaMat2[0][2]*sigmaMat2[3][3]
	   +2*sigmaMat2[0][2]*sigmaMat2[1][1]*sigmaMat2[0][3]*sigmaMat2[2][3]
	   -2*sigmaMat2[0][2]*sigmaMat2[1][3]*sigmaMat2[0][3]*sigmaMat2[1][2]
	   -  sigmaMat2[1][1]*sigmaMat2[0][3]*sigmaMat2[0][3]*sigmaMat2[2][2]
	   +  sigmaMat2[0][3]*sigmaMat2[0][3]*sigmaMat2[1][2]*sigmaMat2[1][2]
	   +  sigmaMat2[0][2]*sigmaMat2[0][2]*sigmaMat2[1][3]*sigmaMat2[1][3]
	   ));
	
// 	emitz= sqrt( ll*dd -ld*ld )* dpm*1000.;
// 	emitz= sqrt( sigmaMat2[4][4]*sigmaMat2[5][5] -sigmaMat2[4][5]*sigmaMat2[4][5] )* sigmaMat[6][5]*1000.;
	emitz= sqrt( sigmaMat2[4][4]*sigmaMat2[5][5] -sigmaMat2[4][5]*sigmaMat2[4][5] )*1000.;
	double tempt =  sigmaMat2[1][1]+sigmaMat2[3][3];
	tempt *= sigmaMat[6][5]*sigmaMat[6][5]*1e6/0.511000;
	eps_dx= sqrt((sigmaMat2[1][5]*sigmaMat2[1][5]*sigmaMat2[1][1] - 2*sigmaMat2[0][1]*sigmaMat2[0][5]*sigmaMat2[1][5] + sigmaMat2[0][5]*sigmaMat2[0][5]*sigmaMat2[0][0])/sigmaMat2[5][5]);
	eps_lx= sqrt((sigmaMat2[1][4]*sigmaMat2[1][4]*sigmaMat2[1][1] - 2*sigmaMat2[0][1]*sigmaMat2[0][4]*sigmaMat2[1][4] + sigmaMat2[0][4]*sigmaMat2[0][4]*sigmaMat2[0][0])/sigmaMat2[4][4]);
// 	eps_dx= sqrt((ddxpp*ddxpp*xpxp - 2*xxp*ddxp*ddxpp + ddxp*ddxp*xx)/dpdp);
// 	eps_lx= sqrt((dlxpp*dlxpp*xpxp - 2*xxp*dlxp*dlxpp + dlxp*dlxp*xx)/dldl);
	eps_dx= sqrt((sigmaMat[1][5]*sigmaMat[1][5]*sigmaMat[1][1] - 2*sigmaMat[0][1]*sigmaMat[0][5]*sigmaMat[1][5] + sigmaMat[0][5]*sigmaMat[0][5]*sigmaMat[0][0])/sigmaMat[5][5]);
	eps_lx= sqrt((sigmaMat[1][4]*sigmaMat[1][4]*sigmaMat[1][1] - 2*sigmaMat[0][1]*sigmaMat[0][4]*sigmaMat[1][4] + sigmaMat[0][4]*sigmaMat[0][4]*sigmaMat[0][0])/sigmaMat[4][4]);


        epsx = sqrt(sigmaMat2[0][0]*sigmaMat2[1][1] - sigmaMat2[0][1]*sigmaMat2[0][1]);
        epsx_n = epsx*betagamma;
        betx =   (sigmaMat2[0][0])/epsx;
        alfx =  -(sigmaMat2[0][1])/epsx;
	envx = sqrt(sigmaMat2[0][0]);
	divx = sqrt(sigmaMat2[1][1]);
	slopex = sigmaMat2[0][1]/envx;	// d/ds sqrt(<x^2>) = <xx'>/sqrt(<x^2>)

        epsy = sqrt(sigmaMat2[2][2]*sigmaMat2[3][3] - sigmaMat2[2][3]*sigmaMat2[2][3]);
        epsy_n = epsy*betagamma;
        bety =   (sigmaMat2[2][2])/epsy;
        alfy =  -(sigmaMat2[2][3])/epsy;
	envy = sqrt(sigmaMat2[2][2]);
	divy = sqrt(sigmaMat2[3][3]);
	slopey = sigmaMat2[2][3]/envy;	// d/ds sqrt(<x^2>) = <xx'>/sqrt(<x^2>)

	roundness = fabs(envx -envy);
	roundnessp = fabs(sigmaMat2[0][1] - sigmaMat2[2][3]);

        epss = sqrt(sigmaMat2[4][4]*sigmaMat2[5][5] - sigmaMat2[4][5]*sigmaMat2[4][5] )/sigmaMat[6][5]; // length in meters;
	double eprod = 0.5*epsx*epsy;
	if(print & 1) 
	{
// 	printf("%7s %10.4f %11.4e %11.4e  %11.4e   %10.4f %10.4f %10.4f  %10.4f %10.4f  %10.4f %10.4f %10.4f\n",
// 			type,avgPosition,emit4d,eprod, emit6d,betx,alfx,epsx,dispx, dispxp, bety,alfy,epsy);
// 			type,avgPosition,emit4d,eprod, ener,betx,alfx,epsx,dispx, dispxp, bety,alfy,epsy);
//
	}
	if(print & 2)
	{
		printf("npart              %d\n", npart);
		printf("avgPosition        %15.7f\n", avgPosition);
		printf("avgTime            %15.7e\n", avgTime);
		printf("orbit              %15.7e     %15.7e     %15.7e    %15.7e\n", sigmaMat[6][0], sigmaMat[6][1], sigmaMat[6][2], sigmaMat[6][3]);
		printf("dispersion         %15.7e     %15.7e     %15.7e    %15.7e\n", dispx, dispxp, dispy, dispyp);
		printf("emittance x y 4 xy %15.7e     %15.7e     %15.7e    %15.7e\n",   epsx,epsy,emit4d/betagamma,eprod);
		printf("normalizd x y 4 xy %15.7e     %15.7e     %15.7e    %15.7e\n",   epsx_n,epsy_n,emit4d,eprod*betagamma);
		printf("e longitudinal     %15.7e\n", epss);
		printf("gamma,beta         %15.7e     %15.7e\n", gamma, sqrt(1-1/(gamma*gamma)));
		printf("disl               %15.7e\n", disl);
		printf("dislp              %15.7e\n", dislp);
		printf("dldl               %15.7e\n", sigmaMat2[4][4]);
		printf("dpdp               %15.7e\n", sigmaMat2[5][5]);
		printf("dpAvg              %15.7e\n", sigmaMat[6][5]);
		printf("phasAvg            %15.7e\n", sigmaMat[6][4]);
		printf("coolsigma          %15.7e\n", coolsigma);
		printf("sigma_p            %15.7e\n", sigma_p);
		printf("sigma_l            %15.7e\n", sigma_l);
		printf("xminmax            %f %f\n", xmin, xmax);
		printf("yminmax            %f %f\n", ymin, ymax);
		printf("zminmax            %f %f\n", zmin, zmax);
		printf("xx, xpxp, yy, ypyp, tempt = %f %f %f %f %f\n", sigmaMat2[0][0], sigmaMat2[1][1], sigmaMat2[2][2], sigmaMat2[3][3], tempt);
		printf("sqrt(xx), sqrt(yy) = %f %f gamma=%f beta=%f\n", sqrt(sigmaMat2[0][0]), sqrt(sigmaMat2[2][2]), gamma, sqrt(1-1/(gamma*gamma)));
		printf("<x>, <y> = %f %f gamma=%f\n", sigmaMat[6][0], sigmaMat[6][2], gamma);
		printf("eps_dx, eps_dl = %e, %e\n", eps_dx, eps_lx);


		printf("ener, = %f %f %f %f\n", ener,ener,ener,ener);
		printf("betx,alfx,epsx = %f %f %e \n", betx,alfx,epsx);
		printf("bety,alfy,epsy = %f %f %e \n", bety,alfy,epsy);
		printf("sigma\t %15.7e   %15.7e  %15.7e  %15.7e\n", sigmaMat2[0][0],   sigmaMat2[0][1],  sigmaMat2[0][2],   sigmaMat2[0][3]  );
		printf("sigma\t %15.7e   %15.7e  %15.7e  %15.7e\n", sigmaMat2[0][1],  sigmaMat2[1][1], sigmaMat2[1][2],  sigmaMat2[1][3] );
		printf("sigma\t %15.7e   %15.7e  %15.7e  %15.7e\n", sigmaMat2[0][2],   sigmaMat2[1][2],  sigmaMat2[2][2],   sigmaMat2[2][3]  );
		printf("sigma\t %15.7e   %15.7e  %15.7e  %15.7e\n", sigmaMat2[0][3],  sigmaMat2[1][3], sigmaMat2[2][3],  sigmaMat2[3][3] );
		printf("\nM=%e\n",  (sigmaMat2[0][3]-sigmaMat2[1][2])*betagamma);
	}
}


double Particles::getValue(const char * what)
{
	if(!strcmp(what,"ener")) return ener;
	if(!strcmp(what,"gammaE")) return gammaE;
	if(!strcmp(what,"phas")) return phas;
	if(!strcmp(what,"blen")) return blen;
	if(!strcmp(what,"disl")) return disl;
	if(!strcmp(what,"dislp"))  return dislp;
	if(!strcmp(what,"dispx"))   return dispx;
	if(!strcmp(what,"dispxp"))  return dispxp;
	if(!strcmp(what,"dispy"))   return dispy;
	if(!strcmp(what,"dispyp"))  return dispyp;
	if(!strcmp(what,"disp+"))   return s_disp[slices-2];
	if(!strcmp(what,"dispp+"))  return s_dispp[slices-2];
	if(!strcmp(what,"disp0"))   return s_disp[(slices+1)/2];
	if(!strcmp(what,"dispp0"))  return s_dispp[(slices+1)/2];
	if(!strcmp(what,"disp-"))   return s_disp[1];
	if(!strcmp(what,"dispp-"))  return s_dispp[1];
	if(!strcmp(what,"orbx"))   return orbx;
	if(!strcmp(what,"orby"))   return orby;
	if(!strcmp(what,"orbxp"))   return orbxp;
	if(!strcmp(what,"orbyp"))   return orbyp;
	if(!strcmp(what,"epsx"))   return epsx;
	if(!strcmp(what,"epsx_n")) return epsx_n;
	if(!strcmp(what,"betx")) return betx;
	if(!strcmp(what,"alfx")) return alfx;
	if(!strcmp(what,"epsy")) return epsy;
	if(!strcmp(what,"epsy_n")) return epsy_n;
	if(!strcmp(what,"bety")) return bety;
	if(!strcmp(what,"alfy")) return alfy;
	if(!strcmp(what,"epss"))   return epss;
	if(!strcmp(what,"epssdl"))   return epss/sigma_l; 	// this is the enegy spread after removing the chirp
	if(!strcmp(what,"magnet")) return magnetization;
	if(!strcmp(what,"emit4d")) return emit4d;
	if(!strcmp(what,"emit4dsum")) return emit4dsum;
	if(!strcmp(what,"emitz")) return emitz;
	if(!strcmp(what,"temper")) return temperature;
	if(!strcmp(what,"temperx")) return temperxp;
	if(!strcmp(what,"tempery")) return temperyp;
	if(!strcmp(what,"ttrans")) return ttrans;
	if(!strcmp(what,"ttranm")) return ttranm;
	if(!strcmp(what,"r_eps")) return r_eps;
	if(!strcmp(what,"r_bet")) return r_bet;
	if(!strcmp(what,"r_alf")) return r_alf;
	if(!strcmp(what,"r_alfObet")) return r_alfObet;
	if(!strcmp(what,"r_rot")) return r_rot;
	if(!strcmp(what,"epsx_c")) return epsx_c;
	if(!strcmp(what,"epsy_c")) return epsy_c;
	if(!strcmp(what,"epsx_u")) return epsx_u;
	if(!strcmp(what,"epsy_u")) return epsy_u;
	if(!strcmp(what,"eps_dx")) return eps_dx;
	if(!strcmp(what,"eps_lx")) return eps_lx;
	if(!strcmp(what,"rl2")) return rl2;
	if(!strcmp(what,"rot")) return rot;
	if(!strcmp(what,"envx")) return envx;
	if(!strcmp(what,"envy")) return envy;
	if(!strcmp(what,"divx")) return divx;
	if(!strcmp(what,"divy")) return divy;
	if(!strcmp(what,"round")) return roundness;
	if(!strcmp(what,"roundp")) return roundnessp;
	if(!strcmp(what,"slopex")) return slopex;
	if(!strcmp(what,"slopey")) return slopey;
	if(!strcmp(what,"amax")) return amax;
	if(!strcmp(what,"xmax")) return xmax;
	if(!strcmp(what,"ymax")) return ymax;
	if(!strcmp(what,"zmax")) return zmax;
	if(!strcmp(what,"tmax")) return tmax;
	if(!strcmp(what,"xmin")) return xmin;
	if(!strcmp(what,"ymin")) return ymin;
	if(!strcmp(what,"zmin")) return zmin;
	if(!strcmp(what,"tmin")) return tmin;
	if(!strcmp(what,"sigma_p")) return sigma_p;
	if(!strcmp(what,"coolsigma")) return coolsigma;
	if(!strcmp(what,"dp_max")) return dp_max;
	if(!strcmp(what,"dp_min")) return dp_min;
	if(!strcmp(what,"sigma_l")) return sigma_l;
	if(!strcmp(what,"sigma_t")) return sigma_t;
	if(!strcmp(what,"dl_max")) return dl_max;
	if(!strcmp(what,"dl_min")) return dl_min;
	if(!strcmp(what,"epss_dl")) return epss/sigma_l;



	fprintf(stderr, "parameter %s is not defined in Pmla::getValue\n", what);
	printf(         "parameter %s is not defined in Pmla::getValue\n", what);
	exit(-1);
}

int Particles::getCurve(const char * what)
{
	if(!strcmp(what,"ener")) 	return 0;
	if(!strcmp(what,"gammaE")) 	return 0;
	if(!strcmp(what,"phas")) 	return 0;
	if(!strcmp(what,"blen")) 	return 0;
	if(!strcmp(what,"disl")) 	return 0;
	if(!strcmp(what,"dislp"))  	return 0;
	if(!strcmp(what,"dispx"))   	return 0;
	if(!strcmp(what,"dispxp"))  	return 0;
	if(!strcmp(what,"dispy"))   	return 0;
	if(!strcmp(what,"dispyp"))  	return 0;
	if(!strcmp(what,"disp+"))   	return 0;
	if(!strcmp(what,"dispp+"))  	return 0;
	if(!strcmp(what,"disp0"))   	return 0;
	if(!strcmp(what,"dispp0"))  	return 0;
	if(!strcmp(what,"disp-"))   	return 0;
	if(!strcmp(what,"dispp-"))  	return 0;
	if(!strcmp(what,"orbx"))   	return 0;
	if(!strcmp(what,"orby"))   	return 0;
	if(!strcmp(what,"orbxp"))   	return 0;
	if(!strcmp(what,"orbyp"))   	return 0;
	if(!strcmp(what,"epsx"))   	return 1;
	if(!strcmp(what,"epsx_n")) 	return 1;
	if(!strcmp(what,"betx")) 	return 0;
	if(!strcmp(what,"alfx")) 	return 0;
	if(!strcmp(what,"epsy")) 	return 1;
	if(!strcmp(what,"epsy_n")) 	return 1;
	if(!strcmp(what,"bety")) 	return 0;
	if(!strcmp(what,"alfy")) 	return 0;
	if(!strcmp(what,"epss"))   	return 1;
	if(!strcmp(what,"epssdl"))   	return 1; 	// this is the enegy spread after removing the chirp
	if(!strcmp(what,"magnet")) 	return 0;
	if(!strcmp(what,"emit4d")) 	return 1;
	if(!strcmp(what,"emit4dsum")) 	return 0;
	if(!strcmp(what,"emitz")) 	return 1;
	if(!strcmp(what,"temper")) 	return 1;
	if(!strcmp(what,"temperx")) 	return 1;
	if(!strcmp(what,"tempery")) 	return 1;
	if(!strcmp(what,"ttrans")) 	return 1;
	if(!strcmp(what,"ttranm")) 	return 0;
	if(!strcmp(what,"r_eps")) 	return 0;
	if(!strcmp(what,"r_bet")) 	return 0;
	if(!strcmp(what,"r_alf")) 	return 0;
	if(!strcmp(what,"r_alfObet")) 	return 0;
	if(!strcmp(what,"r_rot")) 	return 0;
	if(!strcmp(what,"epsx_c")) 	return 0;
	if(!strcmp(what,"epsy_c")) 	return 0;
	if(!strcmp(what,"epsx_u")) 	return 0;
	if(!strcmp(what,"epsy_u")) 	return 0;
	if(!strcmp(what,"eps_dx")) 	return 0;
	if(!strcmp(what,"eps_lx")) 	return 0;
	if(!strcmp(what,"rl2")) 	return 0;
	if(!strcmp(what,"rot")) 	return 0;
	if(!strcmp(what,"envx")) 	return 0;
	if(!strcmp(what,"envy")) 	return 0;
	if(!strcmp(what,"divx")) 	return 0;
	if(!strcmp(what,"divy")) 	return 0;
	if(!strcmp(what,"round")) 	return 0;
	if(!strcmp(what,"roundp")) 	return 0;
	if(!strcmp(what,"slopex")) 	return 0;
	if(!strcmp(what,"slopey")) 	return 0;
	if(!strcmp(what,"amax")) 	return 0;
	if(!strcmp(what,"xmax")) 	return 0;
	if(!strcmp(what,"ymax")) 	return 0;
	if(!strcmp(what,"zmax")) 	return 0;
	if(!strcmp(what,"tmax")) 	return 0;
	if(!strcmp(what,"xmin")) 	return 0;
	if(!strcmp(what,"ymin")) 	return 0;
	if(!strcmp(what,"zmin")) 	return 0;
	if(!strcmp(what,"tmin")) 	return 0;
	if(!strcmp(what,"sigma_p")) 	return 1;
	if(!strcmp(what,"coolsigma")) 	return 1;
	if(!strcmp(what,"dp_max")) 	return 0;
	if(!strcmp(what,"dp_min")) 	return 0;
	if(!strcmp(what,"sigma_l")) 	return 1;
	if(!strcmp(what,"sigma_t")) 	return 1;
	if(!strcmp(what,"dl_max")) 	return 0;
	if(!strcmp(what,"dl_min")) 	return 0;
	if(!strcmp(what,"epss_dl")) 	return 1;



	fprintf(stderr, "parameter %s is not defined in Pmla::getValue\n", what);
	return -1;
}


