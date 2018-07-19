//	readgdf program   J. Kewisch  2016
//
//
// 	the gdf file consists of a header and a list of data items
// 	each data item has:
// 		a name (not unique)
// 		flags telling the data type, and a start/stop for sub groups
// 		and [optionally] the data
// 	That's pretty much all.
//
//	The constructor of the GDF class reads the file and create the list of items (gdfItemList)
//
// 	If the start_dir flag of an item is set all following  data items are sub-data to this item until the end_dir flag is set.
//	The output of the GPT program has the "position" and "time" items, that have a list of sub-items cointaining the the coordinates 
//	of the particles. This information is collected from the gdfItemList and stored in the screenList and toutList.
// 	Both lists contain pointers to the original buff memory, so don't delete buff.
//
//	For LEReC we need to filter which electrons we will count to calculate the beam quality. In screenList and toutList the particles
//	are organized in columns. To filter we copy the data of one snapshot into the dcords array of a Particle object, which is organized in rows.
//	Then we use the scord pointers to sort the electrons according to criteria  using sortz(). 
//
//	the function twiss() may then be used to calculate opticla functions.








#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <error.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <vector>
#include <algorithm>



const double elmass = 0.51104e6;	// eV
const double clight=2.997824562e8;      // m/s


inline double square(double x) { return x*x;}
double dminv( double * m6, int n);


	typedef unsigned int uint;
	typedef unsigned char uchar;

	static const int GDFID = 94325877;
	static const int GDFNAMELEN = 16; 

// Data types 
	static const  int t_ascii    = 0x001;      //  Char array
	static const  int t_s32      = 0x002;      //  Signed long           
	static const  int t_dbl      = 0x003;      //  Double    
	static const  int t_nul      = 0x010;      //  No data               
	
	static const  int t_undef    = 0x000;      //  Data type not defined 
	static const  int t_u8       = 0x020;      //  Unsigned char         
	static const  int t_s8       = 0x030;      //  Signed char           
	static const  int t_u16      = 0x040;      //  Unsigned short        
	static const  int t_s16      = 0x050;      //  Signed short          
	static const  int t_u32      = 0x060;      //  Unsigned long         
	static const  int t_u64      = 0x070;      //  Unsigned 64bit int    
	static const  int t_s64      = 0x080;      //  Signed 64bit int      
	static const  int t_flt      = 0x090;      //  Float                 

// Block types 
	static const  int t_dir      = 0x100;      //  Directory entry start
	static const  int t_edir     = 0x200;      //  Directory entry end 
	static const  int t_param    = 0x400;      //  Parameter        
	static const  int t_data     = 0x800;      //  Data array


class GdfItem
{
public:
	char * startPtr;
	int itemLength;
	char name[GDFNAMELEN];
	int ivalue;		// could have use a union here, but memory is cheap
	float fvalue;
	double dvalue;
	const char * charvalue;
	double * darray;
	float * farray;
	int * iarray;

	int block_type;	
	int block_size;		// array size in bytes
	int block_arrsize;	// array size


	GdfItem()
	{
		block_type=0;
		block_size=0;
		block_arrsize=0;
		ivalue=0; 
		fvalue=0.;
		dvalue=0.;
		charvalue=0;
	}

};

typedef std::vector<GdfItem*> GdfItemList;



class Screen	// used for position and time blocks
{
public:
	double parameter; 	// the value of position or time
	double avgPosition;	// the average position, calculated from all particles; for screen this should be identical with parameter
	double avgTime;		// the average time, calculated from all particles; for tout this is identical with parameter
	int type;	// 1 = tout, 2 = screen, 3 = mr
	int npart; 
	double * dcord[21];
	static bool screenSort(Screen * an, Screen * am)
	{
        	return( an-> parameter < am-> parameter) ;
	}
	static bool allSort(Screen * an, Screen * am)
	{
        	return( an-> avgPosition < am-> avgPosition) ;
	}
	Screen()
	{
		npart=0;
		type=0;
	}
};
typedef std::vector<Screen*> ScreenList;




class GDF
{

	int count;
	uchar readUChar(char  * &p);
	int readInt(char  * &p);
	uint readUInt(char  * &p);
	float readFloat(char  * &p);
	double  readDouble(char  * &p);
        char * buff;
        off_t size;
public:
	GDF(const char * filename);
	virtual ~GDF();
	ScreenList * screenList;
	ScreenList * toutList;
	ScreenList * allList;
	GdfItemList * gdfItemList;
	
	double getDouble(const char * name);
	void writeTimeStep(const char * filename, double time);
	void addGdf(GDF * gdfs);
	void saveGdf();
};

class Particles
{
	static bool scordSort(double  * an, double  * am)
	{
        	return( an[7] < am[7]) ;
	}
public:
        double ** dcord;
        double ** scord;
	int npart;
	double parameter;
	double avgPosition;
	double avgTime;
	double coolenergy;
//  results
        double ener;		// energy avarage
        double ener1;		// energy avarage (test)
        double gammaE;		// gamma avarage 
        double phas;		// phase avarage
        double blen;		// buch length rms
        double disl;		// horizontal dispersion of element;
        double dislp;		// horizontal slope of dispersion at element
        double dispx;		// horizontal dispersion of element;
        double dispxp;		// horizontal slope of dispersion at element
        double dispy;		// vertical dispersion of element;
        double dispyp;		// vertical slope of dispersion at element
        double orbx;		// avarage orbit
        double orbxp;		// avarage orbit
        double orby;		// avarage orbit
        double orbyp;		// avarage orbit
        double epsx;		// horizontal emittance at element
        double epsx_n;		// normalized horizontal emittance at element
        double betx;		// horizontal beta function at element
        double alfx;		// horizontal alpha function at element
        double epsy;		// vertical emittance function at element
        double epsy_n;		// normalized vertical emittance at element
        double bety;		// vertical beta function at element
        double alfy;		// vertical alpha function at element
        double epss;		// longitudinal emittance function at element
        double magnetization;	// <beta*gamma*(xy'-yx')>
        double emit4d;		// 4 d emittance
        double emit4dsum;	// sum of 4 d slice emittance
	double emitz;             // longitudinal emittance in  deg*keV
	double  temperature;	//  sqrt(gamma^2*x'^2 + gamma^2*y'^2 + dp/p^2)
	double  temperxp;	//  sqrt(gamma^2*x'^2)
	double  temperyp;	//  sqrt(gamma^2*y'^2)
        double ttrans;		// transverse temperature
        double ttranm;		// transverse temperature (magnetized)
	int slices;
	double r_eps;		// radial emittance
	double r_bet;		// radial beta
	double r_alf;		// radial alfa
	double r_alfObet;	// radial alfa/beta
	double r_rot;		// radial alfa
	double s_eps[500];	// slice emittace
	double s_bet[500];
	double s_alf[500];
	double s_disp[500];
	double s_dispp[500];
	double s_rot[500];	// slice vortex
        double epsx_c;		// horizontal emittance at elementi after magnetization correction
        double epsy_c;		// horizontal emittance at elementi after magnetization correction
        double epsx_u;		// horizontal emittance at elementi after unrotation
        double epsy_u;		// horizontal emittance at elementi after unrotation
        double eps_dx;		// dispersion contribution to the emittance
        double eps_lx;		// longitudinal position contribution to the emittance
        double rl2;		// avarage square of larmour radii
        double rot;		// rotation speed of the beam = (xy'-yx')/r2
        double envx;		// x enveloppe
        double envy;		// y enveloppe
        double roundness;	// enveloppe roundness
        double divx;		// x' enveloppe
        double divy;		// y' enveloppe
        double roundnessp;	// slope roundness
        double slopex;		// x' enveloppe
        double slopey;		// y' enveloppe
	double amax;
	double xmax;
	double ymax;
	double zmax;
	double tmax;
	double xmin;
	double ymin;
	double zmin;
	double tmin;
        double sigma_p;		// energy spread
        double coolsigma;	// sigma p relative to coolenergy
        double dp_max;		// energy spread
        double dp_min;		// energy spread
        double sigma_l;		// bunch length
        double sigma_t;		// bunch length in seconds
        double dz_max;		// bunch length
        double dz_min;		// bunch length
        double dl_max;		// 
        double dl_min;		// 




// matrix

	double sigmaMat[7][7];		// the raw sigma matrix: the columns are x,x',y,y',phase, momentum and 1 (one), colume 6 containes the averages
	double sigmaMat2[7][7];		// the sigma matrix with the avarage removed, colume 6 containes then zero
	double maxPart[7], minPart[7], absPart[7];

	Particles(Screen * thisScreen);
	virtual ~Particles();
	void sortz(char what);
	void twiss(int print);
	void twiss(int print, double ** pt, int np);
	void calcSigmaMatrix(double ** pt, int np, int print);
	double getValue(const char * what);
	static int getCurve(const char * what);
};


