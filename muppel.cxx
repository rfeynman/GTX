


// Notice: This computer software was prepared by Brookhaven Science
// Associates, (the Contractor), under Federal Contract No. DE-AC02-98CH10886
// (the Contract) with the U.S. Department of Energy (DOE).  All rights in
// this computer software are reserved by DOE on behalf of the United States
// Government and the Contractor as provided in the Contract.  You are
// authorized to use this computer software solely for U.S. Governmental
// purposes.  This software is not to be released or distributed in any form
// or manner except as may be provided in your Software User Acknowledgment.
// NEITHER THE UNITED STATES GOVERNMENT NOR THE CONTRACTOR MAKES ANY
// WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
// THIS SOFTWARE.  This notice including this sentence must appear on any
// copies of this computer software.


// ===================================================================================================================
//
//       GPT FITTING PROGRAM
//
//       Jorg Kewisch, 2016
//
//      syntax:	muppel [-r] [-d] [-p] <GPT input file> 
//      flags:	-r	apply the recovery file. This allows to continue where the previous run ended
//      	-d	start the Xvfb program with display :60 or :61
//      	-p	run GPT only once
//
// 	The program uses all cores of the computer. It finds the number of cores from the /prog/cpuinfo file with the line:
// 	int maxthreads = of->setMaxThreads(0);
// 	On windows systems you have to chage thezero to the number of cores you have.
// 	The program creates a directory "calcdir" in which the calculations are done. Whenever a "better" solution is found
// 	the input file is copied to gpt_best.gpt and the  output files are copied back into the current directory.
// 
//      the input file is a GPT input file that contains the fitting instructions in GPT comment lines
//	the fitting package CONDOR by Ir. Frank Vanden Berghen, Univesite' Libre de Bruxelles is used for fitting
//
//	The optimizer uses instructions starting with #@, which are comments in gpt. These instructions can easily be disabled by changing #@ to #x. 
//	The beam line can be calculated in sections. The optimizer adds a section=# to the gpt command. The
//	#@section <first> <last> ...	defines which sections of gpt should be run. 
//
//	#@var  <parameter name> <initial step size> <lower bound> <upper bound>  	defines a fit parameter. parameter number starts with zero. 
//              <active> equal 1 means this variable is used, 0 means it is ignored
//		<initial step size> guess, make it large, but small enough so that the function is stable +- 2 steps away
//		<lower bound> <upper bound>  are boundaries for the step, not the final value 
//	the parameter must be defined in gpt in the form
//	p=5.567; ....
//	at the begin of a line (ignoring whitespace).
//	the optimizer will replace the line up to the semicolon with a new value.
// 	the optimizer will remember the line number where the parameter is found.
//	Changing the input file while the optimizer is running is not god.
//	if there are multiple lines the optimizer will  use the last ocurence.
//	The optimizer is not aware of any if-then-else constructs.
//
//	#@goal <where> <cut what> <cut howmuch> <what> <goal> <weight> <goodEnough> [ <what> <goal> <weight> <goodEnough> ....]
//	defines the figure of merrit. any number of goals can be used
//		<where> is the name of a gpt variable which has to be written into the GPT output file using the output() function.
//		example:
//		here1 = 15.43;
//		output("here1", here1);
//		screen("wcs", "I", here1 , "wcs");
// 		#@goal here1 sigma_p 0 1 eps_x 0 100 0
// 		There can be mant quaruples of <what> <goal> <weight> <goodEnough> on the line.
//		<what> describes what to optimize. see the getval function in gpt.cxx
//		<goal> is the wanted value for <what>, normally zero.
//		<weight> is the weight of a goal when more than one goal is fitted. One is a good weight.
//		<cut what> is 'p' or 'l' or 'n'
//		<cut howmuch> how mant particle to ignore in the evaluation. 0.3 means that 30 % will be ignored
//
//	Other #@ codes are ignored by this program
//
//
// ===================================================================================================================
//	the program can be aborted in a orderly manner by creating a file named "stop" in the current directory , i.e.  "touch stop"
//	the stop file will be removed o exit
//	the program can be paused by creating a file named "pause" in the current directory. when the user deletes "pause" the optimization continues.
//	A file "plot" will be removed after running ~/RTX/main5. "plot will be removed.
// ===================================================================================================================
//
//	this program calls the condor program, which is licensed under GPL. I guess that makes this programm GPL too.
//
//	condor builds a map using (n+1)*(n+2)/2 calculations and then uses this map to optimize




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>
#include <time.h>
#include <pthread.h>
#include <dirent.h>
#include <time.h>
#include <ctype.h>


#include "gdf.hxx"



#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "Solver.hxx"
#include "tools.hxx"
#include "ObjectiveFunction.hxx"

#ifdef WIN32
#include <crtdbg.h>
#endif

pthread_mutex_t rtxMutex= PTHREAD_MUTEX_INITIALIZER;

// const double clight= 2.99792456210e8;	// in m/s
// const double elmass = 0.51104;
// inline double square(double x) { return x*x;}

const int maxFitVariables = 100;		// max number of fit parameters
const int maxFitGoals = 100;		// max number of fit goals
const int maxSwitches = 100;		// max number of fit goals
const int maxRunStore = 1000;	    	// max number of states in the recobvery buffer;

int whatsystem;
char * ParmelaDir;
char * Wine;
char * preParmela;
char * slenvProg;
char * main5;
pid_t pidXvfb;
const char * thislock ="aaaRTXmainpRuns";

int myexit(int err)
{
	if(pidXvfb > 0) kill(pidXvfb, SIGTERM);
	unlink(thislock);
	fflush(stdout);
	fflush(stderr);
	exit(err);
}


int splitline(char * line, char ** space,  char ** word)
{
	int n=0;
	char * p = line;
	char newSpace[1000], newWord[1000];
	while(1)
	{
		char * q = newSpace;
		while( *p > 0 && *p <= ' ') {*q++=*p++; }
		*q=0;
		if(*p == 0) return n;
		q = newWord;
		while(*p > ' ') {*q++=*p++;}
		*q=0;
		space[n]=strdup(newSpace);
		word [n]=strdup(newWord );
		n++;
	}	
}
int str2words(char * line, char * beg[], int length[], char * after[])
{
	int n=0;
	char * p = line;
	while(1)
	{
		while( *p > 0 && *p <= ' ') p++;
		if(*p == 0) return n;
		beg[n] = p;
		int l=0;
		while(*p > ' ') { l++; p++;}
		length[n]=l;
		after[n]=p;
		n++;
	}	
}


class M2condor : public UnconstrainedObjectiveFunction
{
		int nrun;
		double para[maxRunStore][maxFitVariables], resu[maxRunStore], subresu[maxRunStore][maxFitVariables];
		double cath_size, cath_field;
		double cut;
		int cutWhat;
		int has_beta, has_cath, has_slenv;
		char slenvLatt[300];
		char slenvCavt[300];
		char switchName [100][20];
		int  switchFlag[100];
		int nswitches;
		int oneTime;



	public:
		int sections[50], nsections;
		char inputFileName[300];	// input file
		char fileroot[300];		// input file root
		char cwd[300];			// name of current directory, is used in temp file names
		char fit_what[maxFitGoals][10];		// name of what to optimize, i.e. emit4d..
		char fit_where[maxFitGoals][30];	// name of the GPT variable telling us where to look
		double fit_val[maxFitGoals];		// desired value
		double fit_weight[maxFitGoals];		// weight relative to other goals
		double fit_good[maxFitGoals];		// the "good enough" parameter
		char   fit_cut[maxFitGoals];		// what to cut ( n, l, p)
		double fit_cut_ratio[maxFitGoals];	// how much to cut (0.3 means keep 70 %)
		int ngoals;	// number of goals


		class Var
		{
			public:
			char   name[50];	// name of fit parameter
			double mult;	// initial step size for the optimizer
			double lower;	// boundary
			double upper;	// boundary
			double step;	// step size for the step function
			double initVal;	// value of var in input file
			int    nline;	// line # in input file where var is defined
			char * restOfLine; // rest of line in input file where var is defined
		} var[maxFitVariables];
		int nvar;	// number of active variables

		int ngood;	// the number of good particles, which must not get lost 
		int ninput;	// number of particles specified on "input" lines
		double bestgoal;

		FILE * flog;
		M2condor(char * inpfilename, int recover, int oneTime);
    		double eval(Vector v, int *nerror=NULL, int resource=0);
    		double eval();
		void done();
};




M2condor::M2condor(char * inpfilename, int recover, int oneTime) 
{
	this->oneTime = oneTime;
	time_t tim=time(0);
	char logfile[50];
	sprintf(logfile,"gpt_log.%ld.txt", tim);
	unlink("gpt_log.txt");
	symlink(logfile, "gpt_log.txt");
	flog = fopen(logfile, "w");
	if(!flog)
	{
		fprintf(stderr, "cant open gpt_log.txt\n");
		myexit(-1);
	}
	fprintf(flog, "Input file  %s,  pid %d \n", inpfilename, getpid());
       	FILE * fi = fopen(inpfilename, "r");
	if(!fi)
	{
		fprintf(stderr, "cant open %s\n", inpfilename);
		myexit(-1);
	}
	int ls=strlen(inpfilename);
	strcpy(fileroot, inpfilename);
	fileroot[ls-4]=0;
	strcpy(inputFileName, inpfilename);
	getcwd(cwd, 300);
	fprintf(flog, "dir is %s\n", cwd);
	char ctim[100];
	ctime_r(&tim, ctim);
	fprintf(flog, "time is %s\n", ctim);


	ngood = -1; 	// not set
	ninput =1;	// there is always the reference particle
	ngoals=0;
	nvar=0;
	nvar=0;
	nsections=0;
	nswitches = 0;
	cutWhat = 0;
	has_beta= has_cath=has_slenv=0;

	char line[500];
	char * pl;
	char what[100];
	int sofar;
	while( (pl=fgets(line,499,fi)))
	{
		if(line[0] == '#' && line[1] == '@')
		{
// 			printf("=%s\n", line);
			pl +=2;
			sscanf(pl,"%s%n", what, &sofar);
			pl += sofar;

			// #@var 0 booph   1e-2       -10000 1000 booster phase
			if( ! strcmp(what, "var") ) 
			{
				double p_mult, p_lower, p_upper;
				char varname[50];
				int st=sscanf(pl, "%s%lf%lf%lf", varname, &p_mult, &p_lower, &p_upper);
				if(st != 4)  {fprintf(stderr, "error in var statement !!!\n"); myexit(-1);}
				fprintf(flog, "parameter   %-20s %f %f %f\n", varname, p_mult, p_lower, p_upper);
				if(nvar >= maxFitVariables)  {fprintf(stderr, "too many variables !!!\n"); myexit(-1);}

				strcpy(var[nvar].name, varname);
				var[nvar].mult = p_mult;
				var[nvar].lower = p_lower;
				var[nvar].upper = p_upper;
				var[nvar].restOfLine = 0;
				var[nvar].nline = -1;
				nvar ++;
				continue;
			}


			if( ! strcmp(what, "goal") ) 
			{
				printf("%s", line);
				char where[30];
				char cut[10];   // should be n, l or p
				double howMuch;
				if(sscanf(pl, "%s%s%lf%n",  where, cut, &howMuch, &sofar) != 3)
				{
					fprintf(stderr, "misformed goal: %s\n", line);
					myexit(-1);
				}
				pl += sofar;
				printf("%s", pl);
				while(sscanf(pl, "%s%lf%lf%lf%n",  fit_what[ngoals], fit_val+ngoals, fit_weight+ngoals, fit_good+ngoals, &sofar) == 4)
				{
					pl += sofar;
					fit_cut[ngoals] = cut[0];
					fit_cut_ratio[ngoals] = howMuch;
					if(fit_good[ngoals] < 1e-10) fit_good[ngoals] = 1e-10;
					strcpy(fit_where[ngoals], where);
					int c = Particles::getCurve(fit_what[ngoals]);	// this calls exit() if fit_what is not valid
					if( c < 0) myexit(-1);
					ngoals++;
				}
				continue;
			}



			if( ! strcmp(what, "sections") ) 
			{
				int from, to;
				int n = sscanf(pl, "%d%d", &from, &to);
				if( n == 1) to = from;
				else if( n != 2) 
				{
					fprintf(stderr, "misformed sections: %s\n", line);
					myexit(-1);
				}

				nsections = to -from + 1;
				for(int n= 0; n<nsections; n++) sections[n]= from+n;
				continue;

			}


			fprintf(flog, "ignored %s", line);
		}

		what[0]=0;
		sscanf(line, "%s", what);
		for(char * cp = what; *cp; cp++) *cp = tolower(*cp);
		if( ! strcmp(what, "input")  ) 
		{
			int type, num;
			for(char * p = line; *p; p++) if(*p == ',') *p=' ';
			sscanf(line, "%s%d%d", what, &type, &num);
			ninput += num;
		}


	}
	fclose(fi);
	fprintf(flog, "found %d particles on input lines\n", ninput);

	if(!nsections)
	{
		nsections=1;
		sections[0]=0;
		fprintf(flog, "found no sctions, set section=0\n");
	}
	else
		fprintf(flog, "found %d sections\n", nsections);


	if(nvar <= 0 && (!oneTime))
	{
		fprintf(stderr, "nvar = %d, no parameters specified or used\n", nvar);
		myexit(-1);
	}

	// find the var statements
       	fi = fopen(inpfilename, "r");
	if(!fi)
	{
		fprintf(stderr, "cant open %s\n", inpfilename);
		myexit(-1);
	}
	int nline=0;
	while( (pl=fgets(line,499,fi)))
	{
		char line1[500];
		char *q=line1;
		for( char * p = line; *p; p++) if(*p > ' ') { *q =*p; q++;}
		*q=0;
		for(int ivar=0; ivar<nvar; ivar++)
		{
			int l = strlen(var[ivar].name);
			if(! strncmp(line1, var[ivar].name, l)  && line1[l] == '=')
			{
				double val;
				sscanf(line1+l+1, "%le", &val);
				var[ivar].initVal = val;	// value of var in input file
				free(var[ivar].restOfLine);
				for( char * p = line; *p; p++) if(*p == ';') {var[ivar].restOfLine = strdup(p); break;}
				var[ivar].nline = nline;
			}
		}
		nline++;
	}
	fclose(fi);
	// sort (only a few, stupid sort works just fine)
	for(int i=0;   i<nvar-1; i++)
	for(int j=i+1; j<nvar  ; j++)
	{
		if( var[i].nline > var[j].nline )
		{
			Var t=var[i]; var[i]=var[j]; var[j]=t;
		}
	}


	for(int ivar=0; ivar<nvar; ivar++)
	{
		if(var[ivar].nline < 0)
		{
			fprintf(flog, "variable %s definition not found\n", var[ivar].name);
			fprintf(stderr, "variable %s definition not found\n", var[ivar].name);
			myexit(-1);
		}
		fprintf(flog, "var %d  %-20s %15.7e %5d %s", ivar, var[ivar].name, var[ivar].initVal, var[ivar].nline, var[ivar].restOfLine);
	}

	fprintf(flog, "             where    what   val      weigh   col\n");
	for(int i=0; i< ngoals; i++)
		fprintf(flog, "condition  %s %s %f %f\n", fit_where[i], fit_what[i], fit_val[i], fit_weight[i]);










// ----------------------------------------
	nrun=0; 
	bestgoal=1e33;
	if(recover)
	{
		FILE * fp = fopen("recover.hex", "r");
		if(!fp)
		{
			fprintf(stderr,"cant open %s for read\n", "recover.hex");
			myexit(-1);
		}

		int nnvar, nngoals;
		char rfile[200];
		fscanf(fp, "%d%d%s",  &nnvar, &nngoals, rfile);
		if( strcmp(rfile, inputFileName)  )
		{
			fprintf(stderr,"Input file name %s is different from name in recover.hex %s\n", inputFileName, rfile);
			myexit(-1);
		}
		if(nnvar != nvar || nngoals != ngoals)
		{
			fprintf(stderr," parateter mismatch nvar=%d nnvar=%d ngoals=%d nngoals=%d\n", nvar, nnvar, ngoals, nngoals);
			myexit(-1);
		}
		fprintf(flog, "recovered nrun=%d, nvar=%d, ngoals=%d\n", nrun, nnvar, nngoals);

		char seperator[10];
		double val;
		nrun=0;
		while(1)
		{
			int it[2];
			double * tt=(double*)it;
			for(int j=0; j<nvar; j++)
			{
				if(fscanf(fp, "%x%x", it, it+1) != 2) goto ende;
				para[nrun][j]= *tt;
			}
			for(int j=0; j<ngoals; j++)
			{
				if(fscanf(fp, "%x%x", it, it+1) != 2) goto ende;
				subresu[nrun][j]= *tt;
			}
			if(fscanf(fp, "%x%x%s%le", it, it+1, seperator, &val) != 4) goto ende;
			if(seperator[0] != '=') goto ende;
			resu[nrun] = *tt;
			char * nobest= strdup("   ");
			char * isbest= strdup("###");
			char * showbest=nobest;
			if( resu[nrun] < bestgoal) 
			{
				showbest=isbest;
				bestgoal= resu[nrun];
			}
			fprintf(flog, "%sbest *old* %2d %13.5e this %13.5e <- ", showbest, -1, bestgoal, resu[nrun]);
			for(int j=0; j<nvar; j++)  fprintf(flog, "%13.5e ", para[nrun][j]);
			fprintf(flog, ": ");
			for(int j=0; j<ngoals; j++)  fprintf(flog, "%13.5e ", subresu[nrun][j]);
			fprintf(flog, "\n");
			nrun++;
		}
	ende:   fclose(fp);
		fprintf(flog, "recovered %d calculations\n", nrun);
		fflush(flog);
	}
// -----------------------------------------
    	t=1;	// one quality function
//     	strcpy(name,"mainp");


	if(nvar == 0 && !oneTime)
	{
		fprintf(stderr, "No variables set, exiting\n");
		myexit(-1);
	}
	if(ngoals == 0 && !oneTime)
	{
		fprintf(stderr, "No goals set, exiting\n");
		myexit(-1);
	}
	xOptimal.setSize(nvar); 	// sets the number of variables
	xStart.setSize(nvar);


	valueOptimal=0.0;
	
	for(int i=0; i<nvar; i++)  xStart[i]=0.;
	fprintf(flog, "end constructor\n");
}


void M2condor::done()
{
	fprintf(flog, "+++best ---------------done---------------");
}






double M2condor::eval()
{
	for(int i=0; i<nvar; i++)  xStart[i]=0;
	int error;
	double r= eval(xStart, &error, 0);
	if(error) fprintf(stderr, "eval returned error %d\n", error);
	return r;
}



double M2condor::eval(Vector X, int *nerror, int resource)
{
	if(nerror) *nerror=0;
    	double *p=X;
    	double *start=xStart;
	char resname[4];
	sprintf(resname,"%03d",resource);







	fprintf(flog,   "rtxMutex %d asks 1\n", resource); fflush(flog);
        pthread_mutex_lock( &rtxMutex );
	fprintf(flog,   "rtxMutex %d got it 1\n", resource); fflush(flog);
	fprintf(flog, "run %d,%2d start =", nrun,resource);
	for(int i=0; i<nvar; i++)  fprintf(flog, " %20.10e ", start[i]);
	fprintf(flog, "\n");
	fprintf(flog, "run %d,%2d parms =", nrun,resource);
	for(int i=0; i<nvar; i++)  fprintf(flog, " %20.10e ", p[i]);
	fprintf(flog, "\n");
	double pnow[nvar];
	for(int i=0; i<nvar; i++) pnow[i]=p[i]*var[i].mult;
	fprintf(flog, "run %d,%2d trunc =", nrun,resource);
	for(int i=0; i<nvar; i++)  fprintf(flog, " %20.10e ", pnow[i]);
	fprintf(flog, "\n");
	fflush(flog);
        pthread_mutex_unlock( &rtxMutex );
	fprintf(flog,   "rtxMutex %d release 1\n", resource); fflush(flog);
// 	if(boundErr)
// 	{
// 		if(nerror) *nerror=1;
// 		return 2.e33;
// 	}



//	check if we have done that before
	for(int i=0; i<nrun; i++)
	{
		for(int j=0; j<nvar; j++)
			if(fabs(p[j] -para[i][j]) > 1.e-8) goto nextj;
		fprintf(flog, "done that in step %d, %20.10e\n", i, resu[i]);
		return resu[i];

nextj:	;
	}


	// make temp file names. We have the original input "inputFileName" in the main directory. inputFileName is a class member
	// we make a copy with substitutions "tempInputFileName" in the resource directory
	// This file is named "gpt_dir1_dir2_resource.gpt".
	// dir1 and dir2 are the last two levels of directories
	// this file name is argument to parmela and
	// it is visible in the ps -ef command, so it helps that it includes the dir name
	// we need this file name with and without the calcdirs/resource part
	char tempInputFileName[300];
        char cwd1[300];
        strcpy(cwd1, cwd);
	int le=strlen(cwd1);
	char *p1=cwd1, *p2=cwd1;
	for(char *p3=cwd1; *p3; p3++) if(*p3 == '/') { p1=p2; p2=p3;  *p3='_';}
	sprintf(tempInputFileName,"gpt_%s_%03d_%s", p1+1, resource, inputFileName);
	// we need this file name with and without the calcdirs/resource part
	char dirTempInputFileName[300];
	sprintf(dirTempInputFileName,"calcdirs/%03d/%s", resource, tempInputFileName);

// 	printf("%2d  inputFileName		%s\n", resource, inputFileName);
// 	printf("%2d  tempInputFileName		%s\n", resource, tempInputFileName);
// 	printf("%2d  dirTempInputFileName	%s\n", resource, dirTempInputFileName);
	








        FILE * fi = fopen(inputFileName, "r");
	if(!fi)
	{
		fprintf(stderr, "cant open %s\n", inputFileName);
		myexit(-1);
	}
	

        FILE * fo = fopen(dirTempInputFileName, "w");
	if(!fo)
	{
		fprintf(stderr, "cant open %s\n", dirTempInputFileName);
		myexit(-1);
	}

	char line[500];
	char * pl;
	char what[100];
	int sofar;
	int nline=0;
	int ivar=0;
	while( (pl=fgets(line,499,fi)))
	{
		if(nline == var[ivar].nline)
		{

			fprintf(flog, "old: %s", line);
			double value = var[ivar].initVal + p[ivar]*var[ivar].mult;
			if( value > var[ivar].upper || value < var[ivar].lower)
			{
				if(nerror) *nerror=1;
				return 2.e33;
			}
			fprintf(fo, "%s = %20.10e  %s", var[ivar].name, value, var[ivar].restOfLine);
			fprintf(flog, "new: %s = %20.10e  %s", var[ivar].name, value, var[ivar].restOfLine);
			ivar++;
		}
		else
		{
			fprintf(fo, "%s", line);
		}

		nline++;
	}
	fclose(fi);
	fclose(fo);


	char command[1000];
#if 1
	fprintf(flog, "dir is %s, input is %s\n", cwd, inputFileName);



	for(int nsct =0; nsct < nsections; nsct++)
	{
		sprintf(command,"cd calcdirs/000 ;  gpt -j 1 -o outgpt.out%02d %s section=%d",   sections[nsct], tempInputFileName, sections[nsct] );
		memcpy(command+12, resname, 3);
		fprintf(flog, "execute %s\n", command);
		int ret = system(command);
		if(ret)
		{
			fprintf(stderr, "gpt call section %d did not exit normally\n", sections[nsct]);
			myexit(-1);
		}
	}


	fprintf(flog, "done gpt\n");
#else
	myexit(0);
	sprintf(command,"vi  %s", dirTempInputFileName); // here comments are stripped
	int ret = system(command);
	fprintf(flog, "done pm\n");
	return 0;
#endif




	ScreenList * sc  = new ScreenList;
	GDF * gdf[nsections];
	for(int nsct =0; nsct < nsections; nsct++)
	{
		char out2file[150];
		sprintf(out2file, "calcdirs/%s/outgpt.out%02d", resname, sections[nsct]);
		gdf[nsct] = new GDF(out2file); 
		for(ScreenList::iterator j=gdf[nsct]->screenList->begin(); j != gdf[nsct]->screenList->end(); j++)
		{
			sc->push_back(*j);
		}
	}

	// the very first time: get the number of particles at the start of the beam line.
	double goal=0.;
	double fit_this[maxFitGoals];		// result of the current evaluation
	int fit_this_found[maxFitGoals]; 
	for(int i = 0; i<ngoals; i++) fit_this_found[i] =0;
	for(int i = 0; i<ngoals; i++)
	{
		double where = gdf[0]->getDouble(fit_where[i]);
		for(ScreenList::iterator j=sc->begin(); j != sc->end(); j++)
		{
			Screen * s = *j;
			double fa=fabs(s->parameter - where);
			if( !int( fa*1000. ) )		// match within a mm
			{
				Particles p(*j);
				p.sortz(fit_cut[i]);
				// void Particles::twiss(int print, double ** pt, int np)
// 				fit_cut_ratio[i] =1.;
				p.twiss(0, p.scord, int(  p.npart * (1.-fit_cut_ratio[i]  )  )  );
				fit_this[i] = p.getValue(fit_what[i]);
				fprintf(flog, "value %s %e\n", fit_what[i], fit_this[i]);


				double diff = fabs((fit_val[i]-fit_this[i]));
				double gx = square(diff/fit_good[i]);
				const double a = 0.05;
				double gd = 1. - a/(a+square(gx));
				double dd = diff*fit_weight[i];
				double gg = square(diff*fit_weight[i])*gd;
// 				double gg = fabs(diff*fit_weight[i])*gd;
				goal += gg;

				fprintf(flog, "%s  figure dif %e gx %e gd %e dd %e gg %e \n", fit_what[i],  diff, gx, gd, dd, gg);
				fprintf(flog, "value %s at %s(%f) = %20.10e\n", fit_what[i], fit_where[i], where, fit_this[i]);
				fit_this_found[i] =1;
				break;

			}
		}
	}
	for(int i = 0; i<ngoals; i++)
	{
		if(!fit_this_found[i])
		{
			double where = gdf[0]->getDouble(fit_where[i]);
			fprintf(flog, "goal at %s = %f not found!\n", fit_where[i], where);
			myexit(-1);
		}
	}
	delete sc;
	for(int nsct =0; nsct < nsections; nsct++)
	{
		delete gdf[nsct];
	}





	fprintf(flog,   "rtxMutex %d asks 2\n", resource); fflush(flog);
        pthread_mutex_lock( &rtxMutex );
	fprintf(flog,   "rtxMutex %d got it 2\n", resource); fflush(flog);

	time_t tim = time(0);
	char * timc =ctime(&tim);
	timc[16]=0;
// 	char * nobest= strdup("   ");
// 	char * isbest= strdup("###");
	const char  nobest[]= {' ',' ',' ',0};
	const char  isbest[]= {'#','#','#',0};
	const char * showbest=nobest;
	if(goal < bestgoal)
	{
		showbest=isbest;
		bestgoal=goal;
		if(!oneTime)
		{
			sprintf(command,"mv  -f %s gpt_best.gpt", dirTempInputFileName);
// 			fprintf(flog, "execute %s\n", command);
			system(command);
		}
		sprintf(command,"cp  -f calcdirs/%03d/outgpt.out* .", resource);
		fprintf(flog, "execute %s\n", command);
		system(command);
	}
	fprintf(flog, "%sbest %5s %2d %13.5e this %13.5e <- ", showbest, timc+11, resource, bestgoal, goal);
	for(int i=0; i<nvar; i++)  fprintf(flog, "%13.5e ", p[i]);
	fprintf(flog, ": ");
	for(int i=0; i<ngoals; i++)  fprintf(flog, "%13.5e ", fit_this[i]);
	fprintf(flog, "\n");
	fflush(flog);

	for(int j=0; j<nvar; j++) para[nrun][j] = p[j];
	for(int j=0; j<ngoals; j++) subresu[nrun][j] = fit_this[j];
	resu[nrun++] = goal;

	// update recover  -- needs work
	FILE * fpr = fopen("recover.hex", "w");
	if(!fpr)
	{
		fprintf(stderr,"cant open %s for write\n", "recover.hex");
		myexit(-1);
	}

	fprintf(fpr, "%d %d %s\n",  nvar, ngoals, inputFileName);
	for(int i=0; i<nrun; i++)
	{
		for(int j=0; j<nvar; j++) 
		{
			int * it=(int*) (para[i]+j);
			fprintf(fpr,"%8x %8x ", it[0], it[1]);
		}
		for(int j=0; j<ngoals; j++)
		{
			int * it=(int*) (subresu[i]+j);
			fprintf(fpr,"%8x %8x ", it[0], it[1]);
		}
		int * it=(int*) (resu+i);
		fprintf(fpr,"%8x %8x = %e\n", it[0], it[1],  resu[i]);
	}
	fclose(fpr);


	FILE * fs = fopen("plot", "r");
	if(fs)
	{
		fclose(fs);
		system("rm plot");
		system(main5);
		fprintf(stderr, "****made plot*****\n");
	}

        pthread_mutex_unlock( &rtxMutex );
	fprintf(flog,   "rtxMutex %d release 2\n", resource); fflush(flog);
	
	
	while( (fs = fopen("pause", "r")))
	{
		fclose(fs);
		fprintf(stderr, "pause file found\n");
		sleep(30);
	}
        fs = fopen("stop", "r");
	if(fs)
	{
		fclose(fs);
		system("rm stop");
		fprintf(flog, "+++best             stop file found\n");
		myexit(0);
	}

	// turn of negative switches after the first run
	for(int i=0; i<nswitches; i++)
		if( switchFlag[i] < 0)  switchFlag[i] = 0;
   	updateCounter(goal,X);
	return goal;
}							









int main(int argc, char ** argv, char ** env)
{
	int recover =0;
	int oneTime =0;
	int display =0;
        while(argc > 2)
        {
		if( !strcmp(argv[1], "-r") )
		{
			recover=1;
			argv++; argc--;
		}
		if( !strcmp(argv[1], "-d") )
		{
			display=1;
			argv++; argc--;
		}
		if( !strcmp(argv[1], "-o") )
		{
			oneTime=0;
			argv++; argc--;
		}
		if( !strcmp(argv[1], "-p") )
		{
			oneTime=1;
			argv++; argc--;
		}
	}

        if(argc != 2)
        {
                fprintf(stderr, "usage: %s <input file>\n", argv[0]);
                myexit(-1);
        }


	pid_t myPid = getpid();
	FILE * fpr=fopen(thislock, "r");
	if(fpr)
	{
		fprintf(stderr, "file %s exists, %s may already be running\n", thislock, argv[0]);
		exit(-1);
	}
	fpr=fopen(thislock, "w");
	if(!fpr)
	{
		perror(thislock);
		exit(-1);
	}
	fprintf(fpr, "%d\n", myPid);
	fclose(fpr);

	unlink("stop");
	unlink("pause");
	unlink("plot");

	char inputFileName[300];
	char fileroot[300];		// input file root
        strcpy(inputFileName,argv[1]);


	int ls=strlen(inputFileName);
	if(strcmp(inputFileName+ls-4,".gpt") )
	{
		fprintf(stderr, "input file name must end in .gpt, is %s\n", inputFileName+ls-4);
		myexit(-1);
	}

	strcpy(fileroot, inputFileName);
	fileroot[ls-4]=0;




	char * home = getenv("HOME");
	if(!home)
	{
		fprintf(stderr, "environment variable HOME not set\n");
		myexit(-1);
	}
	
	
	char xfile[300];
	char rfile[300];
	
	
	if(display)
	{
		for(display=60; display < 65; display++)
		{
	
			sprintf(rfile,"/tmp/.RTX%d-lock", display);
			sprintf(xfile,"/tmp/.X%d-lock", display);
			FILE * fpr=fopen(rfile, "r");
			if(fpr)
			{
				pid_t rpid = -1;
				pid_t xpid = -1;
				fscanf(fpr, "%d%d", &rpid, &xpid);
				fclose(fpr);
				if(rpid > 0)
				{
					int err = kill(rpid, 0);
					if(err ==  0) continue;		// display in use
				}
				unlink(rfile);
			}
			goto isset;
		}
		fprintf(stderr, "all displays are in use\n");
		myexit(-1);
				

isset:		fprintf(stderr, "display is %d\n", display);
		FILE * fpx=fopen(xfile, "r");
		if(fpx) 
		{
			pid_t xpid = -1;
			fscanf(fpx, "%d", &xpid);
			fclose(fpx);
			if(xpid > 0)
			{
				int err = kill(xpid, SIGTERM);
				if(err ==  0) unlink(xfile);
				else sleep(4);
			}
		}
	
	
	
		char * arg[6];
		arg[0] = strdup("/usr/bin/Xvfb");
		arg[1] = strdup("-sp");
		arg[2] = new char[200]; strcpy(arg[2], home); strcat(arg[2], "/SecurityPolicy");
		arg[3] = new char[200]; sprintf(arg[3], ":%d", display);
		arg[4]=0;
		setenv("DISPLAY", arg[3],1);
	
    		pidXvfb= vfork();
    		if(pidXvfb == 0)
    		{
         		execve(arg[0], arg, env);
            		_exit(2);
     		}
		

	             
		FILE * fpr=fopen(rfile, "w");
		if(fpr)
		{
			fprintf(fpr, "%d %d\n", myPid, pidXvfb);
			fclose(fpr);
		}
	}




	setenv("OMP_WAIT_POLICY", "PASSIVE", 1);




    	double rhoStart=1e-0, rhoEnd=1e-5;
    	int niter=100000;


    	M2condor *mof = new M2condor(inputFileName, recover, oneTime);
    	ObjectiveFunction *of = mof;

	// set number of threads. argument = 0 means look at /proc/cpuinfo 
	// negative means /proc/cpuinfon but reserve -n threads for other things
	int maxthreads = of->setMaxThreads(0);
	int nvar=mof->nvar;
	int neededThreads = nvar*(nvar-1)/2;
	if(neededThreads < nvar) neededThreads= nvar;
	fprintf(mof->flog, "available threads %d\n", maxthreads);
	fprintf(mof->flog, "parameters %d needed threads %d\n",  mof->nvar, neededThreads);
	if(maxthreads > neededThreads) maxthreads = of->setMaxThreads(neededThreads);
	if(oneTime)
	{
		maxthreads =1;
		maxthreads = of->setMaxThreads(maxthreads);
	}
	fprintf(mof->flog, "will use %d threads\n", maxthreads);

	// create resources
	int iret=  system("rm -rf calcdirs");
	if(iret)
	{
		fprintf(stderr, "can't delete calcdirs directory");
		perror("calcdirs");
		myexit(-1);
	}
	fprintf(mof->flog, "rm -rf calcdirs");
	iret=  mkdir("calcdirs", 0755);
	if(iret)
	{
		fprintf(stderr, "can't make calcdirs directory");
		perror("calcdirs");
		myexit(-1);
	}
// 	system("mkdir calcdirs");
	for(int i=0; i<maxthreads; i++)
	{
		char dest[200];
		sprintf(dest, "calcdirs/%03d", i);
		fprintf(mof->flog, "mkdir %s\n", dest);
		int iret=  mkdir(dest, 0755);
		if(iret)
		{
			fprintf(stderr, "can't make %s directory", dest);
			perror(dest);
			myexit(-1);
		}

// 		//copy initial distribution file
		int sec = mof->sections[0]; // copy only the first file
		if( sec > 0)
		{
			char cmd[200];
			sprintf(cmd, "cp -f outgpt.out%02d  calcdirs/%03d", sec-1, i);
			fprintf(mof->flog, "%s\n", cmd);
			iret = system(cmd);
			if(iret)
			{
				perror(cmd);
				myexit(-1);
			}
		}

	}       



	if(oneTime)
	{
		mof->eval();
	}
	else
	{
		CONDOR(rhoStart, rhoEnd, niter, of);
	}


	mof->done();
    	delete of;
	unlink(thislock);

}
