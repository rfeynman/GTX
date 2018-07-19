
allprog=   readGDF gplot glongPlot timeExtract muppel plotTraj
#allexe := $(allprog:%=%.exe)           # for windows
allexe=$(allprog)


all:	$(allprog)

CC = clang
CCOPT=-g
LDOPT=-lm -lstdc++


readGDF:        readGDF.cxx
	                $(CC) -o readGDF readGDF.cxx $(LDOPT)





muppel:	muppel.cxx gdf.o $(HOME)/jkLib $(HOME)/jkLib/libjkLib.a
	$(CC) -o muppel muppel.cxx  gdf.o -I$(HOME)/condorCommonN $(HOME)/condorCommonN/libcondor.a -I$(HOME)/jkLib $(HOME)/jkLib/libjkLib.a -lpthread $(CCOPT) $(LDOPT)



gplot:	gplot.cxx gdf.o  $(HOME)/jkLib $(HOME)/jkLib/libjkLib.a 
	$(CC) -o gplot gplot.cxx gdf.o  -I$(HOME)/jkLib $(HOME)/jkLib/libjkLib.a  $(CCOPT) $(LDOPT)

plotTraj:	plotTraj.cxx gdf.o  $(HOME)/jkLib $(HOME)/jkLib/libjkLib.a 
	$(CC) -o plotTraj plotTraj.cxx gdf.o  -I$(HOME)/jkLib $(HOME)/jkLib/libjkLib.a  $(CCOPT) $(LDOPT)

glongPlot:	glongPlot.cxx gdf.o  $(HOME)/jkLib $(HOME)/jkLib/libjkLib.a 
	$(CC) -o glongPlot glongPlot.cxx gdf.o  -I$(HOME)/jkLib $(HOME)/jkLib/libjkLib.a  $(CCOPT) $(LDOPT)

timeExtract:	timeExtract.cxx gdf.o  $(HOME)/jkLib $(HOME)/jkLib/libjkLib.a 
	$(CC) -o timeExtract timeExtract.cxx gdf.o  -I$(HOME)/jkLib $(HOME)/jkLib/libjkLib.a  $(CCOPT) $(LDOPT)








gdf.o:	gdf.cxx gdf.hxx
	$(CC) gdf.cxx -c $(CCOPT)


clean:
	rm -f *.o $(allexe)
