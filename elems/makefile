# UNIX Makefile for ELEMS

include ../makeinc.gcc

TARGETS = gpt$(EXEEXT)
all: $(TARGETS)

include makefile.ele
ELEMOBJSO = $(ELEMOBJS:.obj=.o)

CFLAGS += -g -I../kernel -I.

elemlist:
	echo "# List your custom GPT-elements in two columns: name filename (where the filename must be listed WITHOUT .c extension)" > elemlist

gps.c: elemlist
	../bin/makeelems ../kernel/elemkrnl elemlist
	@echo Please rerun make.

gpt: gps.o $(ELEMOBJSO) ../kernel/elems.a ../kernel/gps.a ../kernel/gdf.a ../kernel/utils.a
	$(CC) -g $(LDFLAGS) $(TARGET_ARCH) -o $@ gps.o $(ELEMOBJSO) ../kernel/elems.a ../kernel/gps.a ../kernel/gdf.a ../kernel/utils.a
	$(CP) gpt ../bin

clean:
	$(RM) gpt $(ELEMOBJSO) gps.c
