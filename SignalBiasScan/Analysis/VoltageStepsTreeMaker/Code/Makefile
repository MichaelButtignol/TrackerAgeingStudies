GCC = g++
LD = g++
CXXFLAGS = $(shell root-config --cflags) 
LIBS = $(shell root-config --libs)
SRCS = Dict.C VoltageStepsTreeMaker.C Run_TreeMaker_compiled.C
OBJS = Dict.o VoltageStepsTreeMaker.o Run_TreeMaker_compiled.o
PROGRAM = Run_TreeMaker_compiled
HEADERS = ../../../interface/TreeEvent.h
DEFINES = -DNO_ORCA_CLASSES


all:	Dict compile link

compile:	$(OBJS)

link:	$(OBJS)
	@echo "Linking object files ..."
	$(LD) $(OBJS) $(LIBS) -o $(PROGRAM)

Dict: $(HEADERS) LinkDef.h
	@echo "Generating dictionary Dict.C ..."
	rootcint -f Dict.C -c $(DEFINES) $(HEADERS) LinkDef.h
	@echo "Compiling source files ..."

clean:; @rm -f $(OBJS) Dict.h Dict.C $(PROGRAM) *~ */*~

