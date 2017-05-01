#============================================================================
#  a appeler avec le gnumake (commande gmake). Le prefixe GNU fait
# qu'il n'est utilisable qu'avec gnumake.
# Suivant les machines, root est compile avec des compilateurs differents.
# Ce makefile en tient compte et appelle le bon compilateur selon la machine.
#============================================================================
.SUFFIXES: .cxx .o

ARCH=$(shell uname -s)

CC = g++


.cxx.o:
	$(CC) $(CXXFLAGS) -c $*.cxx

CXXFLAGS      =  -O -g -Wall -fPIC
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)


CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)
RLIBS         = $(ROOTLIBS)  $(ROOTGLIBS)  

ALL=DA

all: $(ALL)


DEPS = include/DA.h include/WaveProcessor.h
OBJ = DA.o WaveProcessor.o typeConvert.o


WaveProcessor.o:	src/WaveProcessor.cxx include/WaveProcessor.h 
	$(CC) -c -o $@ $< $(CXXFLAGS)

typeConvert.o:	src/typeConvert.cxx 
	$(CC) -c -o $@ $< $(CXXFLAGS)

DA.o:	src/DA.cxx include/DA.h 
	$(CC) -c -o $@ $< $(CXXFLAGS)


DA: $(OBJ)
	$(CC) -o $@ $^ $(RLIBS)

clean:
	rm -f *.o *~ corre
	rm -f $(ALL)









