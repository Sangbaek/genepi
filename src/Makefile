CCSRCS     = $(wildcard *.cc)
EXE        = genepi.exe
FSRCS      = $(wildcard *.f)

CXX        = g++
#CXX        = g++-3.4

# CCFLAGS    = -g -Wall -fPIC -Wno-deprecated
CCFLAGS    += `$(ROOTSYS)/bin/root-config --cflags`
CCFLAGS    += -I../include -I. -I$(ROOTSYS)/include

ROOTLIB    = `$(ROOTSYS)/bin/root-config --glibs`
CERNLIBS =  -L/Users/sangbaek/CLAS12/cernlib  -lmathlib  -lpacklib -lkernlib
JLIB       = ../jetset_lib/libjetset.a
#MYLIBS     = -Wl,-rpath, -lm -lg2c -L../jetset_lib -ljetset
MYLIBS     = -Wl,-rpath, -lm -L/usr/local/Cellar/gcc/8.1.0/lib/gcc/8 -lgfortran -L../jetset_lib -ljetset
LIBS       = $(ROOTLIB) $(CERNLIB) $(MYLIBS)

CCOBJS     = $(CCSRCS:%.cc=%.o)

all: $(JLIB) $(EXE)

$(JLIB):
	@echo compiling Jetset lib
	make -C ../jetset_lib

$(EXE): $(CCOBJS)
	$(CXX) $(CCFLAGS) -o $@ $(CCOBJS) $(FOBJS) $(LIBS)
	cp $(EXE) ../bin

.cc.o :
	@echo Compiling $@
	@$(CXX) $(CCFLAGS) -c $< -o $@

clean:
	rm -f $(EXE) $(CCOBJS) ../bin/$(EXE)

cleanall:
	rm -f $(EXE) $(CCOBJS) ../bin/$(EXE) ../jetset_lib/*.o $(JLIB)

