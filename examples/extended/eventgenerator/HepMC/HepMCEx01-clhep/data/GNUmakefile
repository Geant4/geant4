# // CLHEP
CLHEPD    := $(CLHEP_BASE_DIR)
CLHEPINC  := -I$(CLHEPD)/include
CLHEPLIB  := -L$(CLHEPD)/lib -lCLHEP

# // CERN libraries
CLIBD     := /cern/pro/lib
CLIB      := -L$(CLIBD) -lpythia6

CXXFLAGS := $(CLHEPINC)
LOPT := $(CLHEPLIB) $(CLIB) -lg2c

%.o: %.cc
	$(CXX) -c -o $@ $(CXXFLAGS) $*.cc

example_MyPythia : example_MyPythia.o
	$(CXX) -o $@ example_MyPythia.o $(LOPT)

example_MyPythia.o : example_MyPythia.cc

