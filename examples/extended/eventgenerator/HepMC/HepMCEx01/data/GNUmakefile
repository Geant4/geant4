# // CLHEP
CLHEPD    := /usr/local/CLHEP
CLHEPINC  := -I$(CLHEPD)/include
CLHEPLIB  := -L$(CLHEPD)/lib -lCLHEP

# // HepMC
HEPMCD     := /opt/heplib/HepMC/1.21
HEPMCINC   := -I$(HEPMCD)
HEPMCLIB   := -L$(HEPMCD) -lHepMC
#HEPMCLIB   := $(HEPMCD)/libHepMC.so

# // CERN libraries
CLIBD     := /cern/pro/lib
CLIB      := -L$(CLIBD) -lpythia6

CXXFLAGS := $(CLHEPINC) $(HEPMCINC)
LOPT := -Wl,-rpath $(HEPMCD) $(HEPMCLIB) $(CLHEP) $(CLIB) -lg2c

%.o: %.cxx
	$(CXX) -c -o $@ $(CXXFLAGS) $*.cxx

example_MyPythia : example_MyPythia.o
	$(CXX) -o $@ $< $(LOPT)

example_MyPythia.o : example_MyPythia.cxx

