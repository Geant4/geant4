# // CLHEP
CLHEPD    := /usr/local/CLHEP
CLHEPINC  := -I$(CLHEPD)/include
CLHEPLIB  := -L$(CLHEPD)/lib -lCLHEP

# // HepMC
HEPMCD     := /opt/heplib/HepMC/2.06.03
HEPMCINC   := -I$(HEPMCD)/include
HEPMCLIB   := -L$(HEPMCD)/lib -lHepMC -lHepMCfio

# // CERN libraries
CLIB      := $(shell cernlib -v pro pythia6205 pdflib804 mathlib) -lg2c

CXXFLAGS := $(CLHEPINC) $(HEPMCINC)
LOPT := -Wl,-rpath $(HEPMCD) $(HEPMCLIB) $(CLHEP) $(CLIB) -lg2c

%.o: %.cxx
	$(CXX) -c -o $@ $(CXXFLAGS) $*.cxx

example_MyPythia : example_MyPythia.o
	$(CXX) -o $@ $< $(LOPT)

example_MyPythia.o : example_MyPythia.cxx

