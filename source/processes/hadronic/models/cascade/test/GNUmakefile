all: 	# :::              	GEANT4 LOW ENERGY BERTINI CASCADE TESTING 
#       :::                             by aatos.heikkinen@cern.ch
	grep ":::" GNUmakefile
#	::: GNUmakefile collects tools for Geant4 hadronic cascade testing for E < 1 GeV  
#       ::: Do 'make global' first for Geant4 installation
#       ::: Example of usage: make compilation, make run1 
ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -Wall

CPPFLAGS  += \
	-I$(G4BASE)/processes/hadronic/models/generator/include \
	-I$(G4BASE)/global/management/include \
	-I$(G4BASE)/global/HEPRandom/include \
	-I$(G4BASE)/global/HEPGeometry/include \
	-I$(G4BASE)/global/HEPNumerics/include \
	-I$(G4BASE)/intercoms/include \
	-I$(G4BASE)/materials/include \
	-I$(G4BASE)/particles/management/include \
	-I$(G4BASE)/particles/bosons/include \
	-I$(G4BASE)/particles/leptons/include \
	-I$(G4BASE)/particles/hadrons/mesons/include \
	-I$(G4BASE)/particles/hadrons/barions/include \
	-I$(G4BASE)/particles/hadrons/ions/include \
	-I$(G4BASE)/particles/shortlived/include \
	-I$(G4BASE)/processes/general/include \
	-I$(G4BASE)/processes/management/include \
	-I$(G4BASE)/processes/parameterisation/include \
	-I$(G4BASE)/processes/hadronic/cross_sections/include \
	-I$(G4BASE)/processes/hadronic/management/include \
	-I$(G4BASE)/processes/hadronic/models/cascade/cascade/include \
	-I$(G4BASE)/processes/hadronic/processes/include \
	-I$(G4BASE)/processes/hadronic/util/include \
	-I$(G4BASE)/track/include \
	-I$(G4BASE)/tracking/include 

LDLIBS   +=  \
	-L$(G4INSTALL)/lib/$(G4SYSTEM) \
	-lG4digits_hits \
	-lG4event \
	-lG4FR \
	-lG4geometry \
	-lG4global \
	-lG4graphics_reps \
	-lG4intercoms \
	-lG4interfaces \
	-lG4materials \
	-lG4modeling \
	-lG4OpenGL \
	-lG4parmodels \
	-lG4particles \
	-lG4persistency \
	-lG4processes \
	-lG4RayTracer \
	-lG4readout \
	-lG4run \
	-lG4tracking \
	-lG4track \
	-lG4Tree \
	-lG4visHepRep \
	-lG4vis_management \
	-lG4visXXX \
	-lG4VRML \
	-lG4global \
	-lCLHEP-g++.1.8.1.0 \

LOADLIBS := $(G4BIN)/$(G4SYSTEM)/hbfinit.o -l$(CLHEP_LIB)
LOADLIBS +=    `cernlib packlib` $(FCLIBS)

compilation:	# ::: make cascade executable
	rm cascade; make cascade

#	::: Input parameters for cascade : #events, particleType, bulletE [MeV], target A = Z + N, and target Z.
#       ::: particleType: nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10;

#	::: Output file format: enegy [GeV] #events, #protons, #neutrons, #pions+ #pions-, #pions0, #fotons
#	::: # of particles are divided by number of events

run1:	# ::: 1-10 MeV
	rm -r cascade.out;
	./cascade 100 1 1  58 28 > cascade.out;
	./cascade 100 1 2  58 28 >> cascade.out; 
	./cascade 100 1 3  58 28 >> cascade.out; 
	./cascade 100 1 4  58 28 >> cascade.out; 
	./cascade 100 1 5  58 28 >> cascade.out; 
	./cascade 100 1 6  58 28 >> cascade.out; 
	./cascade 100 1 7  58 28 >> cascade.out; 
	./cascade 100 1 8  58 28 >> cascade.out; 
	./cascade 100 1 9  58 28 >> cascade.out; 
	./cascade 100 1 10  58 28 >> cascade.out; 
	cat cascade.out

run10:	# ::: 10-100 MeV
	rm -r cascade.out;
	./cascade 100 1 10   58 28 > cascade.out;
	./cascade 100 1 20   58 28 >> cascade.out; 
	./cascade 100 1 30   58 28 >> cascade.out; 
	./cascade 100 1 40   58 28 >> cascade.out; 
	./cascade 100 1 50   58 28 >> cascade.out; 
	./cascade 100 1 60   58 28 >> cascade.out; 
	./cascade 100 1 70   58 28 >> cascade.out; 
	./cascade 100 1 80   58 28 >> cascade.out; 
	./cascade 100 1 90   58 28 >> cascade.out; 
	./cascade 100 1 100  58 28 >> cascade.out; 
	cat cascade.out


run100:	# ::: 100 - 1000 MeV
	rm -r cascade.out;
	./cascade 100 1 100  58 28 > cascade.out;
	./cascade 100 1 200  58 28 >> cascade.out; 
	./cascade 100 1 300  58 28 >> cascade.out; 
	./cascade 100 1 400  58 28 >> cascade.out; 
	./cascade 100 1 500  58 28 >> cascade.out; 
	./cascade 100 1 600  58 28 >> cascade.out; 
	./cascade 100 1 700  58 28 >> cascade.out; 
	./cascade 100 1 800  58 28 >> cascade.out; 
	./cascade 100 1 900  58 28 >> cascade.out; 
	./cascade 100 1 1000 58 28 >> cascade.out; 
	cat cascade.out




