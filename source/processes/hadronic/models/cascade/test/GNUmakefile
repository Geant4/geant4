all: 	# :::              	GEANT4 CASCADE TESTING 
	grep ":::" GNUmakefile
# 	  ::: GNUmakefile collects tools for Geant4 hadronic cascade testing.  

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
	-lG4globman \
	-lG4hadronic_bert_cascade \
	-lG4hepnumerics \
	-lCLHEP

LOADLIBS := $(G4BIN)/$(G4SYSTEM)/hbfinit.o -l$(CLHEP_LIB)
LOADLIBS +=    `cernlib packlib` $(FCLIBS)

# 	  ::: Parameter list for cascade : nEvents, particleType, bulletE [MeV], N, and Z.

#t1:	# ::: problem 1  (Al)  evaporation missing, SOLVED
#	rm -r cascade.out; cascade 10000 1 256 27 13 > cascade.out; gmake a; 

p2:	# ::: p(roblem) 2 (Al) strong energy tail ::: general feature, study particle history, MAKE NOTE
	rm -r cascade.out; cascade 10000 1 160 27 13 > cascade.out; gmake a; 

p3:	# ::: problem 3 (Zr) weak medium energy ::: unknown reason POSTPONED
	rm -r cascade.out; cascade 10000 1 256 91 40 > cascade.out; gmake a; 

p4:	# ::: problem 4 (Be) too few pi+'s ::: general feature, POSTPONED
	rm -r cascade.out; cascade 10000 1 585 9 4 > cascade.out; gmake a; 

p5:	# ::: problem 5 (Al) too many n's in medium energy ::: hunt using history data, POSTPONED
	rm -r cascade.out; cascade 10000 1 585 27 13 > cascade.out; gmake a;

#t6:	# ::: problem 6 () too few neutrons, SOLVED
#	rm -r cascade.out; cascade 10000 1 800 207 82 > cascade.out; gmake a;

p7:	# ::: problem 7 (H)  missing pions ::: CRITICAL, TO BE DONE
	rm -r cascade.out; cascade 1000 1 585 1 1 > cascade.out; make a; 

H1:	# ::: (H) A = 1 protiun
	rm -r cascade.out cascade.eps; cascade 25 1 585 1 1 > cascade.out; cat cascade.out
#	rm -r cascade.out cascade.eps; cascade 1000 1 585 1 1 > cascade.out; make a;

H2:	# ::: (H) A = 2 deuterium
	rm -r cascade.out cascade.eps; cascade 25 1 585 2 1 > cascade.out; cat cascade.out
#	rm -r cascade.out cascade.eps; cascade 1000 1 585 2 1 > cascade.out; make a;

H3:	# ::: (H) A = 3 tritium
	rm -r cascade.out cascade.eps; cascade 1000 1 585 3 1 > cascade.out; make a;

p8:	# ::: problem 8 (Zr) delta missing ::: general feature in codes, POSTPONED studies
	rm -r cascade.out cascade.eps; cascade 1000 2 800 91 40 > cascade.out; make a;

leak:	# ::: find memory leaks
	valgrind cascade 1

time: 	# :: timing for target t
	time -v gmake t

a: 	# ::: a(nalysis) of data in  cascade.dat with batch run of Root
	rm cascade.eps; root -b -q cascade.C; gv cascade.eps; 

b: 	# ::: b(enchmark) testing
	rm -r benchmarks benchmark.out; make benchmarks; benchmarks > benchmarks.out; diff benchmarks.out benchmarks.data

t:	# ::: make test cascade
	rm cascade; make cascade

eps: 	# ::: generate eps-files for documentation
	cp *.eps $(RELD)/pictures/

clean:  # ::: clean
	rm -f \#*.*\# \#*\# \#* .*~ .#* *.*~ *~ *.o; rm -f core cascade ; rm benchmarks; rm *.eps *.root



