all: 	# :::              	GEANT4 CASCADE TESTING 
	grep ":::" GNUmakefile
# 	  ::: GNUmakefile collects tools for Geant4 hadronic cascade testing.  
#         ::: Tested in CERN lxplus Linux cluster.

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

# 	  ::: Parameter list for cascade : nEvents, particleType, bulletE [MeV], A = Z + N, and Z.

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

d:	# ::: study cascade model parameters
	rm data; make data; rm -f data.out; data > data.out; cat data.out ; # make da;

i:	# ::: test interface
	rm testInterface; make testInterface; testInterface

h:	# ::: compare h2 h1
	make h2; make h1

h1:	# ::: (H) A = 1 protiun 2 * 585 = 1170
#	make t; rm -f cascade.out; cascade  10000 1 585 1 1 > cascade.out; make a;
	make t; rm -f cascade.out; cascade  100 1 1000 1 1 > cascade.out; cat cascade.out
#	make t; rm -f cascade.out; cascade  10000 1 585 1 1 > cascade.out; make a;


h2:	# ::: (H) A = 2 deuterium 
#	rm -f cascade.out; cascade 10000 1 585 2 1 > cascade.out; make a;
	make t; rm -f cascade.out; cascade 1 1 585 2 1 > cascade.out; cat cascade.out

h3:	# ::: (H) A = 3 tritium
#	rm -f cascade.out; cascade 1000 1 585 3 1 > cascade.out; make a;
	rm -f cascade.out; cascade 3 1 585 3 1 > cascade.out; cat cascade.out

p8:	# ::: problem 8 (Zr) delta missing ::: general feature in codes, POSTPONED studies
#	rm -f cascade.out; cascade 1000 2 800 91 40 > cascade.out; make a;
	rm -f cascade.out; cascade 3 2 800 91 40 > cascade.out; cat cascade.out

p9:	# ::: problem 9 (Be) wrong spectrum, maybe too little energy goes to evaporation? (500 MeV ok)
	rm -r cascade.out; cascade 10000 1 130 9 4 > cascade.out; gmake a; 
#	rm -r cascade.out; cascade 10000 1 585 9 4 > cascade.out; gmake a; 

leak:	# ::: find memory leaks
	valgrind -v --leak-check=yes  --show-reachable=yes cascade 1

time: 	# :: timing for target t
	time -v gmake t

da:
	rm -f data.eps; root -b -q data.C; gv data.eps;

a: 	# ::: a(nalysis) of data in  cascade.dat with batch run of Root
	rm -f cascade.eps; root -b -q cascade.C; gv cascade.eps; 

b: 	# ::: b(enchmark) testing
	rm -r benchmarks benchmark.out; make benchmarks; benchmarks > benchmarks.out; diff benchmarks.out benchmarks.data

t:	# ::: make test cascade
	rm cascade; make cascade

dox:	# ::: gnerate HTML and latex manuals of source using doxygen 
	cd $(GI); p clean; cd $(GS); p clean; cd ~/public/html/system/;doxygen doxygen.config; cd $(RELD)/doxygen/latex/;make ps_2on1

eps: 	# ::: generate eps-files for documentation. Naming h1p585k10.eps 
#         ::: general naming skeme for plots: (subject in plots to be defined, sofar jus general collection)
#	  ::: h:::z:::n:::p:::e::::::c::::::.eps A with Z, N, particle type, E in MeV, collisions in (1000x) histype
#	  ::: Example: z001n002p001e000585c010h000.eps
#	  ::: most general plot listing all describing target A with Z = 1, N= 2 
#         ::: collided with bullet type 1 (proton) at 585 MeV 10 thousandtimes 
#
	cp *.eps $(RELD)/pictures/ 

clean:  # ::: clean
	rm -f \#*.*\# \#*\# \#* .*~ .#* *.*~ *~ *.o; rm -f core cascade ; rm benchmarks; rm *.eps *.root



