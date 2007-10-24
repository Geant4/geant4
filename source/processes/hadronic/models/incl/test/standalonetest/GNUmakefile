# $Id: GNUmakefile,v 1.2 2007-10-24 15:45:47 miheikki Exp $
#G4TARGET := cascade
#G4TARGET := benchmarks
G4TARGET := inclStandAlone
#G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/binmake.gmk

# Use ROOT for data analysis
CPPFLAGS += $(shell root-config --cflags)
LDFLAGS  += $(shell root-config --libs)

CPPFLAGS += -DUSEROOT
TESTPROGRAM = $(G4WORKDIR)/bin/$(G4SYSTEM)/inclStandAlone

incl:
	rm -f incl.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/testIncl 1  > incl.out; date; head incl.out; 

t1:
	$(TESTPROGRAM) 27 13 1 256 100000 ../data/testrunAl27Proton256MeV_noabla.root

t2:
	$(TESTPROGRAM) 27 13 1 800 100000 ../data/testrunAl27Proton800MeV_noabla.root

t3:
	$(TESTPROGRAM) 56 26 1 800 100000 ../data/testrunFe56Proton800MeV_noabla.root

demo:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclStandAlone 208 82 9 1000 10000 Pb208Alpha1GeV.root
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclStandAlone 208 82 1 1000 10000 Pb208Proton1GeV.root
	root -l compareDemo.C

testrun:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclStandAlone 208 82 1 1000 10000 test.out
#	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclStandAlone 27 13 1 1000 1 test.out > console.out

vg:
	valgrind $(G4WORKDIR)/bin/$(G4SYSTEM)/inclStandAlone 208 82 1 1000 10 valgrindtest.out 2> vg.output
#	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclStandAlone 27 13 1 1000 1 test.out > console.out
