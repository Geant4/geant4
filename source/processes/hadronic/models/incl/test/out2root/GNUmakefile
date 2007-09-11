# $Id: GNUmakefile,v 1.1 2007-09-11 13:35:42 miheikki Exp $
#G4TARGET := cascade
#G4TARGET := benchmarks
G4TARGET := out2root
#G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/binmake.gmk

# Use ROOT for data analysis
CPPFLAGS += $(shell root-config --cflags)
LDFLAGS  += $(shell root-config --libs)

CPPFLAGS += -DUSEROOT

b1:  # benchmark 1
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/benchmarks 1 >> b1.out; date; cat b1.out 

b2:  # test run for cascade  parameters: runId #events particle energy [MeV] A Z
#	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1 1000 1 1000 27 13 > ../data/b2.out; date; head ../data/b2.out; wc ../data/b2.out 
#	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1 1000 1 3000 27 13 > b2.out; date; head b2.out; 


incl:
	rm -f incl.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/testIncl 1  > incl.out; date; head incl.out; 

testrun:
#	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclStandAlone 208 82 1 1000 1 test.out > console.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclStandAlone 27 13 1 1000 1 test.out > console.out
