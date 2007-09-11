# $Id: GNUmakefile,v 1.1 2007-09-11 13:34:09 miheikki Exp $
#G4TARGET := cascade
#G4TARGET := benchmarks
G4TARGET := testAbla
#G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/binmake.gmk

# ROOTFLAGS := $(shell root-config --cflags)
# ROOTLIB := $(shell root-config --libs) 

# CPPFLAGS += ROOTFLAGS
# LDFLAGS += ROOTLIB


b1:  # benchmark 1
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/benchmarks 1 >> b1.out; date; cat b1.out 

b2:  # test run for cascade  parameters: runId #events particle energy [MeV] A Z
#	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1 1000 1 1000 27 13 > ../data/b2.out; date; head ../data/b2.out; wc ../data/b2.out 
#	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1 1000 1 3000 27 13 > b2.out; date; head b2.out; 


incl:
	rm -f incl.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/testIncl 1  > incl.out; date; head incl.out; 

