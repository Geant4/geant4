# $Id: GNUmakefile,v 1.43 2007-05-20 11:21:56 miheikki Exp $
# $Log: not supported by cvs2svn $
name := cascade
#name := benchmarks

G4TARGET := $(name)
#G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/binmake.gmk


b1:  # benchmark 1
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/benchmarks 1 >> b1.out; date; cat b1.out 


b2:  # test run for cascade  parameters: runId #events particle energy [MeV] A Z
#	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1 1000 1 1000 27 13 > ../data/b2.out; date; head ../data/b2.out; wc ../data/b2.out 
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1 1000 1 3000 27 13 > b2.out; date; head b2.out; 



