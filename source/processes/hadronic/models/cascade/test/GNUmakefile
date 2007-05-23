# $Id: GNUmakefile,v 1.51 2007-05-23 12:13:54 miheikki Exp $
#G4TARGET := testIncl
G4TARGET := cascade
G4TARGET := benchmarks

#G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/binmake.gmk

b1:  # benchmark 1
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/benchmarks 1 >> b1.out; date; cat b1.out 

c3:  # compare with kalbach05 Fig.11 c)
	rm -f c3.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 62 10000 1 62 197 79 > c3.out; date; head c3.out; cp c3.out data/data.out 

c2:  # compare with heikkinen07::: 1(b) Bismut isotope 100% neutron production
	rm -f c2.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 90 10000 1 90 209 83 > c2.out; date; head c2.out; cp c2.out data/data.out 

c1:  # test run for cascade  parameters: runId #events particle energy [MeV] A Z
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1000 1000 1 1000 27 13 > c1.out; date; head c1.out; cp c1.out data/data.out 

