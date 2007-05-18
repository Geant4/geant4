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
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1 100 1 100 27 13 > b2.out; date; cat b2.out 



