name := basicTest
name := bertiniDataTest
name := cascade
name := cascadeTest
name := data
name := evaporationTest
name := hetcCascadeTest
name := testInterface
name := toyModelTest
name := benchmarks

G4TARGET := $(name)
#G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/binmake.gmk


b1:  # benchmark 1
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/benchmarks >> b1.out; date; cat b1.out 