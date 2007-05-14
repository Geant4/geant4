name := basicTest
name := benchmarks
name := bertiniDataTest
name := cascade
name := cascadeTest
name := data
name := evaporationTest
name := hetcCascadeTest
name := testInterface
name := toyModelTest

G4TARGET := $(name)
#G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/binmake.gmk


b1:  # benchmark 1
	$(G4INSTALL)/bin/$(G4SYSTEM)/benchmarks >> b1.out 