
name := TestEm6
G4TARGET := $(name)
G4EXLIB := true

	@echo CPPFLAGS = $(CPPFLAGS)
	@echo LOADLIBS = $(LOADLIBS)

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

#hbu****************** do not define G4NOHIST to get histogramming
#G4NOHIST := true
ifdef G4NOHIST
  CPPFLAGS += -DG4NOHIST
else
  LDFLAGS  += -L$(CERN)/pro/lib
  LOADLIBS += -lpacklib $(FCLIBS)
endif

#hbu add verbose
# CPPFLAGS += -v
# LOADLIBS += -v

#hbu add also this dir to get access to extra header files
# CPPFLAGS += $(G4INSTALL)/examples/extended/electromagnetic/TestEm6
# CPPFLAGS += -I $(HOME)/public/geant4/examples/extended/electromagnetic/TestEm6

#hbu
info:
	@echo CPPFLAGS = $(CPPFLAGS)
	@echo LOADLIBS = $(LOADLIBS)
	@echo HOME = $(HOME)


include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
