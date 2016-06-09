# $Id: GNUmakefile 49035 2011-03-31 10:21:28Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := AnaEx02
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

#### G4ANALYSIS_USE := true

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

histclean:
	rm ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/HistoManager.o


ifdef ROOTSYS
  CPPFLAGS += -I$(shell $(ROOTSYS)/bin/root-config --incdir)
  CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
  LDLIBS += $(shell $(ROOTSYS)/bin/root-config --libs)

   ifeq ($(G4SYSTEM),Linux-g++)
     LDLIBS += -Wl,-rpath,$(shell $(ROOTSYS)/bin/root-config --libdir)
   endif
endif
