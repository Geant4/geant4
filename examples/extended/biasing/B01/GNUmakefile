# $Id: GNUmakefile,v 1.6 2002-10-22 14:09:02 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleB01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif


.PHONY: all
all: lib  bin


include $(G4INSTALL)/config/architecture.gmk


#CPPFLAGS += -I$(HOME)/G4scoring/scoring/include



#LDFLAGS += -L$(G4WORKDIR)/tmp/$(G4SYSTEM)/scoring 





include $(G4INSTALL)/config/binmake.gmk
#MYLIBS :=  -lexampleB01 -lscoring $(LDLIBS) -lXmu -lXm
MYLIBS :=  -lexampleB01 $(LDLIBS) -lXmu -lXm
LDLIBS = $(MYLIBS)

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*




