# $Id: GNUmakefile,v 1.2 2002-06-10 12:00:14 radoone Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := NO3gogdml
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

GDMLINSTALL := $(G4INSTALL)/examples/extended/gdml
GDMLLIBROOT := $(GDMLINSTALL)/libgdml

CPPFLAGS    += \
               -I$(GDMLINSTALL)/libgdml/saxana/include \
               -I$(GDMLINSTALL)/libgdml/schema/include \
               -I$(GDMLINSTALL)/libgdml/management/include \
               -I$(GDMLINSTALL)/libgdml/processes/include \
               -I$(GDMLINSTALL)/libgdml/subscribers/include \
               -I$(XERCESCROOT)/include

EXTRALIBS   += \
              -lGDMLsaxana \
              -lGDMLmanagement \
              -lGDMLschema \
              -lGDMLprocesses \
              -lGDMLsubscribers \
              -L $(XERCESCROOT)/lib -lxerces-c
          
EXTRALIBSSOURCEDIRS += \
                       $(GDMLLIBROOT)/saxana \
                       $(GDMLLIBROOT)/schema \
                       $(GDMLLIBROOT)/processes \
                       $(GDMLLIBROOT)/management \
                       $(GDMLLIBROOT)/subscribers
                       
EXTRA_LINK_DEPENDENCIES := \
                         $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)/libNO3gogdml.a \
                         $(G4WORKDIR)/lib/$(G4SYSTEM)/libGDMLmanagement.a \
                         $(G4WORKDIR)/lib/$(G4SYSTEM)/libGDMLprocesses.a \
                         $(G4WORKDIR)/lib/$(G4SYSTEM)/libGDMLsaxana.a \
                         $(G4WORKDIR)/lib/$(G4SYSTEM)/libGDMLschema.a \
                         $(G4WORKDIR)/lib/$(G4SYSTEM)/libGDMLsubscribers.a


.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
