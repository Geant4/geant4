# $Id: GNUmakefile,v 1.3 2002-06-11 15:48:17 radoone Exp $
# --------------------------------------------------------------
# GNUmakefile for GDML sample module. Radovan Chytracek, 30/05/2002.
# --------------------------------------------------------------

name := gogdml
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
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
                       
# EXTRA_LINK_DEPENDENCIES := \
#                          $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)/libgogdml.a \
#                          $(G4WORKDIR)/lib/$(G4SYSTEM)/libGDMLmanagement.a \
#                          $(G4WORKDIR)/lib/$(G4SYSTEM)/libGDMLprocesses.a \
#                          $(G4WORKDIR)/lib/$(G4SYSTEM)/libGDMLsaxana.a \
#                          $(G4WORKDIR)/lib/$(G4SYSTEM)/libGDMLschema.a \
#                          $(G4WORKDIR)/lib/$(G4SYSTEM)/libGDMLsubscribers.a


.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
