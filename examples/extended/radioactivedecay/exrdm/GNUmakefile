# $Id: GNUmakefile,v 1.2 2004/12/13 16:06:02 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module
# --------------------------------------------------------------

name := exRDM
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

# PhysicsList use hadronic-lists, so it needs below include command.

ifndef G4LISTS_BASE
  EXTRALIBS += -L$(G4LIB)/plists/$(G4SYSTEM)
  G4LISTS_BASE = $(G4INSTALL)/physics_lists
else
  EXTRALIBS += -L$(G4LISTS_BASE)/hadronic/plists/lib/$(G4SYSTEM)
endif

# Select your physics lists to link against.
#
EXTRALIBS += -lQGSP_BERT
EXTRALIBS += -lQGSP_BIC
EXTRALIBS += -lQGSP_HP
EXTRALIBS += -lLHEP_BIC
EXTRALIBS += -lLHEP_BIC_HP
EXTRALIBS += -lLHEP_BERT
EXTRALIBS += -lLHEP_BERT_HP
EXTRALIBS += -lPackaging

CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/Packaging/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP_BERT/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP_BIC/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_BERT/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_BERT_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_BIC/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_BIC_HP/include

include $(G4INSTALL)/config/binmake.gmk

#CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
#CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

#

ifdef G4ANALYSIS_USE
   CPPFLAGS += `aida-config --include`
   LOADLIBS += `aida-config --lib`
endif
