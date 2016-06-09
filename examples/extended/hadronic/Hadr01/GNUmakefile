# $Id: GNUmakefile,v 1.5 2006/06/23 10:46:54 vnivanch Exp $
# GEANT4 tag $Name: geant4-08-01 $
# --------------------------------------------------------------
# GNUmakefile for examples module with AIDA interface
# --------------------------------------------------------------

name := hadr01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../
endif

.PHONY: all
all: lib bin

ifdef G4ANALYSIS_USE
  CPPFLAGS += -DG4ANALYSIS_USE
endif

include $(G4INSTALL)/config/architecture.gmk

EXTRALIBS += -L$(G4LIB)/plists/$(G4SYSTEM)
EXTRALIBS += -lLHEP -lLHEP_BERT -lLHEP_BIC -lLHEP_BIC_HP -lLHEP_BERT_HP
EXTRALIBS += -lLHEP_PRECO -lLHEP_EMV  -lQGSP -lQGSP_EMV -lQGSP_EMX -lFTFC
EXTRALIBS +=  -lQGSC -lQGSP_BERT -lQGSP_BIC -lQGSP_BERT_HP -lQBBC -lPackaging

CPPFLAGS += \
            -I$(G4INSTALL)/physics_lists/hadronic/Packaging/include \
            -I$(G4INSTALL)/physics_lists/hadronic/LHEP/include \
            -I$(G4INSTALL)/physics_lists/hadronic/LHEP_BERT/include \
            -I$(G4INSTALL)/physics_lists/hadronic/LHEP_BIC/include \
            -I$(G4INSTALL)/physics_lists/hadronic/LHEP_BIC_HP/include \
            -I$(G4INSTALL)/physics_lists/hadronic/LHEP_BERT_HP/include \
            -I$(G4INSTALL)/physics_lists/hadronic/QGSC/include \
            -I$(G4INSTALL)/physics_lists/hadronic/QGSP/include \
            -I$(G4INSTALL)/physics_lists/hadronic/FTFC/include \
            -I$(G4INSTALL)/physics_lists/hadronic/LHEP_PRECO/include \
            -I$(G4INSTALL)/physics_lists/hadronic/LHEP_EMV/include \
            -I$(G4INSTALL)/physics_lists/hadronic/QGSP_BERT/include \
            -I$(G4INSTALL)/physics_lists/hadronic/QGSP_BERT_HP/include \
            -I$(G4INSTALL)/physics_lists/hadronic/QGSP_BIC/include \
            -I$(G4INSTALL)/physics_lists/hadronic/QGSP_EMV/include \
            -I$(G4INSTALL)/physics_lists/hadronic/QGSP_EMX/include \
            -I$(G4INSTALL)/physics_lists/hadronic/QBBC/include \

# for the aida-config command see the README file
ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
  LOADLIBS += `aida-config --lib`
endif

include $(G4INSTALL)/config/binmake.gmk

histclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/Histo.o

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

