# $Id: GNUmakefile,v 1.7 2005/11/25 01:56:26 capra Exp $
# --------------------------------------------------------------------
# GNUmakefile for exaple subpackage
# R. Capra, 09 Sep 2005
# --------------------------------------------------------------------

name := Radmon

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4UI_USE.gmk
include $(G4INSTALL)/config/G4VIS_USE.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/interfaces/basic/include \
            -I$(G4BASE)/interfaces/common/include \
            -I$(G4BASE)/run/include \
            -I$(G4BASE)/event/include \
            -I$(G4BASE)/tracking/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/solids/Boolean/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/cuts/include \
            -I$(G4BASE)/digits_hits/hits/include \
            -I$(G4BASE)/digits_hits/digits/include \
            -I$(G4BASE)/visualization/management/include \
            -I$(G4BASE)/graphics_reps/include \
            -I../generalPurpose/include \
            -I../generator/include \
            -I../detector/include \
            -I../materials/include \
            -I../physics/include \
            -I../analysis/include \
            -I../sensitivity/include \
            -I../userActions/include

ifdef G4ANALYSIS_USE  
 CPPFLAGS += `aida-config --include`
 CPPFLAGS += -DG4ANALYSIS_USE
 LDFLAGS  += `aida-config --lib`
 LOADLIBS += `aida-config --lib`
endif

include $(G4INSTALL)/config/common.gmk
