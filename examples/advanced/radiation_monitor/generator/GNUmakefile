# $Id: GNUmakefile,v 1.1 2005/10/25 16:37:59 capra Exp $
# --------------------------------------------------------------------
# GNUmakefile for exaple subpackage
# R. Capra, 09 Sep 2005
# --------------------------------------------------------------------

name := Radmon

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/run/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/processes/cuts/include \
            -I$(G4BASE)/event/include \
            -I../generalPurpose/include

include $(G4INSTALL)/config/common.gmk
