# $Id: GNUmakefile,v 1.2 2005/09/19 19:40:26 capra Exp $
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
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/graphics_reps/include \
            -I../generalPurpose/include

include $(G4INSTALL)/config/common.gmk
