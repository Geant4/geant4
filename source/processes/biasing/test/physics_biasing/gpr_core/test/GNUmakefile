# $Id: GNUmakefile,v 1.1 2007-07-27 22:13:08 tinslay Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../../
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
            -I$(G4BASE)/processes/biasing/test/physics_biasing/gpr_base/include \
            -I$(G4BASE)/processes/biasing/test/physics_biasing/gpr_core/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/intercoms/include

#INCFLAGS := -I$(G4BASE)/processes/biasing/test/physics_biasing/base/include \
#            -I$(G4BASE)/processes/management/include \
#	    -I$(G4BASE)/global/management/include \

include $(G4INSTALL)/config/binmake.gmk
