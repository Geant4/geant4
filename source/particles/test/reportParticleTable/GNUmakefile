# $Id: GNUmakefile,v 1.2 2003-09-21 19:38:02 kurasige Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

# $Id: GNUmakefile,v 1.2 2003-09-21 19:38:02 kurasige Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

name := reportParticles
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

# Override some variables for binmake.gmk.
#
# INCFLAGS :=
# LDLIBS   :=

include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += \
            -I$(G4BASE)/particles/utils/include 
