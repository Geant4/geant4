# $Id: GNUmakefile,v 1.1 1998/10/14 15:18:56 allison Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := exampleE01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

# Override some variables for binmake.gmk.
#
# INCFLAGS :=
# LDLIBS   :=

CPPFLAGS  += \
             -I$(OBJY_DIR)/include \
             -I$(HEP_ODBMS_DIR)/include \
             -I$(HISTOODIR)/code
LDFLAGS   += \
             -L$(OBJY_DIR)/lib \
             -L$(HEP_ODBMS_DIR)/lib \
             -L$(HISTOODIR)/lib
LOADLIBS  += -lHistoxx -lHepODBMS -loo -lm

include $(G4INSTALL)/config/binmake.gmk

