# $Id: GNUmakefile,v 1.1 1999-01-07 16:04:58 gunter Exp $
# ----------------------------------------------------------------
# Generic GNUmakefile
# Set environment variable G4TARGET.
# ----------------------------------------------------------------
# TESTTARGET is created according to G4SYSTEM. G.Cosmo, 11/6/96.

MAKEFLAGS= --no-print-directory

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_USE.gmk

# Override some variables for binmake.gmk.
#
# INCFLAGS :=
# LDLIBS   :=

include $(G4INSTALL)/config/binmake.gmk

.PHONY: clean_all_tildas

clean_all_tildas:
	find $(G4BASE)/ -name '*~' -exec rm -i {} \;

