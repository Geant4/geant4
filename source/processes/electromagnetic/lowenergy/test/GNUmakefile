# $Id: GNUmakefile,v 1.11 2002-06-18 14:08:34 pia Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif
# G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk
#
# ---- Not clear why I need this here... MGP, 18 June 2002
OGLLIBS   := -L$(OGLHOME)/lib -lGLU -lGL
LOADLIBS += $(OGLLIBS)
# ----
#
# NEW: this requires the env-var LHCXX_REL_DIR to be set up correctly,
# see the first lines in G4ComptonTest.cc
#
CPPFLAGS += `aida-config --useHBook --include`
LDFLAGS  += `aida-config --useHBook --lib`

# old:
#LDFLAGS  += -L/cern/pro/lib
#LOADLIBS += -lpacklib $(FCLIBS)

include $(G4INSTALL)/config/binmake.gmk




