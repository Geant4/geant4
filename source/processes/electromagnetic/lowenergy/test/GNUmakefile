# $Id: GNUmakefile,v 1.12 2002-07-31 16:11:47 vnivanch Exp $
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
#      VI - will not run without in AIDA 4.0.4
OGLLIBS   := -L$(OGLHOME)/lib -lGLU -lGL
LOADLIBS += $(OGLLIBS)
# ----
#
# NEW: this requires the env-var LHCXX_REL_DIR to be set up correctly,
# see the first lines in G4ComptonTest.cc
#
CPPFLAGS += `aida-config --useHBook --include`
LDFLAGS  += `aida-config --useHBook --lib`

include $(G4INSTALL)/config/binmake.gmk




