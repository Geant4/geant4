# $Id: GNUmakefile,v 1.13 2002-11-28 13:17:18 pfeiffer Exp $
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
CPPFLAGS += `aida-config --include`
LDFLAGS  += `aida-config --lib`

include $(G4INSTALL)/config/binmake.gmk




