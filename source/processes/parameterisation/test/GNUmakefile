# $Id: GNUmakefile,v 1.1 2002-11-02 00:10:35 mverderi Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk


include $(G4INSTALL)/config/binmake.gmk
