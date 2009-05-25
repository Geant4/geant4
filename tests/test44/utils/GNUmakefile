# $Id: GNUmakefile,v 1.2 2009-05-25 16:32:41 bagoulia Exp $
# --------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := reader_test44
endif

name := ${G4TARGET}

G4EXEC_BUILD  = true

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS  += $(shell $(ROOTSYS)/bin/root-config --glibs)

include $(G4INSTALL)/config/binmake.gmk

