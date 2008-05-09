# $Id: GNUmakefile,v 1.1 2008-05-09 10:33:46 vnivanch Exp $
# --------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := reader_test44
endif

name := ${G4TARGET}

G4EXEC_BUILD  = true

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(ROOTSYS)/include

EXTRALIBS += `root-config --glibs`

include $(G4INSTALL)/config/binmake.gmk

