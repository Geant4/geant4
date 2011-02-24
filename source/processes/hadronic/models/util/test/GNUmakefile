# $Id: GNUmakefile,v 1.2 2008-05-15 14:05:36 gunter Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../../../
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk
	


# ---- For CLHEP histograms

#EXTRALIBDIR := -L/cern/pro/lib
#EXTRALIBS := -lpacklib $(FCLIBS) 
#LDFLAGS += $(EXTRALIBDIR)
#LOADLIBS += $(EXTRALIBS)

include $(G4INSTALL)/config/binmake.gmk
    










