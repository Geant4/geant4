# $Id: GNUmakefile,v 1.1 2004-05-12 09:38:33 gunter Exp $
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

EXTRALIBDIR := -L/cern/pro/lib
EXTRALIBS := -lpacklib $(FCLIBS) 
  ifeq ($(G4SYSTEM),SUN-CC)
    EXTRALIBS := -lpacklib -lw -lgen -lsocket -lnsl -lintl $(FCLIBS) 
  endif
LDFLAGS += $(EXTRALIBDIR)
LOADLIBS += $(EXTRALIBS)

include $(G4INSTALL)/config/binmake.gmk
    










