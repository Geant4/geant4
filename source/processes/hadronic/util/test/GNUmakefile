# $Id: GNUmakefile,v 1.1 1999-01-08 16:33:47 gunter Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

# Override some variables for binmake.gmk.
#
INCFLAGS := \
	-I$(G4BASE)/processes/hadronic/util/include \
	-I$(G4BASE)/global/management/include
	
LDLIBS   := \
	    -lG4hadronic_util \
            -lG4global

LDFLAGS  += $(EXTRALIBDIR)
LOADLIBS += $(EXTRALIBS)

include $(G4INSTALL)/config/binmake.gmk
