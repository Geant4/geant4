# $Id: GNUmakefile,v 1.7 2009-06-15 10:35:12 gunter Exp $
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
	
CPPFLAGS += `root-config --cflags`
LOADLIBS += `root-config --ldflags --libs`


include $(G4INSTALL)/config/binmake.gmk
    
INCFLAGS += -I$(G4BASE)/processes/hadronic/models/im_r_matrix/include
INCFLAGS += -I$(G4BASE)/processes/hadronic/models/binary_cascade/include










