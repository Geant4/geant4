# $Id: GNUmakefile,v 1.1 2000-03-23 16:45:50 hweber Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../../../../
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

# ---- For CLHEP histograms

#EXTRALIBDIR := -L/cern/pro/lib
#EXTRALIBS := -lpacklib $(FCLIBS) 
#  ifeq ($(G4SYSTEM),SUN-CC)
#    EXTRALIBS := -lpacklib -lw -lgen -lsocket -lnsl -lintl $(FCLIBS) 
#  endif
#LDFLAGS += $(EXTRALIBDIR)
#LOADLIBS += $(EXTRALIBS) 

#INCFLAGS := \
#            -I$(G4BASE)/processes/hadronic/models/generator/include \

G4LIBDIR := $(G4LIB)/$(G4SYSTEM)

EXTRA_LINK_DEPENDENCIES := $(G4INSTALL)/lib/SUN-CC/libG4procman.a
#	    -lG4event \
#	    -lG4hadronic_kinetic \
#	    -lG4hadronic_deex \
#	    -lG4hadronic_util_gen \
#	    -lG4hadronic_man_gen \
#            -lG4hadronic_util	 \
#	    -lG4procman \
#            -lG4bosons \
#            -lG4leptons \
#            -lG4shortlived \
#            -lG4baryons \
#            -lG4mesons \
#            -lG4ions \
#            -lG4partman \
#	    -lG4materials \
#	    -lG4track \
#            -lG4tracking \
#            -lG4intercoms \
#	    -lG4globman \
#	    -lG4hepnumerics

include $(G4INSTALL)/config/binmake.gmk






