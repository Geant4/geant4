# ------------------------------------------------------------
# GNUmakefile for transportation library.  G.Folger 10-Dec-97.
# ------------------------------------------------------------

name := G4transportation

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
	    -I$(G4BASE)/global/management/include \
	    -I$(G4BASE)/global/HEPRandom/include \
	    -I$(G4BASE)/global/HEPGeometry/include \
	    -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/navigation/include \
	    -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/biasing/include \
	    -I$(G4BASE)/track/include \
	    -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/cuts/include \
	    -I$(G4BASE)/particles/management/include \
	    -I$(G4BASE)/materials/include \
	    -I$(G4BASE)/intercoms/include

ifdef G4DEBUG_TRANSPORT
  CPPFLAGS += -DG4DEBUG_TRANSPORT
endif
ifdef G4DEBUG_POSTSTEP_TRANSPORT
  CPPFLAGS += -DG4DEBUG_POSTSTEP_TRANSPORT
endif

include $(G4INSTALL)/config/common.gmk
