# $Id: GNUmakefile,v 1.1 1999-01-07 16:10:59 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for  processes/decay library.  G.Folger 9/12/97
# --------------------------------------------------------------

name := G4decay

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
	    -I$(G4BASE)/global/HEPRandom/include \
	    -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/magneticfield/include \
 	    -I$(G4BASE)/track/include \
	    -I$(G4BASE)/processes/management/include \
	    -I$(G4BASE)/particles/management/include \
	    -I$(G4BASE)/intercoms/include \
	    -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk
