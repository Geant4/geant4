# $Id: GNUmakefile,v 1.3 2003-10-08 16:36:00 hpw Exp $
# -----------------------------------------------------------
# GNUmakefile for hadronic library.  Gabriele Cosmo, 18/9/96.
# -----------------------------------------------------------

name := G4hadronic_urqmd

G4DEBUG = 1

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

G4TMPDIR = $(G4TMP)/$(G4SYSTEM)/$(name)

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/STLInterface/g4rw \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/hadronic/management/include/ \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/processes/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
	    -I$(G4BASE)/processes/hadronic/models/de_excitation/include \
	    -I$(G4BASE)/processes/hadronic/models/management/include \
	    -I$(G4BASE)/processes/hadronic/models/util/include \
	    -I$(G4BASE)/processes/hadronic/models/im_r_matrix/include \
	    -I$(G4BASE)/processes/hadronic/models/pre_equilibrium/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/particles/shortlived/include \
            -I$(G4BASE)/materials/include
           
include $(G4INSTALL)/config/common.gmk

