# $Id: GNUmakefile,v 1.2 2000-08-02 09:24:31 sscherer Exp $
# -----------------------------------------------------------
# GNUmakefile for hadronic library.  Gabriele Cosmo, 18/9/96.
# -----------------------------------------------------------

name := G4hadronic_quark_md_body

ifndef G4INSTALL
  G4INSTALL = ../../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

G4TMPDIR = $(G4TMP)/$(G4SYSTEM)/$(name)

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/processes/hadronic/models/quark_molecular_dynamics/util/include \
           
include $(G4INSTALL)/config/common.gmk

