# $Id: GNUmakefile,v 1.8 2001-12-04 13:23:51 flongo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := GammaRayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

#
# NEW: this requires the env-var LHCXX_REL_DIR to be set up correctly,
# see the first lines in G4ComptonTest.cc
#

ifdef G4ANALYSIS_USE

   CPPFLAGS  += -I$(LHCXX_REL_DIR)/include

ifdef G4ANALYSIS_USE_NTUPLE
  CPPFLAGS += -DG4ANALYSIS_USE_NTUPLE
endif

#   C_FLAGS 	 += -I$(LHCXX_REL_DIR)/include
   LDFLAGS  += -L$(LHCXX_REL_DIR)/lib
# general stuff:	
   LOADLIBS += -lHepUtilities -lAIDA_Annotation
# add for Histograms: 
   LOADLIBS += -lAIDA_HBook 
# add for Vectors:
   LOADLIBS += -lVector 
# add for Plotter:
   LOADLIBS += -lAIDA_Plotter -lqp -lqpview -lqt-mt -lqutil -lqt 
# add for Ntuples:
   LOADLIBS += -lAIDA_Ntuple -lNtupleTag -lNtupleTag_HBook 
# ... and finally system stuff and to deal with the "Fortran legacy":
   LOADLIBS += -lCHBook -ldl -lg2c 
# Xm Xt
   LDFLAGS  += -L$/usr/local/lib -lXm -L/usr/X11R6/lib -lXt



endif

include $(G4INSTALL)/config/binmake.gmk





