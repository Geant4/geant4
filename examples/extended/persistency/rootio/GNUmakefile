# $Id: GNUmakefile,v 1.1 2002-12-04 02:44:27 morita Exp $
#
# GNUmakefile: makefile for example ROOT I/O impelmentatoin.
#
#  This makefile creates a granular library in $(FADS_LIB).
#  $(name) is the name of granular library.
#
#  Author: Youhei Morita <youhei.morita@kek.jp>
#

# name of the granular library
name := G4ROOTIO

ifndef G4PEX_DIR
  G4PEX_DIR := $(G4INSTALL)/examples/extended/persistency/rootio
endif

include $(G4PEX_DIR)/config/G4PEXsetup.gmk

CPPFLAGS += -I$(G4INSTALL)/source/global/management/include \
            -I$(G4INSTALL)/source/global/HEPGeometry/include \
            -I$(G4INSTALL)/source/digits+hits/hits/include \
            -I$(G4INSTALL)/source/digits+hits/digits/include \
            -I$(G4INSTALL)/source/particles/management/include \
            -I$(G4INSTALL)/source/run/include \
            -I$(G4INSTALL)/source/event/include \
            -I$(G4INSTALL)/source/track/include \
            -I$(G4INSTALL)/source/tracking/include \
            -I$(G4INSTALL)/source/intercoms/include \
            -I$(G4INSTALL)/source/geometry/management/include \
            -I$(G4INSTALL)/source/persistency/include

# check HepMC
##hepmc_barcode := $(shell $(G4PEX_DIR)/config/hepmctest.csh)
##ifeq ($(hepmc_barcode),false)  
  CXXFLAGS += -DHEPMC_HAS_NO_BARCODE
##endif

include $(G4PEX_DIR)/config/common.gmk

