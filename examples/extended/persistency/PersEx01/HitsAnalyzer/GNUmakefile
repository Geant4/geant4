# $Id: GNUmakefile,v 1.1 2002-12-04 14:51:34 morita Exp $
#
# GNUmakefile: makefile for example ROOT I/O impelmentatoin.
#
#  This makefile creates a granular library in $(FADS_LIB).
#  $(name) is the name of granular library.
#
#  Author: Youhei Morita <youhei.morita@kek.jp>
#

# name of the granular library
name := HitsAnalyzer

ifndef G4PEX_DIR
  G4PEX_DIR := $(G4INSTALL)/examples/extended/persistency/rootio
endif

include $(G4PEX_DIR)/config/G4PEXsetup.gmk

cwd := $(shell pwd)

CINT_INCLUDES := -I$(cwd)/include
CINT_INCLUDES += -I$(G4INCLUDE)
CINT_INCLUDES += -I$(G4PEX_DIR)/include
CINT_INCLUDES += -I$(CLHEP_BASE_DIR)/include
CINT_INCLUDES += -I$(ROOTSYS)/include

CPPFLAGS += -I../include
CPPFLAGS += -I$(G4INCLUDE)
CPPFLAGS += -I$(G4PEX_DIR)/include
CPPFLAGS += -I$(ROOTSYS)/include

EXTLDFLAGS := -Wl,-rpath $(G4PEX_DIR)/lib/$(G4SYSTEM) \
              -L$(G4PEX_DIR)/lib/$(G4SYSTEM) \
              -lG4PersEx01

include $(G4PEX_DIR)/config/binmake.gmk

