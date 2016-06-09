# $Id: GNUmakefile,v 1.2 2005/04/22 16:26:46 gcosmo Exp $
#
# GNUmakefile: makefile for example ROOT I/O implementation.
#
#  $(name) is the name of granular library.
#
#  Author: Youhei Morita <youhei.morita@kek.jp>
#

# name of the granular library
name := G4PersEx01

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

CPPFLAGS += -Iinclude
CPPFLAGS += -I$(G4INCLUDE)
CPPFLAGS += -I$(G4PEX_DIR)/include

include $(G4PEX_DIR)/config/binmake.gmk

