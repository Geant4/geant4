# --------------------------------
# GNUmakefile for examples module.
# --------------------------------

name := examplePar02
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

ifndef XERCESCROOT      
        @echo XERCESCROOT not defined!
endif

include $(G4INSTALL)/config/binmake.gmk
