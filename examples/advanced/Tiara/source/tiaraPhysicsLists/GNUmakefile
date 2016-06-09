# $Id: GNUmakefile,v 1.2 2003/12/04 11:08:48 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile   Michael Dressel
# --------------------------------------------------------------

name := tiaraPhysicsLists
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all tiaraClean

all: lib 

include $(G4INSTALL)/config/architecture.gmk


CPPFLAGS += -I../tiaraPhysicsPackaging/include
LDFLAGS += -L$(G4WORKDIR)/tmp/Linux-g++/tiaraPhysicsPackaging -ltiaraPhysicsPackaging

include $(TIARA_BASE)/source/hadronic_lists.gmk

include $(G4INSTALL)/config/binmake.gmk

tiaraClean:
	@echo Cleaning up $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
	@rm -rf $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)


