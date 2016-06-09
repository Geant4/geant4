# $Id: GNUmakefile,v 1.1.1.1 2003/06/12 13:08:25 dressel Exp $
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

include $(G4INSTALL)/config/binmake.gmk

tiaraClean:
	@echo Cleaning up $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
	@rm -rf $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)


