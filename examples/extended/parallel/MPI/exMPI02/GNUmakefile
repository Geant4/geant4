# $Id: GNUmakefile,v 1.2 2010-01-11 16:55:08 gcosmo Exp $
# $Name: not supported by cvs2svn $
# ===========================================================
#   Makefile for building G4 MPI aplication
# ===========================================================

name := exMPI02
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
include ../mpi_interface/G4MPI.gmk

LDLIBS1 += -lG4UImpi

ifdef ROOTSYS
  CPPFLAGS += -I$(shell $(ROOTSYS)/bin/root-config --incdir)
  LDLIBS   += $(shell $(ROOTSYS)/bin/root-config --libs)
endif

ifeq ($(G4SYSTEM),Linux-g++)
  LDLIBS   += -Wl,-rpath,$(shell $(ROOTSYS)/bin/root-config --libdir)
endif

