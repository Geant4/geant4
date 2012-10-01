# $Id: GNUmakefile,v 1.4 2008-06-11 17:17:54 schaelic Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := Pol01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

plots:	plotResults.class
	@echo calling JAIDA ...
	java -cp .:$(CLASSPATH) plotResults Pol01.aida

plotResults.class: plotResults.java
	@echo compiling $<
	javac $<
