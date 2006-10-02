# $Id: GNUmakefile,v 1.1 2006-10-02 13:48:10 vnivanch Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestPol2
G4TARGET := $(name)
G4EXLIB := true

# compile with polarization library
CPPFLAGS += -I../TestPol1/include
EXTRALIBS += -L$(G4WORKDIR)/tmp/$(G4SYSTEM)/TestPol1 -lTestPol1

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE
  # for the aida-config command see the README file
  CPPFLAGS += `aida-config --include`
  LDLIBS  += `aida-config --lib`  

plots:	plotResults.class
	@echo calling JAIDA ...
	java -cp .:$(CLASSPATH) plotResults testPol2.aida

plotResults.class: plotResults.java
	@echo compiling $<
	javac $<
endif
