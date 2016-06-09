# $Id: GNUmakefile,v 1.4 2003/12/04 11:08:48 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile   Michael Dressel
# --------------------------------------------------------------

name := tiaraPhysicsPackaging
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif
include $(G4INSTALL)/config/architecture.gmk
include $(TIARA_BASE)/source/hadronic_lists.gmk

.PHONY: all tiaraClean

all: $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)/lib$(name)-handMade.so


include $(G4INSTALL)/config/binmake.gmk


# this is required because the normal G4 way of building the library
# does not take all the dependencies into account
libdir = $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
srcFiles = $(wildcard src/*.cc)
TMP = $(subst .cc,.o,$(srcFiles)) 
objectfiles = $(subst src/,$(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)/,$(TMP))

MYLIBS = $(subst -l$(name),,$(LDLIBS)) 	
$(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)/lib$(name)-handMade.so: $(objectfiles)
	@echo Creating shared library: $(libdir)/$@  
	$(CXX) -Wl,-soname,$(@F) -shared $(LDFLAGS) $(MYLIBS) -L$(libdir) $(objectfiles) -lc -o $@ 

tiaraClean:
	@echo Cleaning up $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
	@rm -rf $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)

