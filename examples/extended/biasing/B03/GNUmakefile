# $Id: GNUmakefile,v 1.7 2002-11-18 13:22:48 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
#
# The variable ANAPHE_BASE has to be set before using this make file.
#

name := exampleB03
G4TARGET := $(name)
G4EXLIB := true

# you have to set PYTHON_DIR to the directory where python2.2 is installed
PYTHON_DIR:=/afs/cern.ch/sw/lhcxx/specific/redhat72/gcc-2.95.2/PublicDomainPackages/2.0.0


ifndef G4INSTALL
  G4INSTALL = ../../..
endif


.PHONY: all
all: lib _B03App.so bin
include $(G4INSTALL)/config/architecture.gmk


CPPFLAGS += -I$(PYTHON_DIR)/include/python2.2/



LDFLAGS += \
	   -L$(PYTHON_DIR)/lib/python2.2/config \
	   -L/lib \
	   -lnsl -ldl -lreadline /lib/libtermcap.so.2  -lieee -lpthread -lutil -lpython2.2

include $(G4INSTALL)/config/binmake.gmk

#solibs = $(shell /bin/ls $(G4WORKDIR)/lib/$(G4SYSTEM)/*.so)

libdir = $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
objectfiles = $(shell /bin/ls $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)/*.o)
_B03App.so:
	@echo Creating shared library: $(libdir)/$@  
	$(CXX) -Wl,-soname,$(@F) -shared $(LDFLAGS) $(LDLIBS) -L$(libdir) $(objectfiles) -o $(libdir)/$@ 


