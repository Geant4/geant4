# $Id: GNUmakefile,v 1.9 2002-11-19 15:03:58 mdressel Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
#
# The variable PYTHON_DIR has to be setup. See ../README.
#

name := exampleB03
G4TARGET := $(name)
G4EXLIB := true


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
	   -lnsl -ldl -lreadline /lib/libtermcap.so.2  -lieee -lpthread -lutil -lexampleB03 -lpython2.2

include $(G4INSTALL)/config/binmake.gmk

libdir = $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
objectfiles = $(shell /bin/ls $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)/*.o)
_B03App.so:
	@echo Creating shared library: $(libdir)/$@  
	$(CXX) -Wl,-soname,$(@F) -shared $(LDFLAGS) $(LDLIBS) -L$(libdir) $(objectfiles) -o $(libdir)/$@ 


