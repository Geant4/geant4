# $Id: GNUmakefile,v 1.4 2002-11-08 17:35:17 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleB03
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif


.PHONY: all
all: lib B03Appc.so bin
include $(G4INSTALL)/config/architecture.gmk


CPPFLAGS += -I/home/dressel/python/local/include/python2.2/ 



LDFLAGS += \
	   -L/home/dressel/python/local/lib/python2.2/config \
	   -L/lib \
	   -lpython2.2 \
	   -lnsl -ldl -lreadline /lib/libtermcap.so.2 -lCLHEP-g++ -lieee -lpthread -lutil 

include $(G4INSTALL)/config/binmake.gmk

#solibs = $(shell /bin/ls $(G4WORKDIR)/lib/$(G4SYSTEM)/*.so)

libdir = $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
objectfiles = $(shell /bin/ls $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)/*.o)
B03Appc.so:
	@echo Creating shared library: $(libdir)/$@  
	$(CXX) -Wl,-soname,$(@F) -shared $(LDFLAGS) $(LDLIBS) -L$(libdir) $(objectfiles) -o $(libdir)/$@ 


