# $Id: GNUmakefile,v 1.5 2002-11-13 08:59:21 mdressel Exp $
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
all: lib _B03App.so bin
include $(G4INSTALL)/config/architecture.gmk

ANAPHE_BASE := /opt/Anaphe/5.0.1/specific/redhat72/gcc-2.95.2

CPPFLAGS += -I$(ANAPHE_BASE)/PublicDomainPackages/2.0.0/include/python2.2/



LDFLAGS += \
	   -L$(ANAPHE_BASE)/PublicDomainPackages/2.0.0/lib/python2.2/config \
	   -L/lib \
	   -lpython2.2 \
	   -lnsl -ldl -lreadline /lib/libtermcap.so.2  -lieee -lpthread -lutil 

include $(G4INSTALL)/config/binmake.gmk

#solibs = $(shell /bin/ls $(G4WORKDIR)/lib/$(G4SYSTEM)/*.so)

libdir = $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
objectfiles = $(shell /bin/ls $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)/*.o)
_B03App.so:
	@echo Creating shared library: $(libdir)/$@  
	$(CXX) -Wl,-soname,$(@F) -shared $(LDFLAGS) $(LDLIBS) -L$(libdir) $(objectfiles) -o $(libdir)/$@ 


