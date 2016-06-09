# $Id: GNUmakefile,v 1.1.1.1 2003/06/12 13:08:24 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile   Michael Dressel
# --------------------------------------------------------------


name := CLHEPWrapper
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk


.PHONY: all tiaraClean

SWIGFLAGS = -I${CLHEP_BASE_DIR}/include
include $(TIARA_BASE)/source/swig.gmk


all: $(SWIGSO)



include $(G4INSTALL)/config/binmake.gmk



CPPFLAGS += -I$(PYTHON_INCLUDE_DIR)

LDFLAGS += -L$(PYTHON_LIB_DIR)

LDLIBS +=   -lpython$(PYTHONVERSION) \
	    -lnsl \
	    -ldl \
	    -lreadline \
	    -ltermcap \
	    -lieee \
	    -lpthread \
	    -lutil 


tiaraClean:
	@echo Cleaning up $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
	@rm -rf $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
