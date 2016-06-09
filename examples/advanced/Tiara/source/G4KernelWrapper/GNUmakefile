# $Id: GNUmakefile,v 1.7 2006/06/23 13:15:34 ahoward Exp $
# --------------------------------------------------------------
# GNUmakefile   Michael Dressel
# --------------------------------------------------------------

name := G4KernelWrapper
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif
include $(G4INSTALL)/config/architecture.gmk

SWIGFLAGS = -I./ \
	-I${G4INSTALL}/source/global/management/include \
	-I${G4INSTALL}/source/geometry/biasing/include \
	-I${G4INSTALL}/source/digits_hits/detector/include \
	-I${G4INSTALL}/source/event/include \
	-I${G4INSTALL}/source/processes/transportation/include \
	-I${G4INSTALL}/source/visualization/management/include \
	-I${TIARA_BASE}/source/CLHEPWrapper/include \
	-I${CLHEP_BASE_DIR}/include

include $(TIARA_BASE)/source/swig.gmk

.PHONY: all tiaraClean

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
	@rm -f *.pyc
