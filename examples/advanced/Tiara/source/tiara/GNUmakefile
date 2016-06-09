# $Id: GNUmakefile,v 1.3 2003/06/16 15:36:54 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile   Michael Dressel
# --------------------------------------------------------------

name := tiara
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif
include $(G4INSTALL)/config/architecture.gmk

.PHONY: all tiaraClean

all: lib bin




include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE
CPPFLAGS += `$(ANAPHE_SCRIPTS)/aida-config -i`
endif

CPPFLAGS += -I$(PYTHON_INCLUDE_DIR) \
	    -I$(TIARA_BASE)/source/tiaraPhysicsPackaging/include \
	    -I$(TIARA_BASE)/source/tiaraPhysicsLists/include 

LDFLAGS += -L$(G4WORKDIR)/tmp/$(G4SYSTEM)/tiaraPhysicsPackaging \
	    -L$(G4WORKDIR)/tmp/$(G4SYSTEM)/tiaraPhysicsLists \
	    -L$(PYTHON_LIB_DIR)

LDLIBS1 +=  -ltiaraPhysicsLists \
	    -ltiaraPhysicsPackaging-handMade

ifdef G4ANALYSIS_USE
LDLIBS +=   `$(ANAPHE_SCRIPTS)/aida-config -l`
endif

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
