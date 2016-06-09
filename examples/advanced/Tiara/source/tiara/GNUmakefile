# $Id: GNUmakefile,v 1.5 2003/12/08 17:53:26 gcosmo Exp $
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

EXTRALIBS += -L$(G4WORKDIR)/tmp/$(G4SYSTEM)/tiaraPhysicsPackaging \
	     -L$(G4WORKDIR)/tmp/$(G4SYSTEM)/tiaraPhysicsLists \
             -ltiaraPhysicsLists \
	     -ltiaraPhysicsPackaging-handMade

include $(TIARA_BASE)/source/hadronic_lists.gmk

.PHONY: all tiaraClean

all: lib bin

include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE
  CPPFLAGS += `$(ANAPHE_SCRIPTS)/aida-config -i`
endif

CPPFLAGS += -I$(PYTHON_INCLUDE_DIR) \
	    -I$(TIARA_BASE)/source/tiaraPhysicsPackaging/include \
	    -I$(TIARA_BASE)/source/tiaraPhysicsLists/include 

LDFLAGS += -L$(PYTHON_LIB_DIR)

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
