# $Id: GNUmakefile,v 1.1.1.1 2003-06-12 13:08:24 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile   Michael Dressel
# --------------------------------------------------------------

name := tiara
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all tiaraClean

all: lib bin 


include $(G4INSTALL)/config/architecture.gmk


include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += `$(ANAPHE_SCRIPTS)/aida-config -i` \
	    -I$(PYTHON_INCLUDE_DIR) \
	    -I$(TIARA_BASE)/source/tiaraPhysicsPackaging/include \
	    -I$(TIARA_BASE)/source/tiaraPhysicsLists/include 

LDFLAGS += -L$(G4WORKDIR)/tmp/$(G4SYSTEM)/tiaraPhysicsPackaging \
	    -L$(G4WORKDIR)/tmp/$(G4SYSTEM)/tiaraPhysicsLists \
	    -L$(PYTHON_LIB_DIR)

LDLIBS1 +=  -ltiaraPhysicsLists \
	    -ltiaraPhysicsPackaging-handMade

LDLIBS +=   `$(ANAPHE_SCRIPTS)/aida-config -l` \
	    -lpython$(PYTHONVERSION) \
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
