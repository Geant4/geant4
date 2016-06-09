# $Id: GNUmakefile,v 1.6 2004/06/09 15:04:37 daquinog Exp $
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
 CPPFLAGS += `aida-config --include`
# LDFLAGS  += `aida-config --lib`
endif


CPPFLAGS += -I$(PYTHON_INCLUDE_DIR) \
	    -I$(TIARA_BASE)/source/tiaraPhysicsPackaging/include \
	    -I$(TIARA_BASE)/source/tiaraPhysicsLists/include 

LDFLAGS += -L$(PYTHON_LIB_DIR)

ifdef G4ANALYSIS_USE
 LDLIBS += `aida-config --lib`
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
