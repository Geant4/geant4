# $Id: GNUmakefile,v 1.8 2004/12/14 13:11:07 daquinog Exp $
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

include $(TIARA_BASE)/source/hadronic_lists.gmk

.PHONY: all tiaraClean

all: lib bin

include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE
 CPPFLAGS += `aida-config --include`
# LDFLAGS  += `aida-config --lib`
endif


CPPFLAGS += -I$(PYTHON_INCLUDE_DIR) \
            -I$(G4INSTALL)/physics_lists/hadronic/LHEP_BIC/include \
            -I$(G4INSTALL)/physics_lists/hadronic/LHEP_BIC_HP/include \
	    -I$(G4INSTALL)/physics_lists/hadronic/LHEP/include \
	    -I$(G4INSTALL)/physics_lists/hadronic/LHEP_PRECO/include \
	    -I$(G4INSTALL)/physics_lists/hadronic/LHEP_LEAD/include \
	    -I$(G4INSTALL)/physics_lists/hadronic/LHEP_PRECO_HP/include

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
