# $Id: GNUmakefile,v 1.7 2004/06/09 15:04:35 daquinog Exp $
# --------------------------------------------------------------
# GNUmakefile   Michael Dressel
# --------------------------------------------------------------

name := TiaraWrapper
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif
include $(G4INSTALL)/config/architecture.gmk


SWIGFLAGS = -I../G4KernelWrapper \
	-I$(G4INSTALL)/source/geometry/biasing/include \
	-I$(G4INSTALL)/source/digits_hits/detector/include \
	-I$(G4INSTALL)/source/event/include \
	-I$(G4INSTALL)/source/processes/transportation/include \
	-I$(TIARA_BASE)/source/tiara/include \
	-I$(TIARA_BASE)/source/tiaraPhysicsLists/include \
	-I$(TIARA_BASE)/source/CLHEPWrapper/include \
	-I$(TIARA_BASE)/source/G4KernelWrapper/include \
	-I${CLHEP_BASE_DIR}/include 

#include $(TIARA_BASE)/source/hadronic_lists.gmk
#SWIGFLAGS += ${LSTFLAGS}

libdirTiara = $(G4WORKDIR)/tmp/$(G4SYSTEM)/tiara
SWIG_OBJ_DEP = $(wildcard $(libdirTiara)/*.o)


include $(TIARA_BASE)/source/swig.gmk

tiaraDeps = $(notdir $(wildcard ${G4WORKDIR}/tmp/${G4SYSTEM}/tiara/*.d))
copyDeps = $(addprefix ${G4WORKDIR}/tmp/${G4SYSTEM}/TiaraWrapper/, $(tiaraDeps)) 


.PHONY: all tiaraClean

all: $(copyDeps) $(SWIGSO)



include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE
 CPPFLAGS += `aida-config --include`
# LDFLAGS  += `aida-config --lib`
endif

CPPFLAGS += -I$(TIARA_BASE)/source/tiara/include\
	    -I$(PYTHON_INCLUDE_DIR) \
	    -I$(TIARA_BASE)/source/tiaraPhysicsPackaging/include \
	    -I$(TIARA_BASE)/source/tiaraPhysicsLists/include \
	    -I$(TIARA_BASE)/source/G4KernelWrapper/include

LDFLAGS +=  -L$(G4WORKDIR)/tmp/$(G4SYSTEM)/tiara \
	    -L$(G4WORKDIR)/tmp/$(G4SYSTEM)/tiaraPhysicsPackaging \
	    -L$(G4WORKDIR)/tmp/$(G4SYSTEM)/tiaraPhysicsLists \
	    -L$(PYTHON_LIB_DIR)

LDLIBS1 += -ltiaraPhysicsLists \
	   -ltiaraPhysicsPackaging-handMade

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




${G4WORKDIR}/tmp/${G4SYSTEM}/TiaraWrapper/%.d: ${G4WORKDIR}/tmp/${G4SYSTEM}/tiara/%.d
	@cp  $< $@



tiaraClean:
	@echo Cleaning up $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
	@echo $(tiaraDeps)
	@rm -rf $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
	@rm -f *.pyc
