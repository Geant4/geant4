# $Id: GNUmakefile,v 1.1.1.1 2003-06-12 13:08:24 dressel Exp $
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


SWIGFLAGS = -I$(G4INSTALL)/source/geometry/biasing/include \
	-I$(G4INSTALL)/source/digits+hits/detector/include \
	-I$(G4INSTALL)/source/event/include \
	-I$(G4INSTALL)/source/processes/transportation/include \
	-I$(TIARA_BASE)/source/tiara/include \
	-I$(TIARA_BASE)/source/tiaraPhysicsLists/include \
	-I$(TIARA_BASE)/source/CLHEPWrapper/include \
	-I$(TIARA_BASE)/source/G4KernelWrapper/include 



libdirTiara = $(G4WORKDIR)/tmp/$(G4SYSTEM)/tiara
SWIG_OBJ_DEP = $(wildcard $(libdirTiara)/*.o)

#$(shell /bin/ls $(libdirTiara)/*.o)


include $(TIARA_BASE)/source/swig.gmk

tiaraDeps = $(shell cd  ${G4WORKDIR}/tmp/${G4SYSTEM}/tiara/; ls *.d)
copyDeps = $(addprefix ${G4WORKDIR}/tmp/${G4SYSTEM}/TiaraWrapper/, $(tiaraDeps)) 


.PHONY: all tiaraClean

all: $(copyDeps) $(SWIGSO)



include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += -I$(TIARA_BASE)/source/tiara/include\
	     `$(ANAPHE_SCRIPTS)/aida-config -i` \
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

LDLIBS +=   `$(ANAPHE_SCRIPTS)/aida-config -l` \
	    -lpython$(PYTHONVERSION) \
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
	@rm -rf $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)

