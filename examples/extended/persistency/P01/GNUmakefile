# $Id: GNUmakefile,v 1.4 2010-12-01 14:18:27 witoldp Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleP01
G4TARGET := $(name)
G4EXLIB := true

LOCALDIR := ${PWD}

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

ifndef PLATFORM
  PLATFORM = slc4_ia32_gcc34
endif

CPPFLAGS += -I$(ROOTSYS)/include
#EXTRALIBS := -L$(ROOTSYS)/lib -lCore -lTree -lCint -lReflex -lCintex -ldl
EXTRALIBS := $(shell $(ROOTSYS)/bin/root-config --libs) -lReflex -lCintex

DICTIONARYPATH := ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/dictionary

.PHONY: all
all: dictionary lib bin

dictionary:
ifndef ROOTSYS
	@echo ROOTSYS not defined!
endif
ifndef GCCXMLPATH
	@echo GCCXMLPATH not defined!
endif
	mkdir -p ${DICTIONARYPATH}
	${ROOTSYS}/bin/genreflex ${LOCALDIR}/include/Classes.hh \
	-s ${LOCALDIR}/xml/selection.xml -o ${DICTIONARYPATH}/Classes_rflx.cpp \
	--gccxmlpath ${GCCXMLPATH} -Iinclude -I${G4INCLUDE} -I${CLHEP_INCLUDE_DIR}
#
	gcc -o ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/libClassesDict.so -fpic -shared \
	${DICTIONARYPATH}/Classes_rflx.cpp src/*.cc \
	-Iinclude -I${ROOTSYS}/include -I${G4INCLUDE} -I${CLHEP_INCLUDE_DIR} \
	${EXTRALIBS} -L${CLHEP_LIB_DIR} -l${CLHEP_LIB} -L${G4LIB}/${G4SYSTEM} ${LDLIBS2}
#
	mkdir -p ${G4WORKDIR}/bin/${G4SYSTEM}
	gcc -o ${G4WORKDIR}/bin/${G4SYSTEM}/readHits readHits.cc \
	-Iinclude -I${ROOTSYS}/include -I${G4INCLUDE} -I${CLHEP_INCLUDE_DIR} \
	${EXTRALIBS} -L${CLHEP_LIB_DIR} -l${CLHEP_LIB} -L${G4LIB}/${G4SYSTEM} ${LDLIBS2}

clean_all: clean
	@$(RM) hits.root

include $(G4INSTALL)/config/binmake.gmk
