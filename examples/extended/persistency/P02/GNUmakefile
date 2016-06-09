# $Id: GNUmakefile,v 1.2 2010-01-11 17:03:09 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleP02
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
EXTRALIBS := -L$(ROOTSYS)/lib -lCore -lTree -lCint -lReflex -lCintex -ldl -lNet -lRIO

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
	${ROOTSYS}/bin/genreflex ${LOCALDIR}/include/ExP02Classes.hh \
	-s ${LOCALDIR}/xml/selection.xml -o ${DICTIONARYPATH}/Classes_rflx.cpp \
	--gccxmlpath ${GCCXMLPATH} -Iinclude -I${G4INCLUDE} -I${CLHEP_INCLUDE_DIR}
#
	gcc -o ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/libClassesDict.so -fpic -shared \
	${DICTIONARYPATH}/Classes_rflx.cpp src/*.cc \
	-Iinclude -I${ROOTSYS}/include -I${G4INCLUDE} -I${CLHEP_INCLUDE_DIR} \
	${EXTRALIBS} -L${CLHEP_LIB_DIR} -l${CLHEP_LIB} -L${G4LIB}/${G4SYSTEM} ${LDLIBS2}
#
#	mkdir -p ${G4WORKDIR}/bin/${G4SYSTEM}
#	gcc -o ${G4WORKDIR}/bin/${G4SYSTEM}/readHits readHits.cc \
#	-Iinclude -I${ROOTSYS}/include -I${G4INCLUDE} -I${CLHEP_INCLUDE_DIR} \
#	${EXTRALIBS} -L${CLHEP_LIB_DIR} -l${CLHEP_LIB} -L${G4LIB}/${G4SYSTEM} ${LDLIBS2}

clean_all: clean
	@$(RM) hits.root

include $(G4INSTALL)/config/binmake.gmk
