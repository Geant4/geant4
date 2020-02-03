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
	${ROOTSYS}/bin/genreflex ${LOCALDIR}/include/ExP02Classes.hh \
	-s ${LOCALDIR}/xml/selection.xml -o ${DICTIONARYPATH}/ExP02Classes_rflx.cpp \
	--gccxmlpath ${GCCXMLPATH} -Iinclude -I${G4INCLUDE} -I${CLHEP_INCLUDE_DIR}
#
	$(CXX) -o ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/libExP02ClassesDict.so -fpic -shared \
	${DICTIONARYPATH}/ExP02Classes_rflx.cpp src/*.cc \
	-Iinclude -I${ROOTSYS}/include -I${G4INCLUDE} -I${CLHEP_INCLUDE_DIR} \
	${EXTRALIBS} -L${G4LIB}/${G4SYSTEM} ${LDLIBS2} ${LOADLIBS}
#
#	mkdir -p ${G4WORKDIR}/bin/${G4SYSTEM}
#	gcc -o ${G4WORKDIR}/bin/${G4SYSTEM}/readHits readHits.cc \
#	-Iinclude -I${ROOTSYS}/include -I${G4INCLUDE} -I${CLHEP_INCLUDE_DIR} \
#	${EXTRALIBS} -L${G4LIB}/${G4SYSTEM} ${LDLIBS2} ${LDLIBS3} ${LOADLIBS}

clean_all: clean
	@$(RM) hits.root

include $(G4INSTALL)/config/binmake.gmk
