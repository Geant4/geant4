ARCH          = linux
OBJSuf        = o
SRCSuf        = cc
LIBSuf        = so
DEPSuf        = d

# Linux
ifeq ($(ARCH),linux)
  CXX           = g++
  LD            = g++
 ifdef PRECOMPOUND_DEBUG
    CXXFLAGS      = -g -Wall -fPIC
    LDFLAGS       = -g
 else
    CXXFLAGS      = -O2 -Wall -fPIC
    LDFLAGS       = -O2
  endif
  SOFLAGS       = -shared
endif

ifeq ($(CXX),)
  $(error $(ARCH) invalid architecture)
endif

LOCALINCLUDES := -I./include


ifndef G4INSTALL
  G4INSTALL = ../../../../../../../../..
endif
ifndef G4TMP
  G4TMP = ${G4INSTALL}/tmp
endif

PRECOMPOUND_TMP = ${G4TMP}/${G4SYSTEM}/G4PreCompoundTest/Utilities

DUMMY := $(shell [ ! -d $(PRECOMPOUND_TMP) ] && mkdir -p $(PRECOMPOUND_TMP) )


# root
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLIBS     := $(filter-out -lpthread,$(ROOTLIBS))
ROOTGLIBS    := $(filter-out -lpthread,$(ROOTGLIBS))
ROOTINCLUDES := $(shell root-config --cflags)


SRC_EXPDATA     := src/expdata.$(SRCSuf) \
		   src/expdata_unit.$(SRCSuf) \
                   src/expdata_deunit.$(SRCSuf) \
                   src/expdata_daunit.$(SRCSuf) \
                   src/expdata_ddunit.$(SRCSuf) \
                   src/expdata_ddaunit.$(SRCSuf) \
                   src/expdata_de.$(SRCSuf) \
                   src/expdata_da.$(SRCSuf) \
                   src/expdata_dd.$(SRCSuf) \
                   src/expdata_dda.$(SRCSuf) \
                   src/expdatadict.$(SRCSuf)

OBJ_EXPDATA      := $(patsubst %.${SRCSuf},%.${OBJSuf},$(subst src/,$(PRECOMPOUND_TMP)/,${SRC_EXPDATA})) 

DEP_EXPDATA      := $(patsubst %.${OBJSuf},%.${DEPSuf},${OBJ_EXPDATA})

LIB_EXPDATA      := $(PRECOMPOUND_TMP)/../libexpdata.$(LIBSuf)

.SUFFIXES:
.SUFFIXES: .${DEPSuf} .${SRCSuf} .${OBJSuf} .${LIBSuf}

.PHONY: all print clean clear

all: ${LIB_EXPDATA}


${LIB_EXPDATA}: ${OBJ_EXPDATA}
ifdef CPPVERBOSE
	${LD} ${SOFLAGS} ${LDFLAGS} ${OBJ_EXPDATA} -o $@
else
	@echo "Making shared library $@..."
	@${LD} ${SOFLAGS} ${LDFLAGS} ${OBJ_EXPDATA} -o $@
endif
 
${OBJ_EXPDATA}: ${PRECOMPOUND_TMP}/%.${OBJSuf}: src/%.${SRCSuf}
ifdef CPPVERBOSE
	${CXX} ${CXXFLAGS} ${LOCALINCLUDES} ${ROOTINCLUDES}  -o $@  -c $<
else
	@echo "Compiling $<..."
	@${CXX} ${CXXFLAGS} ${LOCALINCLUDES} ${ROOTINCLUDES}  -o $@  -c $<
endif


${DEP_EXPDATA}: ${PRECOMPOUND_TMP}/%.${DEPSuf}: src/%.${SRCSuf}
ifdef CPPVERBOSE
	${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} $< | sed 's!$*.o!${PRECOMPOUND_TMP}/& $@!' >$@
else
	@echo "Making dependency for file $<..."
	@${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} $< | sed 's!$*.o!${PRECOMPOUND_TMP}/& $@!' >$@
endif

-include ${DEP_EXPDATA}

src/expdatadict.cc: include/expdata.h include/expdata_unit.h include/expdata_deunit.h \
	            include/expdata_daunit.h include/expdata_ddunit.h include/expdata_ddaunit.h \
		    include/expdata_de.h include/expdata_da.h include/expdata_dd.h  \
		    include/expdata_dda.h include/expdataLinkDef.h 
ifdef CPPVERBOSE
	rm -f src/expdatadict.* include/expdatadict.*
	ln -s include/expdata.h .
	ln -s include/expdata.icc .
	ln -s include/expdata_unit.h .
	ln -s include/expdata_unit.icc .
	ln -s include/expdata_de.h .
	ln -s include/expdata_de.icc .
	ln -s include/expdata_da.h .
	ln -s include/expdata_da.icc .
	ln -s include/expdata_dd.h .
	ln -s include/expdata_dd.icc .
	ln -s include/expdata_dda.h .
	ln -s include/expdata_dda.icc .
	ln -s include/expdata_deunit.h .
	ln -s include/expdata_deunit.icc .
	ln -s include/expdata_daunit.h .
	ln -s include/expdata_daunit.icc .
	ln -s include/expdata_ddunit.h .
	ln -s include/expdata_ddunit.icc .
	ln -s include/expdata_ddaunit.h .
	ln -s include/expdata_ddaunit.icc .
	rootcint src/expdatadict.cc -c expdata.h expdata_unit.h \
	 expdata_deunit.h  expdata_daunit.h \
         expdata_ddunit.h expdata_ddaunit.h \
	 expdata_de.h expdata_da.h expdata_dd.h expdata_dda.h \
	 include/expdataLinkDef.h
	rm expdata.h expdata.icc \
	 expdata_unit.h expdata_unit.icc \
	 expdata_deunit.h expdata_deunit.icc \
	 expdata_daunit.h expdata_daunit.icc \
         expdata_ddunit.h expdata_ddunit.icc \
         expdata_ddaunit.h expdata_ddaunit.icc \
	 expdata_de.h expdata_de.icc \
	 expdata_da.h expdata_da.icc \
	 expdata_dd.h expdata_dd.icc \
	 expdata_dda.h expdata_dda.icc
	mv src/expdatadict.h include/expdatadict.h
else
	@rm -f src/expdatadict.* include/expdatadict.*
	@ln -s include/expdata.h .
	@ln -s include/expdata.icc .
	@ln -s include/expdata_unit.h .
	@ln -s include/expdata_unit.icc .
	@ln -s include/expdata_de.h .
	@ln -s include/expdata_de.icc .
	@ln -s include/expdata_da.h .
	@ln -s include/expdata_da.icc .
	@ln -s include/expdata_dd.h .
	@ln -s include/expdata_dd.icc .
	@ln -s include/expdata_dda.h .
	@ln -s include/expdata_dda.icc .
	@ln -s include/expdata_deunit.h .
	@ln -s include/expdata_deunit.icc .
	@ln -s include/expdata_daunit.h .
	@ln -s include/expdata_daunit.icc .
	@ln -s include/expdata_ddunit.h .
	@ln -s include/expdata_ddunit.icc .
	@ln -s include/expdata_ddaunit.h .
	@ln -s include/expdata_ddaunit.icc .
	@rootcint src/expdatadict.cc -c expdata.h expdata_unit.h \
	expdata_deunit.h  expdata_daunit.h \
        expdata_ddunit.h expdata_ddaunit.h \
	expdata_de.h expdata_da.h expdata_dd.h expdata_dda.h \
	include/expdataLinkDef.h        
	@rm expdata.h expdata.icc \
	 expdata_unit.h expdata_unit.icc \
	 expdata_deunit.h expdata_deunit.icc \
	 expdata_daunit.h expdata_daunit.icc \
	 expdata_ddunit.h expdata_ddunit.icc \
	 expdata_ddaunit.h expdata_ddaunit.icc \
	 expdata_de.h expdata_de.icc \
	 expdata_da.h expdata_da.icc \
	 expdata_dd.h expdata_dd.icc \
	 expdata_dda.h expdata_dda.icc
	@mv src/expdatadict.h include/expdatadict.h
endif


clear:
	@rm -f ./*~
	@rm -f include/*~
	@rm -f src/*~
	@rm -f ./core

clean: 
	@rm -f ${OBJ_EXPDATA}
	@rm -f ${DEP_EXPDATA}
	@rm -f ${LIB_EXPDATA}
	@rm -f include/expdatadict.h
	@rm -f src/expdatadict.cc

print:
	@echo "----------------------------------------------------------"
	@echo "ARCH     = ${ARCH}"
	@echo "CXX      = ${CXX}"
	@echo "CXXFLAGS = ${CXXFLAGS}"
	@echo "LD       = ${LD}"
	@echo "LDFLAGS  = ${LDFALGS}"
	@echo "SOFLAGS  = ${SOFLAGS}"
	@echo "----------------------------------------------------------"
	@echo "SRCSuf = ${SRCSuf}"
	@echo "OBJSuf = ${OBJSuf}"
	@echo "DEPSuf = ${DEPSuf}"
	@echo "LIBSuf = ${LIBSuf}"
	@echo "----------------------------------------------------------"
	@echo "ROOTCFLAGS    = ${ROOTCFLAGS}"
	@echo "ROOTLIBS      = ${ROOTLIBS}"
	@echo "ROOTGLIBS     = ${ROOTGLIBS}"
	@echo "ROOTLIBS      = ${ROOTLIBS}"
	@echo "ROOTGLIBS     = ${ROOTGLIBS}"
	@echo "ROOTINCLUDES  = ${ROOTINCLUDES}"
	@echo "LOCALINCLUDES = ${LOCALINCLUDES}"
	@echo "----------------------------------------------------------"
	@echo "SRC_EXPDATA = ${SRC_EXPDATA}"
	@echo "OBJ_EXPDATA = ${OBJ_EXPDATA}"
	@echo "DEP_EXPDATA = ${DEP_EXPDATA}"
	@echo "LIB_EXPDATA = ${LIB_EXPDATA}"
	@echo "----------------------------------------------------------"
