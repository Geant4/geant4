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

ifdef ROOTSYS
  ROOTCONFIG = $(ROOTSYS)/bin/root-config
  ROOTCINT = $(ROOTSYS)/bin/rootcint
else
  ROOTCONFIG = root-config
  ROOTCINT = rootcint 
endif

# root
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
ROOTLIBS     := $(filter-out -lpthread,$(ROOTLIBS))
ROOTGLIBS    := $(filter-out -lpthread,$(ROOTGLIBS))
ROOTINCLUDES := $(shell $(ROOTCONFIG) --cflags)


SRC_PRECO        := src/precoreaction.$(SRCSuf) \
                    src/precofragment_base.$(SRCSuf) \
                    src/precodict.$(SRCSuf)

OBJ_PRECO        := $(patsubst %.${SRCSuf},%.${OBJSuf},$(subst src/,$(PRECOMPOUND_TMP)/,${SRC_PRECO}))
DEP_PRECO        := $(patsubst %.${OBJSuf},%.${DEPSuf},${OBJ_PRECO})
LIB_PRECO        := $(PRECOMPOUND_TMP)/../libPreco.$(LIBSuf)


.SUFFIXES:
.SUFFIXES: .${DEPSuf} .${SRCSuf} .${OBJSuf} .${LIBSuf}

.PHONY: all print clean clear

all: ${LIB_PRECO}


${LIB_PRECO}: ${OBJ_PRECO}
ifdef CPPVERBOSE
	${LD} ${SOFLAGS} ${LDFLAGS} ${OBJ_PRECO} -o ${LIB_PRECO}
else
	@echo "Making shared library $@..."	
	@${LD} ${SOFLAGS} ${LDFLAGS} ${OBJ_PRECO} -o ${LIB_PRECO}
endif

${OBJ_PRECO}: ${PRECOMPOUND_TMP}/%.o: src/%.cc
ifdef CPPVERBOSE
	${CXX} ${CXXFLAGS} ${LOCALINCLUDES} ${ROOTINCLUDES}  -o $@  -c $<
else
	@echo "Compiling $<..."
	@${CXX} ${CXXFLAGS} ${LOCALINCLUDES} ${ROOTINCLUDES}  -o $@  -c $<
endif

${DEP_PRECO}: ${PRECOMPOUND_TMP}/%.d: src/%.cc
ifdef CPPVERBOSE
	${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} $< | sed 's!$*.o!${PRECOMPOUND_TMP}/& $@!' >$@
else
	@echo "Making dependency for file $<..."
	@${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} $< | sed 's!$*.o!${PRECOMPOUND_TMP}/& $@!' >$@
endif

-include ${DEP_PRECO}


src/precodict.cc: include/precoreaction.h include/precofragment.h \
	                  include/precofragment_base.h include/precolinkdef.h
ifdef CPPVERBOSE
	rm -f src/precodict.* include/precodict.*
	ln -s include/precoreaction.h .
	ln -s include/precofragment.h .
	ln -s include/precofragment_base.h .
	$(ROOTCINT) src/precodict.cc -c precofragment_base.h precofragment.h precoreaction.h include/precolinkdef.h
	rm precoreaction.h precofragment.h precofragment_base.h
	mv src/precodict.h include/precodict.h
else 
	@echo "Generating dictionary $@..."
	@rm -f src/precodict.* include/precodict.*
	@ln -s include/precoreaction.h .
	@ln -s include/precofragment.h .
	@ln -s include/precofragment_base.h .
	@$(ROOTCINT) src/precodict.cc -c precofragment_base.h precofragment.h precoreaction.h include/precolinkdef.h
	@rm precoreaction.h precofragment.h precofragment_base.h
	@mv src/precodict.h include/precodict.h
endif

clear:
	@rm -f ./*~
	@rm -f include/*~
	@rm -f src/*~
	@rm -f ./core

clean:
	@rm -f ${OBJ_PRECO}
	@rm -f ${DEP_PRECO}
	@rm -f ${LIB_PRECO}
	@rm -f include/precodict.h 
	@rm -f src/precodict.cc


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
	@echo "SRC_PRECO = ${SRC_PRECO}"
	@echo "OBJ_PRECO = ${OBJ_PRECO}"
	@echo "DEP_PRECO = ${DEP_PRECO}"
	@echo "LIB_PRECO = ${LIB_PRECO}"
	@echo "----------------------------------------------------------"
