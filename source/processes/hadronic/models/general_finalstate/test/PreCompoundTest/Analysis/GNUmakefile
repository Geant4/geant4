
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
    CXXFLAGS      = -O2 -march=i686 -mcpu=i686 -Wall -fPIC
    LDFLAGS       = -O2
  endif
  SOFLAGS       = -shared
endif

ifeq ($(CXX),)
  $(error $(ARCH) invalid architecture)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../../../../../..
endif
ifndef G4TMP
  G4TMP = ${G4INSTALL}/tmp
endif

PRECOMPOUND_TMP = ${G4TMP}/${G4SYSTEM}/G4PreCompoundTest/Analysis

DUMMY := $(shell [ ! -d $(PRECOMPOUND_TMP) ] && mkdir -p $(PRECOMPOUND_TMP) )
DUMMY := $(shell [ ! -d ./bin ] && mkdir -p ./bin )

# root
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLIBS     := $(filter-out -lpthread,$(ROOTLIBS))
ROOTGLIBS    := $(filter-out -lpthread,$(ROOTGLIBS))
ROOTINCLUDES := $(shell root-config --cflags)

LOCALINCLUDES     := -I../Utilities/gzio/include -I../Utilities/preco/include -I../Utilities/expdata/include
LOCALLIBSDIR      := -L${PRECOMPOUND_TMP}/..
LOCALLIBPRECO     := -lPreco
LOCALLIBGZ        := -lz
LOCALLIBEXPDATA   := -lexpdata
LOCALLIBPLOTS     := -lplots
LOCALLIBMAKEPLOTS := -lmakeplots

LIB_PRECO   := ${PRECOMPOUND_TMP}/../libPreco.so
LIB_EXPDATA := ${PRECOMPOUND_TMP}/../libexpdata.so

#################### PLOTS LIB #############################################

SRC_PLOTS       := src/TranslateFragment.$(SRCSuf) \
                   src/mkplexpfilemanager.$(SRCSuf) \
                   src/mkplsimfilemanager.$(SRCSuf) \
                   src/mkplcomparisonhistograms.$(SRCSuf) \
                   src/mkpltesthistograms.$(SRCSuf) \
                   src/mkplsimmanager.$(SRCSuf) \
                   src/mkplexpdatamanager.$(SRCSuf) \
                   src/mkplmainframe.$(SRCSuf) \
                   src/mkplmakeplots.$(SRCSuf) \
                   src/mkplenergyselector.$(SRCSuf) \
                   src/mkplangleselector.$(SRCSuf) \
                   src/mkplrangeselector.$(SRCSuf) \
                   src/mkplcomparisonselector.$(SRCSuf) \
                   src/mkplvhistogram.$(SRCSuf) \
                   src/mkplde.$(SRCSuf) \
                   src/mkplda.$(SRCSuf) \
                   src/mkpldd.$(SRCSuf) \
                   src/mkpldda.$(SRCSuf) 


OBJ_PLOTS       := $(patsubst %.${SRCSuf},%.${OBJSuf},$(subst src/,${PRECOMPOUND_TMP}/,${SRC_PLOTS}))

DEP_PLOTS       := $(patsubst %.${OBJSuf},%.${DEPSuf},${OBJ_PLOTS})

LIB_PLOTS       := ${PRECOMPOUND_TMP}/../libmakeplots.$(LIBSuf)

################## MAKEPLOTS ###################################################

SRC_MAKEPLOTS   := makeplots.$(SRCSuf)

OBJ_MAKEPLOTS   := $(addprefix ${PRECOMPOUND_TMP}/,$(patsubst %.${SRCSuf},%.${OBJSuf},${SRC_MAKEPLOTS}))

DEP_MAKEPLOTS   := $(patsubst %.${OBJSuf},%.${DEPSuf},${OBJ_MAKEPLOTS})

EXE_MAKEPLOTS   := $(addprefix bin/,$(notdir $(basename ${OBJ_MAKEPLOTS})))

################## CONVERSION ##################################################

SRC_CONVERSION   := rconversion.$(SRCSuf)

OBJ_CONVERSION   := $(addprefix ${PRECOMPOUND_TMP}/,$(patsubst %.${SRCSuf},%.${OBJSuf},${SRC_CONVERSION}))

DEP_CONVERSION   := $(patsubst %.${OBJSuf},%.${DEPSuf},${OBJ_CONVERSION})

EXE_CONVERSION   := $(addprefix bin/,$(notdir $(basename ${OBJ_CONVERSION})))


################### SHOWDATA ###################################################

SRC_SHOWDATA      := showdata.$(SRCSuf)

OBJ_SHOWDATA      := $(addprefix ${PRECOMPOUND_TMP}/,$(patsubst %.${SRCSuf},%.${OBJSuf},${SRC_SHOWDATA}))

DEP_SHOWDATA      := $(patsubst %.${OBJSuf},%.${DEPSuf},${OBJ_SHOWDATA})

EXE_SHOWDATA      := $(addprefix bin/,$(notdir $(basename ${OBJ_SHOWDATA})))

################### GSHOWDATA ##################################################

SRC_GSHOWDATA     := gshowdata.$(SRCSuf)

OBJ_GSHOWDATA     := $(addprefix ${PRECOMPOUND_TMP}/,$(patsubst %.${SRCSuf},%.${OBJSuf},${SRC_GSHOWDATA}))

DEP_GSHOWDATA     := $(patsubst %.${OBJSuf},%.${DEPSuf},${OBJ_GSHOWDATA})

EXE_GSHOWDATA     := $(addprefix bin/,$(notdir $(basename ${OBJ_GSHOWDATA})))

################## GMAKEPLOTS ##################################################

SRC_GMAKEPLOTS   := gmakeplots.$(SRCSuf)

OBJ_GMAKEPLOTS   := $(addprefix ${PRECOMPOUND_TMP}/,$(patsubst %.${SRCSuf},%.${OBJSuf},${SRC_GMAKEPLOTS}))

DEP_GMAKEPLOTS   := $(patsubst %.${OBJSuf},%.${DEPSuf},${OBJ_GMAKEPLOTS})

EXE_GMAKEPLOTS   := $(addprefix bin/,$(notdir $(basename ${OBJ_GMAKEPLOTS})))


#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .${DEPSuf} .${SRCSuf} .${OBJSuf} .${LIBSuf}

.PHONY: all print clean clear clean

#all: ${LIB_PLOTS} ${EXE_MAKEPLOTS} ${EXE_CONVERSION}  ${EXE_SHOWDATA} ${EXE_GSHOWDATA} ${EXE_GMAKEPLOTS}
all: ${LIB_PLOTS}  ${EXE_CONVERSION} ${EXE_GSHOWDATA} ${EXE_GMAKEPLOTS}

# ------------------------------------------------------------------------------------

${LIB_PLOTS}: ${OBJ_PLOTS} ${OBJ_PLOTS_C}
ifdef CPPVERBOSE
	@echo "Making shared library $@..."
	${LD} ${SOFLAGS} ${LDFLAGS} ${OBJ_PLOTS} -o $@
else
	@echo "Making shared library $@..."
	@${LD} ${SOFLAGS} ${LDFLAGS} ${OBJ_PLOTS} -o $@
endif

${OBJ_PLOTS}: ${PRECOMPOUND_TMP}/%.${OBJSuf}: src/%.${SRCSuf}
ifdef CPPVERBOSE
	${CXX} ${CXXFLAGS} ${LOCALINCLUDES} -I./include ${ROOTINCLUDES}  -o $@  -c $<
else
	@echo "Compiling $<..."
	@${CXX} ${CXXFLAGS} ${LOCALINCLUDES} -I./include ${ROOTINCLUDES}  -o $@  -c $<
endif


${DEP_PLOTS}: ${PRECOMPOUND_TMP}/%.${DEPSuf}: src/%.${SRCSuf}
ifdef CPPVERBOSE
	${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} -I./include $< | sed 's!$*.o!${PRECOMPOUND_TMP}/& $@!' >$@
else
	@echo "Making dependency for file $<..."
	@${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} -I./include $< | sed 's!$*.o!${PRECOMPOUND_TMP}/& $@!' >$@
endif



-include ${DEP_PLOTS}

#-----------------------------------------------------------------------------------

# ${EXE_MAKEPLOTS}: ${LIB_PRECO} ${LIB_EXPDATA} ${LIB_PLOTS} ${OBJ_MAKEPLOTS}
# ifdef CPPVERBOSE
# 	$(LD) $(LDFLAGS) $(ROOTGLIBS) $(LOCALLIBSDIR) ${LOCALLIBGZ} ${LOCALLIBPRECO} ${LOCALLIBEXPDATA} -lmakeplots ${OBJ_MAKEPLOTS} -o $@
# else
# 	@echo "Linking $@..."
# 	@$(LD) $(LDFLAGS) $(ROOTGLIBS) $(LOCALLIBSDIR) ${LOCALLIBGZ} ${LOCALLIBPRECO} ${LOCALLIBEXPDATA} -lmakeplots ${OBJ_MAKEPLOTS} -o $@
# endif


# ${OBJ_MAKEPLOTS}: ${SRC_MAKEPLOTS}
# ifdef CPPVERBOSE
# 	${CXX} ${CXXFLAGS} ${LOCALINCLUDES} -I./include ${ROOTINCLUDES}  -o $@  -c $<
# else
# 	@echo "Compiling $<..."
# 	@${CXX} ${CXXFLAGS} ${LOCALINCLUDES} -I./include ${ROOTINCLUDES} -o $@  -c $<
# endif

# ${DEP_MAKEPLOTS}: ${SRC_MAKEPLOTS}
# ifdef CPPVERBOSE
# 	${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} -I./include $< | sed 's!$*.o!${PRECOMPOUND_TMP}/& $@!' >$@
# else
# 	@echo "Making dependency for file $<..."
# 	@${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} -I./include $< | sed 's!$*.o!${PRECOMPOUND_TMP}/& $@!' >$@
# endif

# -include ${DEP_MAKEPLOTS}



# ------------------------------------------------------------------------------

${EXE_CONVERSION}: ${LIB_PRECO} ${OBJ_CONVERSION}
ifdef CPPVERBOSE
	$(LD) $(LDFLAGS) $(ROOTLIBS) $(LOCALLIBSDIR) ${LOCALLIBGZ} ${LOCALLIBPRECO} ${OBJ_CONVERSION} -o $@
else
	@echo "Linking $@..."
	@$(LD) $(LDFLAGS) $(ROOTLIBS) $(LOCALLIBSDIR) ${LOCALLIBGZ} ${LOCALLIBPRECO} ${OBJ_CONVERSION} -o $@
endif

${OBJ_CONVERSION}: ${SRC_CONVERSION}
ifdef CPPVERBOSE
	${CXX} ${CXXFLAGS} ${LOCALINCLUDES} ${ROOTINCLUDES}  -o $@  -c $<
else
	@echo "Compiling $<..."
	@${CXX} ${CXXFLAGS} ${LOCALINCLUDES} ${ROOTINCLUDES}  -o $@  -c $<
endif

${DEP_CONVERSION}: ${SRC_CONVERSION}
ifdef CPPVERBOSE
	${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} $< | sed 's!$*.o!${PRECOMPOUND_TMP}/& $@!' >$@
else
	@echo "Making dependency for file $<..."
	@${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} $< | sed 's!$*.o!${PRECOMPOUND_TMP}/& $@!' >$@
endif

-include ${DEP_CONVERSION}

#-----------------------------------------------------------------------------------


# ${EXE_SHOWDATA}: ${LIB_EXPDATA} ${OBJ_SHOWDATA}
# ifdef CPPVERBOSE
# 	$(LD)  $(LDFLAGS) $(ROOTLIBS) $(LOCALLIBSDIR) $(LOCALLIBEXPDATA) ${OBJ_SHOWDATA} -o $@
# else
# 	@echo "Linking $@..."
# 	@$(LD)  $(LDFLAGS) $(ROOTLIBS) $(LOCALLIBSDIR) $(LOCALLIBEXPDATA) ${OBJ_SHOWDATA} -o $@
# endif


# ${OBJ_SHOWDATA}: ${SRC_SHOWDATA}
# ifdef CPPVERBOSE
# 	${CXX} ${CXXFLAGS} ${LOCALINCLUDES} ${ROOTINCLUDES} -o $@  -c $<
# else
# 	@echo "Compiling $<..."
# 	@${CXX} ${CXXFLAGS} ${LOCALINCLUDES} ${ROOTINCLUDES} -o $@  -c $<
# endif

# ${DEP_SHOWDATA}: ${SRC_SHOWDATA}
# ifdef CPPVERBOSE
# 	${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} $< | sed 's!$*.o!tmp/& $@!' >$@
# else
# 	@echo "Making dependency for file $<..."
# 	@${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} $< | sed 's!$*.o!tmp/& $@!' >$@
# endif

# -include ${DEP_SHOWDATA}


#-----------------------------------------------------------------------------------


${EXE_GSHOWDATA}: ${LIB_EXPDATA} ${LIB_PLOTS} ${OBJ_GSHOWDATA}
ifdef CPPVERBOSE
	$(LD)  $(LDFLAGS) $(ROOTGLIBS) $(LOCALLIBSDIR) $(LOCALLIBEXPDATA) ${OBJ_GSHOWDATA} \
	-lmakeplots -lexpdata -lPreco -o $@
else
	@echo "Linking $@..."
	@$(LD)  $(LDFLAGS) $(ROOTGLIBS) $(LOCALLIBSDIR) $(LOCALLIBEXPDATA) ${OBJ_GSHOWDATA} \
	-lmakeplots -lexpdata -lPreco -o $@
endif


${OBJ_GSHOWDATA}: ${SRC_GSHOWDATA}
ifdef CPPVERBOSE
	${CXX} ${CXXFLAGS} ${LOCALINCLUDES} ${ROOTINCLUDES} -I./include -o $@  -c $<
else
	@echo "Compiling $<..."
	@${CXX} ${CXXFLAGS} ${LOCALINCLUDES} ${ROOTINCLUDES} -I./include -o $@  -c $<
endif

${DEP_GSHOWDATA}: ${SRC_GSHOWDATA}
ifdef CPPVERBOSE
	${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} -I./include $< | sed 's!$*.o!tmp/& $@!' >$@
else
	@echo "Making dependency for file $<..."
	@${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} -I./include $< | sed 's!$*.o!tmp/& $@!' >$@
endif

-include ${DEP_GSHOWDATA}

#-----------------------------------------------------------------------------------

${EXE_GMAKEPLOTS}: ${LIB_PRECO} ${LIB_EXPDATA} ${LIB_PLOTS} ${OBJ_GMAKEPLOTS}
ifdef CPPVERBOSE
	$(LD) $(LDFLAGS) $(ROOTGLIBS) $(LOCALLIBSDIR) ${LOCALLIBGZ} ${LOCALLIBPRECO} ${LOCALLIBEXPDATA} -lmakeplots ${OBJ_GMAKEPLOTS} -o $@
else
	@echo "Linking $@..."
	@$(LD) $(LDFLAGS) $(ROOTGLIBS) $(LOCALLIBSDIR) ${LOCALLIBGZ} ${LOCALLIBPRECO} ${LOCALLIBEXPDATA} -lmakeplots ${OBJ_GMAKEPLOTS} -o $@
endif


${OBJ_GMAKEPLOTS}: ${SRC_GMAKEPLOTS}
ifdef CPPVERBOSE
	${CXX} ${CXXFLAGS} ${LOCALINCLUDES} -I./include ${ROOTINCLUDES}  -o $@  -c $<
else
	@echo "Compiling $<..."
	@${CXX} ${CXXFLAGS} ${LOCALINCLUDES} -I./include ${ROOTINCLUDES} -o $@  -c $<
endif

${DEP_GMAKEPLOTS}: ${SRC_GMAKEPLOTS}
ifdef CPPVERBOSE
	${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} -I./include $< | sed 's!$*.o!${PRECOMPOUND_TMP}/& $@!' >$@
else
	@echo "Making dependency for file $<..."
	@${CXX} -MM ${LOCALINCLUDES} ${ROOTINCLUDES} -I./include $< | sed 's!$*.o!${PRECOMPOUND_TMP}/& $@!' >$@
endif

-include ${DEP_GMAKEPLOTS}



#-----------------------------------------------------------------------------------

clear:
	@rm -f ./*~
	@rm -f include/*~
	@rm -f src/*~
	@rm -f ./core
	@rm -f ./core.*
	@rm -f ./bin/core
	@rm -f ./bin/core.*

#clean: converclean makeplotsclean plotsclean showdataclean gshowdataclean gmakeplotsclean

clean: converclean  plotsclean  gshowdataclean gmakeplotsclean


converclean: 
	@rm -f ${OBJ_CONVERSION}
	@rm -f ${DEP_CONVERSION}
	@rm -f ${EXE_CONVERSION}


makeplotsclean:
	@rm -f ${OBJ_MAKEPLOTS}
	@rm -f ${DEP_MAKEPLOTS}
	@rm -f ${EXE_MAKEPLOTS}

plotsclean:
	@rm -f ${OBJ_PLOTS}
	@rm -f ${DEP_PLOTS}
	@rm -f ${LIB_PLOTS}

showdataclean:
	@rm -f ${OBJ_SHOWDATA}
	@rm -f ${DEP_SHOWDATA}
	@rm -f ${EXE_SHOWDATA}

gshowdataclean:
	@rm -f ${OBJ_GSHOWDATA}
	@rm -f ${DEP_GSHOWDATA}
	@rm -f ${EXE_GSHOWDATA}

gmakeplotsclean:
	@rm -f ${OBJ_GMAKEPLOTS}
	@rm -f ${DEP_GMAKEPLOTS}
	@rm -f ${EXE_GMAKEPLOTS}


print:
	@echo "----------------------------------------------------------"
	@echo "ARCH     = ${ARCH}"
	@echo "CXX      = ${CXX}"
	@echo "CXXFLAGS = ${CXXFLAGS}"
	@echo "LD       = ${LD}"
	@echo "LDFLAGS  = ${LDFALGS}"
	@echo "SOFLAGS  = ${SOFLAGS}"
	@echo "----------------------------------------------------------"
	@echo "ROOTCFLAGS    = ${ROOTCFLAGS}"
	@echo "ROOTLIBS      = ${ROOTLIBS}"
	@echo "ROOTGLIBS     = ${ROOTGLIBS}"
	@echo "ROOTLIBS      = ${ROOTLIBS}"
	@echo "ROOTGLIBS     = ${ROOTGLIBS}"
	@echo "ROOTINCLUDES  = ${ROOTINCLUDES}"
	@echo "LOCALINCLUDES = ${LOCALINCLUDES}"
	@echo "----------------------------------------------------------"
	@echo "SRCSuf = ${SRCSuf}"
	@echo "OBJSuf = ${OBJSuf}"
	@echo "DEPSuf = ${DEPSuf}"
	@echo "LIBSuf = ${LIBSuf}"
	@echo "----------------------------------------------------------"
	@echo "SRC_PLOTS = ${SRC_PLOTS}"
	@echo "OBJ_PLOTS = ${OBJ_PLOTS}"
	@echo "DEP_PLOTS = ${DEP_PLOTS}"
	@echo "LIB_PLOTS = ${LIB_PLOTS}"
	@echo "----------------------------------------------------------"
	@echo "SRC_MAKEPLOTS = ${SRC_MAKEPLOTS}"
	@echo "OBJ_MAKEPLOTS = ${OBJ_MAKEPLOTS}"
	@echo "DEP_MAKEPLOTS = ${DEP_MAKEPLOTS}"
	@echo "EXE_MAKEPLOTS = ${EXE_MAKEPLOTS}"
	@echo "----------------------------------------------------------"
	@echo "SRC_SHOWDATA = ${SRC_SHOWDATA}"
	@echo "OBJ_SHOWDATA = ${OBJ_SHOWDATA}"
	@echo "DEP_SHOWDATA = ${DEP_SHOWDATA}"
	@echo "EXE_SHOWDATA = ${EXE_SHOWDATA}"
	@echo "----------------------------------------------------------"
	@echo "SRC_GSHOWDATA = ${SRC_GSHOWDATA}"
	@echo "OBJ_GSHOWDATA = ${OBJ_GSHOWDATA}"
	@echo "DEP_GSHOWDATA = ${DEP_GSHOWDATA}"
	@echo "EXE_GSHOWDATA = ${EXE_GSHOWDATA}"
	@echo "----------------------------------------------------------"
	@echo "SRC_GMAKEPLOTS = ${SRC_GMAKEPLOTS}"
	@echo "OBJ_GMAKEPLOTS = ${OBJ_GMAKEPLOTS}"
	@echo "DEP_GMAKEPLOTS = ${DEP_GMAKEPLOTS}"
	@echo "EXE_GMAKEPLOTS = ${EXE_GMAKEPLOTS}"
	@echo "----------------------------------------------------------"
