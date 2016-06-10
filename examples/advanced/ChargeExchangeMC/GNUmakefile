name := ChargeExchangeMC
G4TARGET := $(name)
G4EXLIB := true

CPPFLAGS += -DCEXMC_PROG_NAME=\"$(name)\"

# if CEXMC_USE_PERSISTENCY is 'yes' then run and events data can be read and
# written; requires boost::serialize headers and library
CEXMC_USE_PERSISTENCY := no
# if CEXMC_USE_CUSTOM_FILTER is 'yes' then Custom filter can be used for
# existing events data; requires boost::spirit 2.x headers. Notice: if
# CEXMC_USE_PERSISTENCY is not 'yes' then Custom Filter will not be used anyway
CEXMC_USE_CUSTOM_FILTER := no
# if CEXMC_DEBUG_CUSTOM_FILTER is 'yes' then AST trees will be printed out
CEXMC_DEBUG_CUSTOM_FILTER := no
# if CEXMC_USE_HISTOGRAMING is 'yes' then ROOT histograming framework will be
# compiled. Notice: if ROOT CERN is not installed in your system then the
# histograming module won't compile anyway
CEXMC_USE_HISTOGRAMING := yes
# if CEXMC_USE_QGSP_BERT is 'yes' then QGSP_BERT will be used as basic physics,
# otherwise - FTFP_BERT or QGSP_BIC_EMY
CEXMC_USE_QGSP_BERT := no
# if CEXMC_USE_QGSP_BIC_EMY is 'yes' then QGSP_BIC_EMY will be used as basic
# physics, otherwise - FTFP_BERT or QGSP_BERT
CEXMC_USE_QGSP_BIC_EMY := no
# if CEXMC_USE_GENBOD is 'yes' then original FORTRAN routine GENBOD() will be
# used as phase space generator
CEXMC_USE_GENBOD := no
# if CEXMC_DEBUG_TP is 'yes' then additional info will be printed on track
# points data
CEXMC_DEBUG_TP := no


ifndef G4INSTALL
  G4INSTALL = ../../..
endif

ifeq ($(CEXMC_USE_GENBOD),yes)
  CPPFLAGS += -DCEXMC_USE_GENBOD
  EXTRALIBS += `cernlib geant321 phtools packlib kernlib`
  GCC_VERSION := $(shell gcc --version | head -1 | awk '{ printf $$3 }' | \
                         awk -F"." '{ printf $$1 }')
  ifdef CEXMC_FORTRAN_LIB
    EXTRALIBS += $(CEXMC_FORTRAN_LIB)
  else
    # try to setup fortran lib automatically
    # WARNING: the following is not robust check because cernlib can be built
    # against libg2c even when using gcc-4 series
    # Please define CEXMC_FORTRAN_LIB if the check fails
    ifeq ($(GCC_VERSION),3)
      EXTRALIBS += -lg2c
    else
      EXTRALIBS += -lgfortran
    endif
  endif
endif

ifdef BOOST_INCLUDE_PATH
  CPPFLAGS += -I$(BOOST_INCLUDE_PATH)
endif

ifdef BOOST_LIBRARY_PATH
  EXTRALIBS += -L$(BOOST_LIBRARY_PATH)
endif

ifeq ($(CEXMC_USE_PERSISTENCY),yes)
  EXTRALIBS += -lboost_serialization
  CPPFLAGS += -DCEXMC_USE_PERSISTENCY
  ifeq ($(CEXMC_USE_CUSTOM_FILTER),yes)
    CPPFLAGS += -DCEXMC_USE_CUSTOM_FILTER
    ifeq ($(CEXMC_DEBUG_CUSTOM_FILTER),yes)
      CPPFLAGS += -DCEXMC_DEBUG_CF
    endif
  endif
endif

ifeq ($(CEXMC_USE_HISTOGRAMING),yes)
  # try to determine if ROOT will be used automatically
  USE_ROOT := $(shell which root-config 2>/dev/null)
  ifneq ($(USE_ROOT),)
    CPPFLAGS += -I`root-config --incdir`
    EXTRALIBS += `root-config --libs`
    CPPFLAGS += -DCEXMC_USE_ROOT
    # try to determine if ROOT-Qt binding will be used automatically
    USE_ROOTQT := $(shell root-config --features | grep qt)
    ifneq ($(USE_ROOTQT),)
      EXTRALIBS += -lGQt
      CPPFLAGS += -DCEXMC_USE_ROOTQT
    endif
  endif
endif

ifeq ($(CEXMC_USE_QGSP_BERT),yes)
  CPPFLAGS += -DCEXMC_USE_QGSP_BERT
else
  ifeq ($(CEXMC_USE_QGSP_BIC_EMY),yes)
    CPPFLAGS += -DCEXMC_USE_QGSP_BIC_EMY
  endif
endif

ifeq ($(CEXMC_DEBUG_TP),yes)
  CPPFLAGS += -DCEXMC_DEBUG_TP
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

