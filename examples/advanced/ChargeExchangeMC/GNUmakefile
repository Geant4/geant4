name := ChargeExchangeMC
G4TARGET := $(name)
G4EXLIB := true

# if CEXMC_USE_GENBOD is 'yes' then original FORTRAN routine GENBOD() will be
# used as phase space generator
CEXMC_USE_GENBOD := no

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

ifdef BOOST_INCLUDE_PATH
  CPPFLAGS += -I$(BOOST_INCLUDE_PATH)
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

ifdef BOOST_LIBRARY_PATH
  EXTRALIBS += -L$(BOOST_LIBRARY_PATH)
endif
EXTRALIBS += -lboost_serialization

# try to determine if ROOT will be used automatically
USE_ROOT := $(shell which root-config 2>/dev/null)

ifneq ($(USE_ROOT)),)
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

#CPPFLAGS += -DCEXMC_USE_QGSP_BIC_EMY
#CPPFLAGS += -DCEXMC_USE_CUSTOM_FILTER
#CPPFLAGS += -DCEXMC_DEBUG_CF
CPPFLAGS += -DCEXMC_DEBUG_TP

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
