# $Id: GNUmakefile,v 1.2 2000-11-15 20:27:37 flongo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := GammaRayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

#
# AIDA_EXAMPLES GNUmakefile
#
# Object files are stored in $(OS)
# 
# Deleting object files : gmake clean
#

# in order to set up the environment properly, define
# the following variables (examples are shown for  CERN)
# 
# setenv LHCXXTOP /afs/cern.ch/sw/lhcxx
# setenv PLATF redhat61
# 
# and execute the following command:
# source $LHCXXTOP/share/LHCXX/3.0.0/install/sharedstart.csh
#

# use the afs installed version of LHCXX

#ifeq "$(OS)" "Solaris"
# include $(LHCXX_REL_DIR)/etc/ospace.mk
# C_INCLUDES += -I$(OSPACE_DIR)/ospace/stl 
# C_FLAGS    += -features=rtti 
#endif
#ifeq "$(OS)" "Solaris-7"
# C_FLAGS += -DAIDA_USE_STD
#endif


ifdef G4HIS_USE_AIDA
  export USE_HEPODBMS_CLIENT=1 
  include $(LHCXX_REL_DIR)/etc/lhcxx.mk
  CPPFLAGS += -I. -I$(LHCXX_REL_DIR)/include 
# for now ... qplotter-includes need to be fixed ...
  CPPFLAGS += -I$(LHCXXTOP)/specific/$(PLATF)/Qplotter/1.0.0.0 \
	      -I$(LHCXXTOP)/specific/$(PLATF)/Qt/2.1.0/include
  LOADLIBS += $(LHCXX_ALL_LIBS) \
            -L/afs/cern.ch/sw/geant4/dev/Mesa/Linux-g++/lib/ -lGLU
endif







