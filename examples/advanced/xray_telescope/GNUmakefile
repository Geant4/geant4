# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.3 2000-11-06 18:02:31 nartallo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := XrayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

G4VIS_USE = 1
G4VIS_USE_DAWNFILE = 1	
G4VIS_USE_DAWN = 1	
G4VIS_USE_VRMLFILE = 1	
G4VIS_USE_VRML = 1

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

###export USE_HEPODBMS_CLIENT=1
###include $(LHCXX_REL_DIR)/etc/lhcxx.mk

#ifeq "$(OS)" "Solaris"
# include $(LHCXX_REL_DIR)/etc/ospace.mk
# C_INCLUDES += -I$(OSPACE_DIR)/ospace/stl 
# C_FLAGS    += -features=rtti 
#endif

#ifeq "$(OS)" "Solaris-7"
# C_FLAGS += -DAIDA_USE_STD
#endif

###CPPFLAGS += -I. -I$(LHCXX_REL_DIR)/include 

# for now ... qplotter-includes need to be fixed ...
###CPPFLAGS += -I$(LHCXXTOP)/specific/$(PLATF)/Qplotter/1.0.0.0 \
	      -I$(LHCXXTOP)/specific/$(PLATF)/Qt/2.1.0/include

###LOADLIBS += $(LHCXX_ALL_LIBS)

