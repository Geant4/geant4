#  $Id: GNUmakefile,v 1.5 2003/12/08 17:53:25 gcosmo Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-06-00 $
# --------------------------------------------------------------



TIARA_SRC_DIRS =  tiaraPhysicsPackaging tiaraPhysicsLists tiara 
TIARA_SRC_SWIG_DIRS = CLHEPWrapper G4KernelWrapper TiaraWrapper

include $(G4INSTALL)/config/architecture.gmk

.PHONY: Tiara tiaraClean swigClean all
all: Tiara

Tiara:
	@for dir in $(TIARA_SRC_DIRS); do ( cd source/$$dir && $(MAKE)); done;:
	@for dir in $(TIARA_SRC_SWIG_DIRS); do ( cd source/$$dir && $(MAKE)); done;:

tiaraClean:
	@for dir in $(TIARA_SRC_DIRS); do ( cd source/$$dir && $(MAKE) tiaraClean); done;:
	@for dir in $(TIARA_SRC_SWIG_DIRS); do ( cd source/$$dir && $(MAKE) tiaraClean); done;:

swigClean:
	@for dir in $(TIARA_SRC_SWIG_DIRS); do ( cd source/$$dir && $(MAKE) swigClean); done;:

