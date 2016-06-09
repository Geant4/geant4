#  $Id: GNUmakefile,v 1.2 2003/06/16 17:06:43 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-05-02 $
# --------------------------------------------------------------



TIARA_SRC_DIRS =  tiaraPhysicsPackaging tiaraPhysicsLists tiara 
TIARA_SRC_SWIG_DIRS = CLHEPWrapper G4KernelWrapper TiaraWrapper

include $(G4INSTALL)/config/architecture.gmk

.PHONY: Tiara tiaraClean swigClean

Tiara:
	@for dir in $(TIARA_SRC_DIRS); do ( cd source/$$dir && $(MAKE)); done;:
	@for dir in $(TIARA_SRC_SWIG_DIRS); do ( cd source/$$dir && $(MAKE)); done;:

tiaraClean:
	@for dir in $(TIARA_SRC_DIRS); do ( cd source/$$dir && $(MAKE) tiaraClean); done;:
	@for dir in $(TIARA_SRC_SWIG_DIRS); do ( cd source/$$dir && $(MAKE) tiaraClean); done;:

swigClean:
	@for dir in $(TIARA_SRC_SWIG_DIRS); do ( cd source/$$dir && $(MAKE) swigClean); done;:
