#  $Id: GNUmakefile,v 1.3 2003-12-02 13:30:50 gcosmo Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: not supported by cvs2svn $
# --------------------------------------------------------------



TIARA_SRC_DIRS =  tiaraPhysicsPackaging tiaraPhysicsLists tiara 
TIARA_SRC_SWIG_DIRS = CLHEPWrapper G4KernelWrapper TiaraWrapper

include $(G4INSTALL)/config/architecture.gmk

.PHONY: Tiara tiaraClean swigClean all
all: lib bin

Tiara:
	@for dir in $(TIARA_SRC_DIRS); do ( cd source/$$dir && $(MAKE)); done;:
	@for dir in $(TIARA_SRC_SWIG_DIRS); do ( cd source/$$dir && $(MAKE)); done;:

tiaraClean:
	@for dir in $(TIARA_SRC_DIRS); do ( cd source/$$dir && $(MAKE) tiaraClean); done;:
	@for dir in $(TIARA_SRC_SWIG_DIRS); do ( cd source/$$dir && $(MAKE) tiaraClean); done;:

swigClean:
	@for dir in $(TIARA_SRC_SWIG_DIRS); do ( cd source/$$dir && $(MAKE) swigClean); done;:

include $(G4INSTALL)/config/binmake.gmk

