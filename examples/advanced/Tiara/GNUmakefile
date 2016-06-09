#  $Id: GNUmakefile,v 1.7 2004/12/08 15:37:13 daquinog Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-07-01 $
# --------------------------------------------------------------



TIARA_SRC_DIRS = tiara
TIARA_SRC_SWIG_DIRS = CLHEPWrapper G4KernelWrapper TiaraWrapper

include $(G4INSTALL)/config/architecture.gmk

.PHONY: Tiara tiaraClean swigClean clean_all all
all: Tiara

#ifdef G4ANALYSIS_USE
# CPPFLAGS += `aida-config --include`
# LDFLAGS  += `aida-config --lib`
# LDLIBS += `aida-config --lib`
#endif

Tiara:
	@for dir in $(TIARA_SRC_DIRS); do ( cd source/$$dir && $(MAKE)); done;:
	@for dir in $(TIARA_SRC_SWIG_DIRS); do ( cd source/$$dir && $(MAKE)); done;:

clean_all: swigClean tiaraClean
	@rm -f source/py_modules/*.pyc
	@rm -rf run/simData

tiaraClean:
	@for dir in $(TIARA_SRC_DIRS); do ( cd source/$$dir && $(MAKE) tiaraClean); done;:
	@for dir in $(TIARA_SRC_SWIG_DIRS); do ( cd source/$$dir && $(MAKE) tiaraClean); done;:

swigClean:
	@for dir in $(TIARA_SRC_SWIG_DIRS); do ( cd source/$$dir && $(MAKE) swigClean); done;:
