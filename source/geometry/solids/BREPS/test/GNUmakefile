# $Id: GNUmakefile,v 1.6 2000-02-28 15:38:27 gcosmo Exp $
# ------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

# Override some variables from GNUmakefile.architecture.
INCFLAGS := \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
	    -I$(G4BASE)/geometry/solids/STEP/include \
	    -I$(G4BASE)/geometry/solids/Boolean/include \
	    -I$(G4BASE)/geometry/solids/STEPinterface/include \
	    -I$(G4BASE)/geometry/solids/BREPS/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include

LDLIBS   := \
	    -lG4stepinterface \
	    -lG4brep \
	    -lG4geomBoolean \
	    -lG4step \
	    -lG4csg \
	    -lG4geometrymng \
            -lG4graphics_reps \
            -lG4globman

include $(G4INSTALL)/config/binmake.gmk

.PHONY: clean_all

clean_all:
	@for i in CurveTest.cc G4*.cc; do \
	(target=`basename $$i .cc`; \
	echo  "Removing $$target ... "; \
	$(MAKE) clean G4SKIP_DEPEND=1 G4TARGET=$$target); done
	@echo  "Removing STEPTest ... "
	@$(MAKE) clean G4SKIP_DEPEND=1 G4TARGET=STEPTest
	@echo  "Wiping output files ... "
	@$(RM) -rf *.newout*
