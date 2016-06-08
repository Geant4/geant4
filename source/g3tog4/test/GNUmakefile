# $Id: GNUmakefile,v 1.3 1999/05/01 21:42:04 lockman Exp $
# ----------------------------------------------------------------

ifndef G4TARGET
G4TARGET = clGeometry
endif

name = g3tog4_example
 
## following flag indicates that we will build a temporary library named
## lib$(name).a containing .o files from source files in src/

G4EXLIB := true

## build a separate CFRONT template repository for the code in this example

G4EXEC_BUILD = true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4UI_USE.gmk
include $(G4INSTALL)/config/G4VIS_USE.gmk
include $(G4INSTALL)/config/interactivity.gmk

## Override some variables for binmake.gmk.

INCFLAGS  := \
        -I$(G4BASE)/global/management/include \
        -I$(G4BASE)/global/HEPRandom/include \
        -I$(G4BASE)/global/HEPGeometry/include \
        -I$(G4BASE)/intercoms/include \
        -I$(G4BASE)/geometry/management/include \
        -I$(G4BASE)/geometry/volumes/include \
        -I$(G4BASE)/geometry/magneticfield/include \
        -I$(G4BASE)/geometry/solids/CSG/include \
        -I$(G4BASE)/geometry/solids/BREPS/include \
        -I$(G4BASE)/geometry/solids/STEP/include \
        -I$(G4BASE)/tracking/include \
        -I$(G4BASE)/track/include \
        -I$(G4BASE)/event/include \
        -I$(G4BASE)/materials/include \
        -I$(G4BASE)/particles/management/include \
        -I$(G4BASE)/particles/bosons/include \
        -I$(G4BASE)/particles/leptons/include \
        -I$(G4BASE)/particles/hadrons/mesons/include \
        -I$(G4BASE)/particles/hadrons/barions/include \
        -I$(G4BASE)/particles/hadrons/ions/include \
        -I$(G4BASE)/processes/management/include \
        -I$(G4BASE)/processes/decay/include \
        -I$(G4BASE)/processes/electromagnetic/standard/include \
        -I$(G4BASE)/processes/electromagnetic/muons/include \
        -I$(G4BASE)/processes/electromagnetic/utils/include \
        -I$(G4BASE)/processes/electromagnetic/xrays/include \
        -I$(G4BASE)/graphics_reps/include \
        -I$(G4BASE)/digits+hits/detector/include \
        -I$(G4BASE)/digits+hits/digits/include \
        -I$(G4BASE)/digits+hits/hits/include \
        -I$(G4BASE)/run/include \
        -I$(G4BASE)/g3tog4/include \
        -I$(G4BASE)/g3tog4/test/include

### LDLIBS must be specified since -lG3toG4 is not part of the standard 
### libraries

LDLIBS   := -l$(name) $(VISLIBS) $(UILIBS) \
        -lG4run \
        -lG4event \
        -lG3toG4 \
        -lG4tracking \
        -lG4processes \
        -lG4particles \
        -lG4track \
        -lG4materials \
        -lG4digits+hits \
        -lG4geometry \
        -lG4graphics_reps \
        -lG4intercoms \
        -lG4global \
	$(INTYLIBS)

.PHONY: first lib bin test1

first: $(G4INSTALL)/lib/$(G4SYSTEM)/libG3toG4.a lib bin

$(G4INSTALL)/lib/$(G4SYSTEM)/libG3toG4.a : $(G4BASE)/g3tog4/src/*.cc
	cd $(G4BASE)/g3tog4; $(MAKE) lib

##############################################################################
include $(G4INSTALL)/config/binmake.gmk

test1: 
	echo LDLIBS: $(LDLIBS)










