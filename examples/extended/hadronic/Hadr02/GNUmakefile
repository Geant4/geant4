# $Id: GNUmakefile 77519 2013-11-25 10:54:57Z gcosmo $
# GEANT4 tag $Name: 
# --------------------------------------------------------------
# GNUmakefile for examples module with DPMJET interface
# --------------------------------------------------------------

name := Hadr02
G4TARGET := $(name)
G4EXLIB := true

DPMJDIR = dpmjet2_5
URQMDDIR= urqmd1_3
HIJINGDIR= hijing1_383

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

.PHONY: all
all: dpm urqmd hijing lib bin

dpm:
ifdef G4_USE_DPMJET
	(cd ${DPMJDIR} && $(MAKE));
	(mv ${DPMJDIR}/*.o ${G4TMPDIR});
endif

ifdef G4_USE_DPMJET
  CPPFLAGS += -DG4_USE_DPMJET \
              -DG4DPMJET25GDSHDEBUG

  EXTRALIBS = -lgfortran -lgmp -lmpfr \
              -L${CERNLIB}/lib -lmathlib -lkernlib -lpacklib 
endif


urqmd:
ifdef G4_USE_URQMD
	(cd ${URQMDDIR} && \
         cp GNUmakefile urqmd-1.3cr && cp *.f urqmd-1.3cr &&\
	 cd urqmd-1.3cr && ${MAKE} TYPE="G4INTERFACE");
	( mv ${URQMDDIR}/urqmd-1.3cr/obj_G4INTERFACE/*.o ${G4TMPDIR} );
endif

ifdef G4_USE_URQMD
  CPPFLAGS += -DG4_USE_URQMD

  EXTRALIBS = -lgfortran -lgmp -lmpfr \
              -L${CERNLIB}/lib -lmathlib -lkernlib -lpacklib 
endif

hijing:
ifdef G4_USE_HIJING
	(cd ${HIJINGDIR} && ${MAKE});
	(mv ${HIJINGDIR}/obj_Linux/*.o ${G4TMPDIR});
endif

ifdef G4_USE_HIJING
  CPPFLAGS += -DG4_USE_HIJING
  EXTRALIBS = -lgfortran 
#-lgmp -lmpfr \
#		-L${CERNROOT}/cernlib -lmathlib -lkernlib -lpacklib -lpawlib
endif


include $(G4INSTALL)/config/binmake.gmk

histclean:
	@rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/Histo.o

visclean:
	@rm -f g4*.prim g4*.eps g4*.wrl
	@rm -f .DAWN_*

#
