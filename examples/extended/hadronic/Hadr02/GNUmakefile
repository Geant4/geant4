# ----------------------------------------------------------------
# GNUmakefile for examples module with external Fortran generators
# ----------------------------------------------------------------

name := Hadr02
G4TARGET := $(name)
G4EXLIB := true

URQMDDIR= urqmd1_3
HIJINGDIR= hijing1_383
CRMCDIR= crmc_interface

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

.PHONY: all
all: urqmd hijing crmc lib bin


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
              -lgmp -lmpfr \
	      -L${CERNROOT}/lib -lmathlib -lkernlib -lpacklib -lpawlib
endif

crmc:
ifdef G4_USE_CRMC
	(cp ${CRMCDIR}/Build/lib/*.so ${G4TMPDIR});
	(cp ${CRMCDIR}/Build/crmc.param .);
endif

ifdef G4_USE_CRMC
  CPPFLAGS += -DG4_USE_CRMC -I${CRMCDIR}/Build/src -I${CRMCDIR}/src
  EXTRALIBS = -lgfortran -lgmp -lmpfr \
              -L${CRMCDIR}/Build/lib/ -lCrmc 
endif

include $(G4INSTALL)/config/binmake.gmk

dclean:
	@rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/G4DPMJET2_5Model.o

histclean:
	@rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/Histo.o

visclean:
	@rm -f g4*.prim g4*.eps g4*.wrl
	@rm -f .DAWN_*
#
