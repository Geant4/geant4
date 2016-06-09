# $Id: GNUmakefile,v 1.4 2009/03/06 09:51:42 gbarrand Exp $
# --------------------------------------------------------------

name := Brachy
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all

ifdef G4DIANE_USE
G4LIB_BUILD_SHARED = 1
all: lib bin
	@echo "Linking DIANE libraries..."
	@g++ -shared -o${G4WORKDIR}/tmp/${G4SYSTEM}/${name}/libDIANE_${name}.so -L${G4WORKDIR}/tmp/${G4SYSTEM}/${name} -l${name} ${LDFLAGS} ${LDLIBS} 
else
all: lib bin
endif

include $(G4INSTALL)/config/binmake.gmk

