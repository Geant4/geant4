# $Id: GNUmakefile,v 1.2 2004-05-25 12:39:07 guatelli Exp $
# --------------------------------------------------------------

name := Brachy
G4TARGET := $(name)
G4EXLIB := true
G4LIB_BUILD_SHARED = 1

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all

ifdef G4DIANE_USE
all: lib bin
	@echo "Linking DIANE libraries..."
	@g++ -shared -o${G4WORKDIR}/tmp/${G4SYSTEM}/${name}/libDIANE_${name}.so -L${G4WORKDIR}/tmp/${G4SYSTEM}/${name} -l${name} ${LDFLAGS} ${LDLIBS} 
else
all: lib bin
endif

include $(G4INSTALL)/config/binmake.gmk
ifdef G4ANALYSIS_USE  
 CPPFLAGS += `aida-config --include`
 LDFLAGS  += `aida-config --lib`
 LOADLIBS += `aida-config --lib` 
endif


