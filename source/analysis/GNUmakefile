# $Id: GNUmakefile 87171 2014-11-26 09:37:08Z gcosmo $
# --------------------------------------------------------------
#

name := G4analysis

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS = libG4intercoms.lib libG4global.lib
ifdef G4LIB_BUILD_EXPAT
  GLOBLIBS += libG4expat.lib
endif
ifdef G4LIB_BUILD_ZLIB
  GLOBLIBS += libG4zlib.lib
endif

include $(G4INSTALL)/config/architecture.gmk

SUBDIRS  = g4tools management hntools csv root xml
SUBLIBS  = G4tools G4analysismng G4hntools G4csv G4root G4xml

include $(G4INSTALL)/config/globlib.gmk

includes:: 
	@for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) $@ ); done
