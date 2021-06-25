# --------------------------------------------------------------
#

name := G4analysis

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS = libG4intercoms.lib libG4global.lib libG4tools.lib
ifdef G4LIB_BUILD_EXPAT
  GLOBLIBS += libG4expat.lib
endif
ifdef G4LIB_BUILD_ZLIB
  GLOBLIBS += libG4zlib.lib
endif

include $(G4INSTALL)/config/architecture.gmk

SUBDIRS  = management hntools csv root xml factory accumulables
SUBLIBS  = G4analysismng G4hntools G4csv G4root G4xml G4analysisfac G4accumulables

include $(G4INSTALL)/config/globlib.gmk

includes::
	@for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) $@ ); done
