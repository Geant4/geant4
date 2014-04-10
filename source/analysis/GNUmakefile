# $Id$
# --------------------------------------------------------------
#

name := G4analysis

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS = libG4intercoms.lib libG4global.lib
ifdef G4LIB_BUILD_ZLIB
  GLOBLIBS += libG4zlib.lib
endif

include $(G4INSTALL)/config/architecture.gmk

SUBDIRS  = management g4tools hntools csv root xml
SUBLIBS  = G4analysismng G4tools G4hntools G4csv G4root G4xml

include $(G4INSTALL)/config/globlib.gmk

includes:: 
	@for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) $@ ); done
	cd g4tools; make includes
