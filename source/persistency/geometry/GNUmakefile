# --------------------------------------------------------------
# GNUmakefile for persistency geometry subcaterogy.
#                                      Youhei Morita, 09/07/98.
# --------------------------------------------------------------

# order of SUBDIRS are important for creating header files.
# don't mess it up!
SUBDIRS = global solids/CSG solids/specific solids/Boolean management 
SUBLIBS = G4pgeomGlobal G4pcsg G4pspecsolids G4pgeomBoolean G4pgeomn 

ifndef G4INSTALL
  G4INSTALL = ../../..
endif
include $(G4INSTALL)/config/globlib.gmk

# --------------------------------------------------------------

glob global:
	@echo Error - gmake global is not supported at this directory.
	exit 1

