# --------------------------------------------------------------
# GNUmakefile for persistency geometry subcaterogy.
#                                      Youhei Morita, 09/07/98.
# --------------------------------------------------------------

# order of SUBDIRS are important for creating header files.
# don't mess it up!
SUBDIRS = management solids/CSG volumes
SUBLIBS = G4pgeomn G4pcsg G4pvol

ifndef G4INSTALL
  G4INSTALL = ../../..
endif
include $(G4INSTALL)/config/globlib.gmk

# --------------------------------------------------------------

glob global:
	@echo Error - gmake global is not supported at this directory.
	exit 1

