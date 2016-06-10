# $Id: GNUmakefile 67969 2013-03-13 09:44:42Z gcosmo $
# ----------------------------------------------------------------------
# GNUmakefile for physicslists library, Gunter Folger 25-October 2006
# ----------------------------------------------------------------------

name := G4physicslists

SUBDIRS  = lists 
SUBDIRS += builders
SUBDIRS += constructors/decay
SUBDIRS += constructors/electromagnetic
SUBDIRS += constructors/factory
SUBDIRS += constructors/gamma_lepto_nuclear
SUBDIRS += constructors/hadron_elastic
SUBDIRS += constructors/hadron_inelastic
SUBDIRS += constructors/ions
SUBDIRS += constructors/limiters
SUBDIRS += constructors/stopping
SUBDIRS += util

SUBLIBS  = G4phys_lists 
SUBLIBS += G4phys_builders 
SUBLIBS += G4phys_ctor_decay
SUBLIBS += G4phys_ctor_em
SUBLIBS += G4phys_ctor_factory
SUBLIBS += G4phys_ctor_glnuclear
SUBLIBS += G4phys_ctor_helastic
SUBLIBS += G4phys_ctor_hinelastic
SUBLIBS += G4phys_ctor_ions
SUBLIBS += G4phys_ctor_limiters
SUBLIBS += G4phys_ctor_stopping
SUBLIBS += G4phys_util

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk
