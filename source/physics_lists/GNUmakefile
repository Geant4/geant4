# $Id: GNUmakefile 79957 2014-03-27 14:29:51Z gcosmo $
# ----------------------------------------------------------------------
# GNUmakefile for physicslists library, Gunter Folger 25-October 2006
# ----------------------------------------------------------------------

name := G4physicslists

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS  = libG4run.lib libG4event.lib
GLOBLIBS += libG4tracking.lib libG4processes.lib libG4digits_hits.lib
GLOBLIBS += libG4track.lib libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4graphics_reps.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib
GLOBLIBS += libG4expat.lib

include $(G4INSTALL)/config/architecture.gmk

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
SUBLIBS += G4physlist_util

include $(G4INSTALL)/config/globlib.gmk
