# $Id: GNUmakefile 77567 2013-11-26 09:40:31Z gcosmo $
# ----------------------------------------------------------------------
# GNUmakefile for physicslists library, Gunter Folger 25-October 2006
# ----------------------------------------------------------------------

name := G4physlist_ctor

SUBDIRS  = decay
SUBDIRS += electromagnetic
SUBDIRS += factory
SUBDIRS += gamma_lepto_nuclear
SUBDIRS += hadron_elastic
SUBDIRS += hadron_inelastic
SUBDIRS += ions
SUBDIRS += limiters
SUBDIRS += stopping
 


SUBLIBS  = G4phys_ctor_decay
SUBLIBS += G4phys_ctor_em
SUBLIBS += G4phys_ctor_factory
SUBLIBS += G4phys_ctor_glnuclear
SUBLIBS += G4phys_ctor_helastic
SUBLIBS += G4phys_ctor_hinelastic
SUBLIBS += G4phys_ctor_ions
SUBLIBS += G4phys_ctor_limiters
SUBLIBS += G4phys_ctor_stopping

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/globlib.gmk
