// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventorX.hh,v 1.4 1999-12-15 14:54:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  27th March 1996
// OpenInventor graphics system factory.

#if defined (G4VIS_BUILD_OIX_DRIVER) || defined (G4VIS_USE_OIX)

#ifndef G4OPENINVENTORX_HH
#define G4OPENINVENTORX_HH

#include "G4OpenInventor.hh"

class G4OpenInventorX: public G4OpenInventor {
public:
  G4OpenInventorX ();
  virtual ~G4OpenInventorX ();
};

#endif

#endif
