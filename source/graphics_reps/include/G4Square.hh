// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Square.hh,v 1.4 1999-12-15 14:50:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  17/11/96.

// Class Description:
// G4Square is a kind of 3D-position marker. 
// Its shape is 2-dimensional right square.
// It inherits G4VMarker.  See G4VMarker.hh for more details.
// Class Description - End:

#ifndef G4SQUARE_HH
#define G4SQUARE_HH

#include "G4VMarker.hh"

class G4Square: public G4VMarker {

public: // With description

  G4Square ();
  G4Square (const G4Point3D& pos);
  G4Square (const G4VMarker& marker);
  virtual ~G4Square ();

  //////////////////////////////////////////////////////
  // Assignment...
  virtual G4Visible&  operator = (const G4Visible& right);
  virtual G4VVisPrim& operator = (const G4VVisPrim& right);
  virtual G4VMarker&  operator = (const G4VMarker& right);
  virtual G4Square&   operator = (const G4Square& right);

};

#include "G4Square.icc"

#endif
