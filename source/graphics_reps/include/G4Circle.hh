// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Circle.hh,v 1.3 1999-11-17 07:39:20 stanaka Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  17/11/96.

// Class Description:
// G4Circle is a kind of 3D-position marker. 
// Its shape is 2-dimensional circle.
// It inherits G4VMarker.  See G4VMarker.hh for more details.
// Class Description - End:

#ifndef G4CIRCLE_HH
#define G4CIRCLE_HH

#include "G4VMarker.hh"

class G4Circle: public G4VMarker {

public: // With description

  G4Circle ();
  G4Circle (const G4Point3D& pos);
  G4Circle (const G4VMarker& marker);
  virtual ~G4Circle ();

  //////////////////////////////////////////////////////
  // Assignment...
  virtual G4Visible&  operator = (const G4Visible& right);
  virtual G4VVisPrim& operator = (const G4VVisPrim& right);
  virtual G4VMarker&  operator = (const G4VMarker& right);
  virtual G4Circle &  operator = (const G4Circle& right);

 };

#include "G4Circle.icc"

#endif
