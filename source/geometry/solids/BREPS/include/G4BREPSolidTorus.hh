// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidTorus.hh,v 1.4 2000-11-08 14:21:59 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BREPSolidTorus
//
// Class description:
// 
//  Definition of a generic BREP torus:
//
//  G4BREPSolidTorus(const G4String& name,
//                   const G4ThreeVector& origin,
//                   const G4ThreeVector& axis,
//                   const G4ThreeVector& direction,
//                         G4double MinRadius,
//                         G4double MaxRadius)

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BREPSolidTorus
#define __G4BREPSolidTorus

#include "G4BREPSolid.hh"

class G4BREPSolidTorus : public G4BREPSolid
{

 public:  // with description

  G4BREPSolidTorus(const G4String& name,
		   const G4ThreeVector& origin,
		   const G4ThreeVector& axis,
		   const G4ThreeVector& direction,
		   G4double MinRadius,
		   G4double MaxRadius);
    // Constructor defining the torus through its surfaces.

  ~G4BREPSolidTorus();
    // Empty destructor.

 private:

  G4BREPSolidTorus(const G4BREPSolidTorus&);
  G4BREPSolidTorus& operator=(const G4BREPSolidTorus&);
    // Private copy constructor and assignment operator.

};

#endif
