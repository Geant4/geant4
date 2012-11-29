//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
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

  G4VSolid* Clone() const;
    // Returns a pointer of a dynamically allocated copy of the solid.

  std::ostream& StreamInfo(std::ostream& os) const;
    // Streams solid contents to output stream.

 public:  // without description

  G4BREPSolidTorus(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  G4BREPSolidTorus(const G4BREPSolidTorus& rhs);
  G4BREPSolidTorus& operator=(const G4BREPSolidTorus& rhs);
    // Copy constructor and assignment operator.

 private:

  void InitializeTorus();

 private:

  struct G4BREPTorusParams
  {
    G4ThreeVector origin;
    G4ThreeVector axis;
    G4ThreeVector direction;
    G4double      MinRadius;
    G4double      MaxRadius;
  } constructorParams;
  
};

#endif
