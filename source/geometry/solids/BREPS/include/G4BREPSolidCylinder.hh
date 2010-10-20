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
// $Id: G4BREPSolidCylinder.hh,v 1.12 2010-10-20 09:14:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BREPSolidCylinder
//
// Class description:
//
//  Definition of a generic BREP cylinder:
//
//  G4BREPSolidCylinder(const G4String& name,
//                      const G4ThreeVector& origin,
//                      const G4ThreeVector& axis,
//                      const G4ThreeVector& direction,		      
//                            G4double radius,
//                            G4double length)

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BREPSolidCylinder
#define __G4BREPSolidCylinder

#include "G4BREPSolid.hh"

class G4BREPSolidCylinder : public G4BREPSolid
{

 public: // with description

  G4BREPSolidCylinder(const G4String& name,
		      const G4ThreeVector& origin,
		      const G4ThreeVector& axis,
		      const G4ThreeVector& direction,		      
		      G4double radius,
		      G4double length);
    // Constructor defining the cylinder through planes
    // and circular surfaces.

  ~G4BREPSolidCylinder();
    // Empty destructor.

  G4VSolid* Clone() const;
    // Returns a pointer of a dynamically allocated copy of the solid.

  std::ostream& StreamInfo(std::ostream& os) const;
    // Streams solid contents to output stream.

 public:  // without description

  G4BREPSolidCylinder(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  G4BREPSolidCylinder(const G4BREPSolidCylinder& rhs);
  G4BREPSolidCylinder& operator=(const G4BREPSolidCylinder& rhs);
    // Copy constructor and assignment operator.

 private:

  void InitializeCylinder();

 private:

  struct G4BREPCylinderParams
  {
    G4ThreeVector origin;
    G4ThreeVector axis;
    G4ThreeVector direction;
    G4double      length;
    G4double      radius;
  } constructorParams;
};

#endif
