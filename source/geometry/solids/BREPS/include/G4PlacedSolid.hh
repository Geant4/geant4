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
// Class G4PlacedSolid
//
// Class description:
//
// Class for a generic BREP solid placed in a 3D space.
// Attributes are: the solid, a traslation vector and a
// rotation matrix.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef G4PLACEDSOLID_HH
#define G4PLACEDSOLID_HH

#include "G4BREPSolid.hh"
#include "G4RotationMatrix.hh"

class G4PlacedSolid
{

public:  // with description

  G4PlacedSolid();
  G4PlacedSolid(G4BREPSolid*, G4Axis2Placement3D* =0);
  ~G4PlacedSolid();
    // Constructors & destructor

  inline G4bool operator==(const G4PlacedSolid& ps) const;
    // Equality operator.

  inline G4VSolid*         GetSolid() const;
  inline G4RotationMatrix* GetRotation() const;
  inline G4ThreeVector*    GetTranslation() const;
    // Accessors.

  std::ostream& StreamInfo(std::ostream& os) const;
    // Streams solid contents to output stream.

private:

  G4PlacedSolid(const G4PlacedSolid&);
  G4PlacedSolid& operator=(const G4PlacedSolid&);
    // Private copy constructor and assignment operator.

private:

  G4BREPSolid*      solid;
  G4RotationMatrix* solidRotation;
  G4ThreeVector*    solidTranslation;

};

#include "G4PlacedSolid.icc"

#endif
