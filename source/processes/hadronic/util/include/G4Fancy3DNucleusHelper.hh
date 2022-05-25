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
#ifndef G4Fancy3DNucleusHelper_h
#define G4Fancy3DNucleusHelper_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4Fancy3DNucleusHelper ----------------
//             by Gunter Folger, May 1998.
// Simple container class to support sorting nucleon configuration
// ------------------------------------------------------------
// 20110805  M. Kelsey -- Extracted from G4Fancy3DNucleus.cc

#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Fancy3DNucleusHelper {
public:
  G4Fancy3DNucleusHelper(): Vector(0,0,0), Size(0), Index(0) {}

  G4Fancy3DNucleusHelper(const G4ThreeVector &vec, G4double size,
			 G4int index)
    : Vector(vec), Size(size), Index(index) {}

  void Fill(const G4ThreeVector &vec, G4double size, G4int index) {
    Vector = vec;
    Size = size;
    Index = index;
  }

  // Equality only true for same object (tests by address)
  G4bool operator==(const G4Fancy3DNucleusHelper &right) const
  {
    return (this == &right);
  }

  // Sorting on size value
  G4bool operator<(const G4Fancy3DNucleusHelper &right) const
  {
    return (Size < right.Size);
  }

  // Data members are public for simplicity of sorting
  G4ThreeVector Vector;
  G4double Size;
  G4int Index;
};

#endif	/* G4Fancy3DNucleusHelper_h */
