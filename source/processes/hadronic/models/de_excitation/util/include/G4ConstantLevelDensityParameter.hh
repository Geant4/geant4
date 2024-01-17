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
//
// Hadronic Process: Nuclear De-excitations (photon evaporation)
// by C. Dallapiccola (Nov 1998) 
//

#ifndef G4ConstantLevelDensityParameter_h
#define G4ConstantLevelDensityParameter_h 1

#include "G4VLevelDensityParameter.hh"

class G4ConstantLevelDensityParameter : public G4VLevelDensityParameter
{
public:

  G4ConstantLevelDensityParameter() = default;
  ~G4ConstantLevelDensityParameter() = default;

  G4double LevelDensityParameter(G4int A, G4int, G4double) const final
  { 
    return 0.1*A;
  }

  G4ConstantLevelDensityParameter(const G4ConstantLevelDensityParameter &right) = delete;
  const G4ConstantLevelDensityParameter& operator=
  (const G4ConstantLevelDensityParameter &right) = delete;
  G4bool operator==(const G4ConstantLevelDensityParameter &right) const = delete;
  G4bool operator!=(const G4ConstantLevelDensityParameter &right) const = delete;
};

#endif
