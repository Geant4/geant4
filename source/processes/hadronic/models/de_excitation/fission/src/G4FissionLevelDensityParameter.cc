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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// J.M.Quesada (July 2009):  fission level density parameter tuned for spallation data.
// V.Ivanchenko (July 2009): fixed logic and remove wrong usage of MeV
// J.M.Quesada (30.10.09):   retuning for IAEA spallation data

#include "G4FissionLevelDensityParameter.hh"
#include "G4NuclearLevelData.hh"

G4FissionLevelDensityParameter::G4FissionLevelDensityParameter()
{
  fNucData = G4NuclearLevelData::GetInstance();
}

G4FissionLevelDensityParameter::~G4FissionLevelDensityParameter()
{}

G4double G4FissionLevelDensityParameter::
LevelDensityParameter(G4int A, G4int Z, G4double U) const 
{
  G4double EvapLDP = fNucData->GetLevelDensity(Z, A, U); 

  if(Z >= 89)      { EvapLDP *= 1.05; }
  else if(Z <= 85) { EvapLDP *= 1.03; }
  else             { EvapLDP *= (1.03 + 0.005*(Z - 85)); }

  return EvapLDP;
}
