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
/// \file GB05SD.cc
/// \brief Implementation of the GB05SD class

#include "GB05SD.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB05SD::GB05SD(G4String name)
: G4VSensitiveDetector(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool GB05SD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto track        = step->GetTrack();
  auto preStepPoint = step->GetPreStepPoint();

  // -- simply prints few particle characteristics:
  G4cout << std::setw(14)
         << track->GetParticleDefinition()->GetParticleName()
         << ", kinetic energy (MeV) = "
         << std::setw(12)
         << preStepPoint->GetKineticEnergy()/MeV
         << ", position (cm) = "
         << preStepPoint->GetPosition()/cm
         << ",\t weight = "
         << preStepPoint->GetWeight()
         << G4endl;
  
  return true;
}
