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
/// \file Par02PrimaryParticleInformation.cc
/// \brief Implementation of the Par02PrimaryParticleInformation class

#include "Par02PrimaryParticleInformation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02PrimaryParticleInformation::Par02PrimaryParticleInformation( 
  G4int aPartID, G4int aPDG, G4ThreeVector aMomentum ) : 
  fPartID( aPartID ), fPDG( aPDG ), fMomentumMC( aMomentum ), 
  fMomentumTracker( 0 ), fResolutionTracker( 0 ), fEfficiencyTracker( 0 ), 
  fPositionEMCal( 0 ), fEnergyEMCal( 0 ), fResolutionEMCal( 0 ), fEfficiencyEMCal( 0 ),
  fPositionHCal( 0 ), fEnergyHCal( 0 ), fResolutionHCal( 0 ), fEfficiencyHCal( 0 ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02PrimaryParticleInformation::~Par02PrimaryParticleInformation() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02PrimaryParticleInformation::Print() const {
  G4cout << "Par02PrimaryParticleInformation: PDG code " << fPDG << G4endl
         << "Particle unique ID:  " << fPartID << G4endl
         << "MC momentum: " << fMomentumMC << G4endl
         << "Tracker momentum: " << fMomentumTracker << G4endl
         << "Tracker resolution: " << fResolutionTracker << G4endl
         << "Tracker efficiency: " << fEfficiencyTracker << G4endl
         << "EMCal energy: " << fEnergyEMCal << " at " << fPositionEMCal << G4endl
         << "EMCal resolution: " << fResolutionEMCal << G4endl
         << "EMCal efficiency: " << fEfficiencyEMCal << G4endl
         << "HCal energy: " << fEnergyHCal << " at "<< fPositionHCal << G4endl
         << "HCal resolution: " << fResolutionHCal << G4endl
         << "HCal efficiency: " << fEfficiencyHCal << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

