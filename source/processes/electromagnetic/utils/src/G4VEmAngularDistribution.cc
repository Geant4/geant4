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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4VEmAngularDistribution
//
// Author:        V. Ivanchenko using design of existing 
//                interface G4VBremAngularDistribution
// 
// Creation date: 13 October 2010
//
// Modifications: 
//
// Class Description: 
//
// Abstract base class for polar angle sampling
//
// Class Description: End 

// -------------------------------------------------------------------
//
//    

#include "G4VEmAngularDistribution.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmAngularDistribution::G4VEmAngularDistribution(const G4String& name) 
  : fName(name)
{
  fLocalDirection.set(0.0,0.0,1.0);
  fPolarisation = G4EmParameters::Instance()->EnablePolarisation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmAngularDistribution::~G4VEmAngularDistribution() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector& G4VEmAngularDistribution::SampleDirectionForShell(
                                         const G4DynamicParticle* dp,
                                         G4double finalTotalEnergy,
                                         G4int Z, G4int,
                                         const G4Material* mat)
{
  return SampleDirection(dp, finalTotalEnergy, Z, mat); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmAngularDistribution::SamplePairDirections(const G4DynamicParticle* dp,
				                    G4double, G4double,
				                    G4ThreeVector& dirElectron,
				                    G4ThreeVector& dirPositron,
				                    G4int, const G4Material*)
{
  dirElectron = dp->GetMomentumDirection();
  dirPositron = dp->GetMomentumDirection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmAngularDistribution::PrintGeneratorInformation() const
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
