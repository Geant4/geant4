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
//
// File name:     G4ModifiedMephi
//
// Author:        V. Ivanchenko
// 
// Creation date: 27 October 2020
//
// Modifications: 
//
//
// -------------------------------------------------------------------
//

#include "G4ModifiedMephi.hh"
#include "Randomize.hh"
#include "G4Log.hh"
#include <CLHEP/Units/PhysicalConstants.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ModifiedMephi::G4ModifiedMephi(const G4String&)
  : G4VEmAngularDistribution("ModifiedMephi")
{}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ModifiedMephi::~G4ModifiedMephi() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector& 
G4ModifiedMephi::SampleDirection(const G4DynamicParticle* dp,
                                 G4double gEnergy, G4int, 
                                 const G4Material*)
{
  // Sample gamma angle (Z - axis along the parent particle).
  G4double cost = SampleCosTheta(dp->GetKineticEnergy(), gEnergy, 
                                 dp->GetDefinition()->GetPDGMass());
  G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));
  G4double phi  = CLHEP::twopi*G4UniformRand(); 

  fLocalDirection.set(sint*std::cos(phi), sint*std::sin(phi), cost);
  fLocalDirection.rotateUz(dp->GetMomentumDirection());

  return fLocalDirection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ModifiedMephi::SampleCosTheta(G4double primKinEnergy, 
                                         G4double gEnergy,
                                         G4double mass)
{
  G4double gam  = 1.0 + primKinEnergy/mass;
  G4double rmax = gam*CLHEP::halfpi*std::min(1.0, gam*mass/gEnergy - 1.0);
  G4double rmax2= rmax*rmax;
  G4double x = G4UniformRand()*rmax2/(1.0 + rmax2);

  return std::cos(std::sqrt(x/(1.0 - x))/gam);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ModifiedMephi::SamplePairDirections(const G4DynamicParticle* dp,
					   G4double elecKinEnergy,
					   G4double posiKinEnergy,
					   G4ThreeVector& dirElectron,
					   G4ThreeVector& dirPositron,
					   G4int, const G4Material*)
{
  G4double phi  = CLHEP::twopi * G4UniformRand();
  G4double sinp = std::sin(phi);
  G4double cosp = std::cos(phi);

  G4double ekin = dp->GetKineticEnergy();
  G4double etwo = elecKinEnergy + posiKinEnergy;
  G4double mass = dp->GetDefinition()->GetPDGMass();

  G4double cost = SampleCosTheta(ekin, etwo, mass);
  G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));

  dirElectron.set(sint*cosp, sint*sinp, cost);
  dirElectron.rotateUz(dp->GetMomentumDirection());

  cost = SampleCosTheta(ekin, etwo, mass);
  sint = std::sqrt((1.0 - cost)*(1.0 + cost));

  dirPositron.set(-sint*cosp, -sint*sinp, cost);
  dirPositron.rotateUz(dp->GetMomentumDirection());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ModifiedMephi::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Angular Generator is Modified Mephi" << G4endl;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
