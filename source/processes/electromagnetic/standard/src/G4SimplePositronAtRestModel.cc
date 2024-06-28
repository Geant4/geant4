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
// GEANT4 Class file
//
//
// File name:     G4SimplePositronAtRestModel
//
// Author:        Vladimir Ivanchenko
// 
// Creation date: 14 May 2024
//
// -------------------------------------------------------------------
//

#include "G4SimplePositronAtRestModel.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4RandomDirection.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4SimplePositronAtRestModel::G4SimplePositronAtRestModel()
  : G4VPositronAtRestModel("Simple")
{}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4SimplePositronAtRestModel::SampleSecondaries(
             std::vector<G4DynamicParticle*>& secParticles,
             G4double&, const G4Material*) const
{
  G4ThreeVector dir1 = G4RandomDirection();
  auto aGamma1 = new G4DynamicParticle(G4Gamma::Gamma(), dir1,
				       CLHEP::electron_mass_c2);
  G4double phi = CLHEP::twopi * G4UniformRand();
  G4double cosphi = std::cos(phi);
  G4double sinphi = std::sin(phi);
  G4ThreeVector pol1(cosphi, sinphi, 0.0);
  pol1.rotateUz(dir1);
  aGamma1->SetPolarization(pol1);
  secParticles.push_back(aGamma1);

  G4ThreeVector dir2 = -dir1;
  auto aGamma2 = new G4DynamicParticle(G4Gamma::Gamma(), dir2,
				       CLHEP::electron_mass_c2);
  G4ThreeVector pol2(-sinphi, cosphi, 0.0);
  pol2.rotateUz(dir1);
  aGamma2->SetPolarization(pol2);
  secParticles.push_back(aGamma2);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4SimplePositronAtRestModel::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Simple AtRest positron 2-gamma annihilation model" << G4endl;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
