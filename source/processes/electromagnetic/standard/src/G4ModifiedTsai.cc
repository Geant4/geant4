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
// File name:     G4ModifiedTsai
//
// Author:        Andreia Trindade (andreia@lip.pt)
//                Pedro Rodrigues  (psilva@lip.pt)
//                Luis Peralta     (luis@lip.pt)
// 
// Creation date: 21 March 2003
//
// Modifications: 
// 21 Mar 2003 A.Trindade First implementation acording with new design
// 24 Mar 2003 A.Trindade Fix in Tsai generator in order to prevent theta 
//                        generation above pi
// 13 Oct 2010  V.Ivanchenko  Moved to standard and improved comment
// 26.04.2011   V.Grichine    Clean-up of PolarAngle method 
//
// Class Description: 
//
// Bremsstrahlung Angular Distribution Generation 
// suggested by L.Urban (Geant3 manual (1993) Phys211)
// Derived from Tsai distribution (Rev Mod Phys 49,421(1977))
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#include "G4ModifiedTsai.hh"
#include "Randomize.hh"
#include "G4Log.hh"
#include <CLHEP/Units/PhysicalConstants.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ModifiedTsai::G4ModifiedTsai(const G4String&)
  : G4VEmAngularDistribution("ModifiedTsai")
{}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ModifiedTsai::~G4ModifiedTsai() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector& 
G4ModifiedTsai::SampleDirection(const G4DynamicParticle* dp,
                                G4double, G4int, const G4Material*)
{
  // Sample gamma angle (Z - axis along the parent particle).
  G4double cost = SampleCosTheta(dp->GetKineticEnergy());
  G4double sint = std::sqrt((1 - cost)*(1 + cost));
  G4double phi  = CLHEP::twopi*G4UniformRand(); 

  fLocalDirection.set(sint*std::cos(phi), sint*std::sin(phi), cost);
  fLocalDirection.rotateUz(dp->GetMomentumDirection());

  return fLocalDirection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ModifiedTsai::SampleCosTheta(G4double kinEnergy)
{
  // Universal distribution suggested by L. Urban (Geant3 manual (1993) 
  // Phys211) derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  G4double uMax = 2*(1. + kinEnergy/CLHEP::electron_mass_c2);   

  static const G4double a1     = 1.6;
  static const G4double a2     = a1/3.;
  static const G4double border = 0.25;
  G4double u;
  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();

  do {
    G4double uu = -G4Log(rndmEngine->flat()*rndmEngine->flat());
    u = (border > rndmEngine->flat()) ? uu*a1 : uu*a2;
    
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while(u > uMax);

  return 1.0 - 2.0*u*u/(uMax*uMax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ModifiedTsai::SamplePairDirections(const G4DynamicParticle* dp,
					  G4double elecKinEnergy,
					  G4double posiKinEnergy,
					  G4ThreeVector& dirElectron,
					  G4ThreeVector& dirPositron,
					  G4int, const G4Material*)
{
  G4double phi  = CLHEP::twopi * G4UniformRand();
  G4double sinp = std::sin(phi);
  G4double cosp = std::cos(phi);

  G4double cost = SampleCosTheta(elecKinEnergy);
  G4double sint = std::sqrt((1. - cost)*(1. + cost));

  dirElectron.set(sint*cosp, sint*sinp, cost);
  dirElectron.rotateUz(dp->GetMomentumDirection());

  cost = SampleCosTheta(posiKinEnergy);
  sint = std::sqrt((1. - cost)*(1. + cost));

  dirPositron.set(-sint*cosp, -sint*sinp, cost);
  dirPositron.rotateUz(dp->GetMomentumDirection());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ModifiedTsai::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Bremsstrahlung Angular Generator is Modified Tsai" << G4endl;
  G4cout << "Distribution suggested by L.Urban (Geant3 manual (1993) Phys211)" 
         << G4endl;
  G4cout << "Derived from Tsai distribution (Rev Mod Phys 49,421(1977)) \n" 
         << G4endl;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
