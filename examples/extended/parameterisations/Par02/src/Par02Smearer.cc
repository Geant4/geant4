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
/// \file Par02Smearer.cc
/// \brief Implementation of the Par02Smearer class

#include "Par02Smearer.hh"
#include "Par02PrimaryParticleInformation.hh"
#include "G4PrimaryParticle.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include <ctime>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02Smearer* Par02Smearer::fPar02Smearer = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02Smearer::Par02Smearer() {
  time_t seed = time( NULL );
  fRandomEngine = new CLHEP::HepJamesRandom( static_cast< long >( seed ) );
  fRandomGauss = new CLHEP::RandGauss( fRandomEngine );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02Smearer::~Par02Smearer() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02Smearer* Par02Smearer::Instance() {
  if ( ! fPar02Smearer ) {
    fPar02Smearer = new Par02Smearer();
  }
  return fPar02Smearer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector Par02Smearer::SmearMomentum( const G4Track* aTrackOriginal, 
                                           G4double aResolution ) {
  return SmearGaussian( aTrackOriginal, aResolution );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Par02Smearer::SmearEnergy( const G4Track* aTrackOriginal, 
                                    G4double aResolution ) {
  G4double newE = -1.0;
  while ( newE < 0.0 ) {  // To ensure that the resulting value is not negative
                          // (vital for energy smearing, does not change direction
                          // for momentum smearing)
    if ( aResolution != -1.0 ) {
      newE = aTrackOriginal->GetKineticEnergy() * Gauss( 1.0, aResolution );
    } else {
      newE = aTrackOriginal->GetKineticEnergy();
    }
  }
  return newE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector Par02Smearer::SmearGaussian( const G4Track* aTrackOriginal, 
                                           G4double aResolution ) {
  G4ThreeVector originP = aTrackOriginal->GetMomentum();
  G4ThreeVector originPos = aTrackOriginal->GetPosition();
  G4double rdm = Gauss( 1.0, aResolution );
  G4ThreeVector smearedMom( originP.x()*rdm, originP.y()*rdm, originP.z()*rdm );
  return smearedMom;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Par02Smearer::Gauss( G4double aMean, G4double aStandardDeviation ) {
  return fRandomGauss->fire( aMean, aStandardDeviation );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

