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
#include "G4LEPTSIonisationModel.hh"
#include "CLHEP/Units/PhysicalConstants.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4LEPTSIonisationModel::G4LEPTSIonisationModel(const G4String& modelName) 
  : G4VLEPTSModel( modelName )
{
  SetDeexcitationFlag(true);
  fParticleChangeForGamma = 0;
  theXSType = XSIonisation; 

} // constructor

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4LEPTSIonisationModel::~G4LEPTSIonisationModel() 
{
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4LEPTSIonisationModel::Initialise(const G4ParticleDefinition* aParticle, 
                          const G4DataVector&)
{
  Init();
  BuildPhysicsTable( *aParticle );
  fParticleChangeForGamma = GetParticleChangeForGamma();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4LEPTSIonisationModel::CrossSectionPerVolume(const G4Material* mate,
                                         const G4ParticleDefinition* aParticle,
                                         G4double kineticEnergy,
                                         G4double,
                                         G4double)
{
  return 1./GetMeanFreePath( mate, aParticle, kineticEnergy );

}


void G4LEPTSIonisationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                 const G4MaterialCutsCouple* mateCuts,
                                 const G4DynamicParticle* aDynamicParticle,
                                 G4double,
                                 G4double)
{
  G4double P0KinEn = aDynamicParticle->GetKineticEnergy();

  G4double Edep=0;
  G4double Energylost=0;
  G4ThreeVector P0Dir = aDynamicParticle->GetMomentumDirection();

  const G4Material* aMaterial = mateCuts->GetMaterial();
  if(P0KinEn < theIonisPot[aMaterial]) {
    theIonisPot[aMaterial] = P0KinEn;
  }
  Energylost = SampleEnergyLoss(aMaterial, theIonisPot[aMaterial], P0KinEn);
  G4ThreeVector P1Dir = SampleNewDirection(aMaterial, P0Dir, P0KinEn/CLHEP::eV, Energylost/CLHEP::eV);
  G4double P1KinEn = std::max(0., P0KinEn - Energylost);
  fParticleChangeForGamma->ProposeMomentumDirection( P1Dir);
  fParticleChangeForGamma->SetProposedKineticEnergy( P1KinEn);
#ifdef DEBUG_LEPTS
  G4cout << " G4LEPTSIonisationModel::SampleSecondaries( SetProposedKineticEnergy " << P1KinEn << " " << P0KinEn << " - " << Energylost << G4endl;
#endif

  G4double P2KinEn;

  if( Energylost < theIonisPotInt[aMaterial]) {  // External Ionisation
    //-    SetModelName("Ionisation");
    Edep = theIonisPot[aMaterial];
    P2KinEn = std::max(0.001*CLHEP::eV, (Energylost - theIonisPot[aMaterial]) );
  }
  else {                   // Auger    
    //-    SetModelName("IonisAuger");
    Edep = 35*CLHEP::eV;   
    P2KinEn = std::max(0.0, (Energylost - theIonisPotInt[aMaterial]) );
    G4double P3KinEn = std::max(0.0, theIonisPotInt[aMaterial] - Edep);

    G4ThreeVector P3Dir;
    P3Dir.setX( G4UniformRand() );
    P3Dir.setY( G4UniformRand() );
    P3Dir.setZ( G4UniformRand() );
    P3Dir /= P3Dir.mag();

    G4DynamicParticle* e3 = new G4DynamicParticle(G4Electron::Electron(), P3Dir, P3KinEn);
    fvect->push_back(e3);
  }

  fParticleChangeForGamma->ProposeLocalEnergyDeposit(Edep);

  if( P2KinEn > theLowestEnergyLimit) {
    G4double cp0 = std::sqrt(P0KinEn*(P0KinEn + 2.*CLHEP::electron_mass_c2));
    G4double cp1 = std::sqrt(P1KinEn*(P1KinEn + 2.*CLHEP::electron_mass_c2));
    G4ThreeVector P2Momentum = cp0*P0Dir -cp1*P1Dir;
    G4ThreeVector P2Dir = P2Momentum / P2Momentum.mag();
    P2Dir.rotateUz(P0Dir);
    G4DynamicParticle* e2 = new G4DynamicParticle(G4Electron::Electron(), P2Dir, P2KinEn);
    fvect->push_back(e2);
  }
 
}
