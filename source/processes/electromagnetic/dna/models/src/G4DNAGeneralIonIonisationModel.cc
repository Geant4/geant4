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

#include "G4DNAGeneralIonIonisationModel.hh"
#include "G4DNARuddIonisationExtendedModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAGeneralIonIonisationModel::G4DNAGeneralIonIonisationModel(
	    const G4ParticleDefinition*, const G4String& nam)
  : G4VEmModel(nam)
{
  fRuddIonisation = new G4DNARuddIonisationExtendedModel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAGeneralIonIonisationModel::Initialise(const G4ParticleDefinition* part,
                                                const G4DataVector& v)
{
  if (part != G4GenericIon::GenericIon()) {
    G4ExceptionDescription ed;
    ed << "Wrong particle type <" << part->GetParticleName()
       << "> - only G4GenericIon is allowed";  
    G4Exception("G4DNAGeneralIonIonisationModel::Initialise(...)",
                 "em2001", FatalException, ed);
  }
  
  // stationary model and tracking cut may be defined between runs
  auto param = G4EmParameters::Instance();
  fLowestEnergy = param->LowestMuHadEnergy();

  // this pointer should be defined once before initialisation of
  // concrete models
  if (nullptr == fParticleChangeForGamma) {
    fParticleChangeForGamma = GetParticleChangeForGamma();

    // should be set once for each model
    fRuddIonisation->SetParticleChange(fParticleChangeForGamma);
  }

  // initialisation of concrete models - put extra model here
  fRuddIonisation->Initialise(part, v);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAGeneralIonIonisationModel::StartTracking(G4Track* track)
{
  fDynParticle = track->GetDynamicParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAGeneralIonIonisationModel::CrossSectionPerVolume(const G4Material* mat,
                                                               const G4ParticleDefinition* p,
							       G4double ekin,
							       G4double, G4double)
{
  // apply tracking cut
  if (ekin <= fLowestEnergy) { return DBL_MAX; }

  // select model
  fCurrentModel = fRuddIonisation;

  // G4int Z = GetAtomicNumber();
  // G4int Q = G4lrint(fDynParticle->GetCharge()/CLHEP::eplus);

  // compute concrete cross section
  return fCurrentModel->CrossSectionPerVolume(mat, p, ekin);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAGeneralIonIonisationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                       const G4MaterialCutsCouple* couple,
                                                       const G4DynamicParticle* dp,
                                                       G4double, G4double)
{
  G4double ekin = dp->GetKineticEnergy();

  // tracking cut
  if (ekin <= fLowestEnergy) {
    fParticleChangeForGamma->SetProposedKineticEnergy(0.);
    fParticleChangeForGamma->ProposeTrackStatus(fStopButAlive);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(ekin);
    return;
  }

  // sample secondaries by the selected model
  fCurrentModel->SampleSecondaries(fvect, couple, dp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
