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
// File name:   G4EmMultiModel
//
// Author:        Vladimir Ivanchenko
// 
// Creation date: 03.05.2004
//
// Modifications: 
// 15-04-05 optimize internal interface (V.Ivanchenko)
// 04-07-10 updated interfaces according to g4 9.4 (V.Ivanchenko)
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmMultiModel.hh"
#include "Randomize.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmMultiModel::G4EmMultiModel(const G4String& nam)
  : G4VEmModel(nam)
{
  model.clear();
  cross_section.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmMultiModel::~G4EmMultiModel() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmMultiModel::AddModel(G4VEmModel* p)
{
  cross_section.push_back(0.0);
  model.push_back(p);
  ++nModels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmMultiModel::Initialise(const G4ParticleDefinition* p, 
                                const G4DataVector& cuts)
{
  G4EmParameters* param = G4EmParameters::Instance();
  G4int verb = IsMaster() ? param->Verbose() : param->WorkerVerbose();
  if(verb > 0) {
    G4cout << "### Initialisation of EM MultiModel " << GetName()
	   << " including following list of " << nModels << " models:" << G4endl;
  }
  for(G4int i=0; i<nModels; ++i) {
    G4cout << "    " << (model[i])->GetName();
    (model[i])->SetParticleChange(pParticleChange, GetModelOfFluctuations());
    (model[i])->Initialise(p, cuts);
  }
  if(verb > 0) { G4cout << G4endl; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmMultiModel::ComputeDEDXPerVolume(const G4Material* mat,
					      const G4ParticleDefinition* p,
					      G4double kineticEnergy,
					      G4double cutEnergy)
{
  G4double dedx = 0.0;
  for(G4int i=0; i<nModels; ++i) {
    dedx += (model[i])->ComputeDEDXPerVolume(mat, p, cutEnergy, kineticEnergy);
  } 

  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmMultiModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition* p,
                                                    G4double kinEnergy,
                                                    G4double Z,
                                                    G4double A,
                                                    G4double cutEnergy,
                                                    G4double maxEnergy)
{
  G4double cross = 0.0;
  for(G4int i=0; i<nModels; ++i) {
    (model[i])->SetCurrentCouple(CurrentCouple());
    cross += (model[i])->ComputeCrossSectionPerAtom(p, kinEnergy, Z, A, 
						    cutEnergy, maxEnergy);
  } 
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmMultiModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
                                       const G4MaterialCutsCouple* couple,
                                       const G4DynamicParticle* dp,
                                       G4double minEnergy,
                                       G4double maxEnergy)
{
  SetCurrentCouple(couple);
  if(nModels > 0) {
    G4int i;
    G4double cross = 0.0;
    for(i=0; i<nModels; ++i) {
      cross += (model[i])->CrossSection(couple, dp->GetParticleDefinition(), 
                                        dp->GetKineticEnergy(), minEnergy, maxEnergy);
      cross_section[i] = cross;
    }

    cross *= G4UniformRand();

    for(i=0; i<nModels; ++i) {
      if(cross <= cross_section[i]) {
        (model[i])->SampleSecondaries(vdp, couple, dp, minEnergy, maxEnergy);
        return;
      }
    }
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

