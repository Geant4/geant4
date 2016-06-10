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

#include "ExExChParticleUserInfo.hh"

ExExChParticleUserInfo::ExExChParticleUserInfo(){
    bHasBeenUnderCoherentEffect = 0;
    
    fNucleiDensity = 1.0;
    fNucleiDensityPreviousStep = 1.0;
    
    fElectronDensity = 1.0;
    fElectronDensityPreviousStep = 1.0;
    
    fNumberOfDechanneling = 0;
    
    fMomentumInChanneling = G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX);
    fPositionInChanneling = G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX);
    
    fMomentumInChannelingInitial = G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX);
    fPositionInChannelingInitial = G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX);
    fInTheCrystal = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChParticleUserInfo::~ExExChParticleUserInfo(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChParticleUserInfo::SetCoherentEffect(G4int flag){
    bHasBeenUnderCoherentEffect = flag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ExExChParticleUserInfo::HasBeenUnderCoherentEffect(){
    return bHasBeenUnderCoherentEffect;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChParticleUserInfo::SetNucleiDensity(G4double density){
    fNucleiDensity = density;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ExExChParticleUserInfo::GetNucleiDensity(){
    return fNucleiDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChParticleUserInfo::SetElectronDensity(G4double density){
    fElectronDensity = density;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ExExChParticleUserInfo::GetElectronDensity(){
    return fElectronDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ExExChParticleUserInfo::GetNucleiDensityPreviousStep(){
    return fNucleiDensityPreviousStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ExExChParticleUserInfo::GetElectronDensityPreviousStep(){
    return fElectronDensityPreviousStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChParticleUserInfo::StoreDensityPreviousStep(){
    fElectronDensityPreviousStep = fElectronDensity;
    fNucleiDensityPreviousStep = fNucleiDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector ExExChParticleUserInfo::GetMomentumChanneled(){
    return fMomentumInChanneling;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChParticleUserInfo::SetMomentumChanneled(
                                        G4ThreeVector momentum){
    fMomentumInChanneling = momentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector ExExChParticleUserInfo::GetPositionChanneled(){
    return fPositionInChanneling;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChParticleUserInfo::SetPositionChanneled(
                                        G4ThreeVector position){
    fPositionInChanneling = position;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector ExExChParticleUserInfo::GetMomentumChanneledInitial(){
    return fMomentumInChannelingInitial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChParticleUserInfo::SetMomentumChanneledInitial(
                                        G4ThreeVector momentum){
    fMomentumInChannelingInitial = momentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector ExExChParticleUserInfo::GetPositionChanneledInitial(){
    return fPositionInChannelingInitial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChParticleUserInfo::SetPositionChanneledInitial(
                                        G4ThreeVector position){
    fPositionInChannelingInitial = position;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ExExChParticleUserInfo::GetNumberOfDechanneling(){
    return fNumberOfDechanneling;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChParticleUserInfo::IncreaseNumberOfDechanneling(){
    fNumberOfDechanneling++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
