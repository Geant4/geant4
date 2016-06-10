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

#include "ExExChPhysicsList.hh"

// general
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessVector.hh"
#include "G4RunManager.hh"

// physics lists
#include "ExExChPhysListEmStandardSS.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "ExExChHadronPhysicsQGSP_BIC.hh"
#include "ExExChHadronElasticPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "ExExChIonPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4EmExtraPhysics.hh"

// channeling
#include "ExExChProcessChanneling.hh"

#include "XVCrystalCharacteristic.hh"
#include "XVCrystalIntegratedDensity.hh"

#include "XCrystalPlanarMoliereTempPotential.hh"
#include "XCrystalPlanarMoliereElectricField.hh"
#include "XCrystalPlanarNucleiDensity.hh"
#include "XCrystalPlanarMoliereElectronDensity.hh"
#include "XCrystalCharacteristicArray.hh"
#include "XCrystalIntegratedDensityPlanar.hh"
#include "XCrystalIntegratedDensityHub.hh"


ExExChPhysicsList::ExExChPhysicsList():  G4VModularPhysicsList(){

    fFilePotentialName = "";
    fTimeStepMin = 2.E2 * CLHEP::angstrom;
    fTransverseVariationMax = 2.E-2 * CLHEP::angstrom;    
    fParticleList = new G4DecayPhysics();
    
    fDecayList = new G4RadioactiveDecayPhysics();

    fRaddecayList = new G4DecayPhysics();
    
    fEmPhysicsList = new ExExChPhysListEmStandardSS();

    fHadronInelasticPhysicsList = new ExExChHadronPhysicsQGSP_BIC();

    fHadronElasticPhysicsList = new ExExChHadronElasticPhysics();

    fStoppingPhysics = new G4StoppingPhysics();

    fIonPhysics = new ExExChIonPhysics();

    fNeutronTrackingCut = new G4NeutronTrackingCut();

    fEmExtraPhysics = new G4EmExtraPhysics();

    fMessenger = new ExExChPhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChPhysicsList::~ExExChPhysicsList(){
    delete fMessenger;
    delete fEmPhysicsList;
    delete fDecayList;
    delete fRaddecayList;
    delete fHadronElasticPhysicsList;
    delete fHadronInelasticPhysicsList;
    delete fStoppingPhysics;
    delete fIonPhysics;
    delete fNeutronTrackingCut;
    delete fEmExtraPhysics;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChPhysicsList::ConstructProcess(){
    
    AddTransportation();

    fParticleList->ConstructProcess();
    
    fEmPhysicsList->ConstructProcess();
    
    fDecayList->ConstructProcess();
    
    fRaddecayList->ConstructProcess();
    
    fHadronElasticPhysicsList->ConstructProcess();

    fHadronInelasticPhysicsList->ConstructProcess();

    fStoppingPhysics->ConstructProcess();

    fIonPhysics->ConstructProcess();

    fNeutronTrackingCut->ConstructProcess();
    
    fEmExtraPhysics->ConstructProcess();

    AddChanneling();

    G4cout << "### ExExChPhysicsList::ConstructProcess is done" << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChPhysicsList::ConstructParticle(){
    fParticleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChPhysicsList::AddChanneling(){
    
    XVCrystalCharacteristic* vPotentialEnergy =
        new XCrystalPlanarMoliereTempPotential();
    XVCrystalCharacteristic* vElectricField =
        new XCrystalPlanarMoliereElectricField();
    XVCrystalCharacteristic* vNucleiDensity =
        new XCrystalPlanarNucleiDensity();
    XVCrystalCharacteristic* vElectronDensity =
        new XCrystalPlanarMoliereElectronDensity();
    
    XVCrystalIntegratedDensity* vIntegratedDensityNuclei =
        new XCrystalIntegratedDensityPlanar();
    vIntegratedDensityNuclei->SetPotential(vPotentialEnergy);
    vIntegratedDensityNuclei->SetDensity(vNucleiDensity);
    
    XVCrystalIntegratedDensity* vIntegratedDensityElectron =
        new XCrystalIntegratedDensityPlanar();
    vIntegratedDensityElectron->SetPotential(vPotentialEnergy);
    vIntegratedDensityElectron->SetDensity(vElectronDensity);
    
    XCrystalIntegratedDensityHub* vIntegratedDensityHub =
        new XCrystalIntegratedDensityHub();
    vIntegratedDensityHub->SetPotential(vPotentialEnergy);
    vIntegratedDensityHub->SetDensityNuclei(vNucleiDensity);
    vIntegratedDensityHub->SetDensityElectron(vElectronDensity);
    
    for(G4int i=-3;i<=+3;i++){
        if(i!=0){
            vIntegratedDensityHub->SetIntegratedDensityNuclei(
                                new XCrystalIntegratedDensityPlanar(),i);
            vIntegratedDensityHub->SetIntegratedDensityElectron(
                                new XCrystalIntegratedDensityPlanar(),i);
        }
    }
    
    ExExChProcessChanneling* channeling =
        new ExExChProcessChanneling();
    channeling->SetPotential(vPotentialEnergy);
    channeling->SetIntegratedDensity(vIntegratedDensityHub);
    channeling->SetElectricField(vElectricField);
    channeling->SetNucleiDensity(vNucleiDensity);
    channeling->SetElectronDensity(vElectronDensity);
    
    channeling->SetTransverseVariationMax(fTransverseVariationMax);
    channeling->SetTimeStepMin(fTimeStepMin);
    if(fFilePotentialName != ""){
        channeling->SetFileCharacteristicsName(fFilePotentialName);
    }
    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pManager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        
        if(particle->GetPDGCharge() != 0)
        {
            pManager->AddDiscreteProcess(channeling);
        }
    }
    G4cout << "\nPhysicsList::AddChanneling: Channeling process added...\n";
    G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void ExExChPhysicsList::SetCuts()
{
    // These values are used as the default production thresholds
    // for the world volume.
    SetCutsWithDefault();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChPhysicsList::SetFilePotentialName(const G4String& vFilename){
    if(fFilePotentialName != vFilename){
        fFilePotentialName = vFilename;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String ExExChPhysicsList::GetFilePotentialName(){
    return fFilePotentialName;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

