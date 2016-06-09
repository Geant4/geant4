//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "exrdmPhysicsList.hh"
#include "exrdmPhysicsListMessenger.hh"

#include "exrdmPhysListParticles.hh"
#include "exrdmPhysListGeneral.hh"
#include "exrdmPhysListEmStandard.hh"
#include "exrdmPhysListEmLowEnergy.hh"
#include "exrdmPhysListHadron.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"

#include "HadronPhysicsQGSP_BERT.hh"
#include "HadronPhysicsQGSP_BIC.hh"
#include "HadronPhysicsQGSP_HP.hh"
#include "HadronPhysicsLHEP_BERT.hh"
#include "HadronPhysicsLHEP_BERT_HP.hh"
#include "HadronPhysicsLHEP_BIC.hh"
#include "HadronPhysicsLHEP_BIC_HP.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmPhysicsList::exrdmPhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  DetectorCuts = 0;
  TargetCuts   = 0;

  pMessenger = new exrdmPhysicsListMessenger(this);

  SetVerboseLevel(1);
  //default physics
   // Particles
  particleList = new exrdmPhysListParticles("particles");

  // General exrdmPhysics
  generalPhysicsList = new exrdmPhysListGeneral("general");

  // EM physics
  emPhysicsList = new exrdmPhysListEmStandard("em-standard");
  
  // Had physics (no hadron as default)
  hadPhysicsList =  0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmPhysicsList::~exrdmPhysicsList()
{
  delete pMessenger;
  delete generalPhysicsList;
  delete emPhysicsList;
  if (hadPhysicsList) delete hadPhysicsList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmPhysicsList::ConstructParticle()
{
  particleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmPhysicsList::ConstructProcess()
{
  AddTransportation();
  //general
   generalPhysicsList->ConstructProcess();
  // em
   emPhysicsList->ConstructProcess();
  // had
   if (hadPhysicsList) hadPhysicsList->ConstructProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmPhysicsList::SelectPhysicsList(const G4String& name)
{
  if (verboseLevel>1) {
    G4cout << "exrdmPhysicsList::SelectPhysicsList: <" << name << ">" << G4endl;
  }
  // default  Had physics
  if (name == "Hadron") {
    if (!hadPhysicsList) delete hadPhysicsList;
    hadPhysicsList = new exrdmPhysListHadron("hadron");
  } else if (name == "QGSP_BERT") {
    if (!hadPhysicsList) delete hadPhysicsList;
    hadPhysicsList = new HadronPhysicsQGSP_BERT("hadron");
  } else if (name == "QGSP_BIC") {
    if (!hadPhysicsList) delete hadPhysicsList;
    hadPhysicsList = new HadronPhysicsQGSP_BIC("hadron");
  } else if (name == "QGSP_HP") {
    if (!hadPhysicsList) delete hadPhysicsList;
    hadPhysicsList = new HadronPhysicsQGSP_HP("hadron");
  } else if (name == "LHEP_BERT") {
    if (!hadPhysicsList) delete hadPhysicsList;
    hadPhysicsList = new HadronPhysicsLHEP_BERT("hadron");
  } else if (name == "LHEP_BERT_HP") {
    if (!hadPhysicsList) delete hadPhysicsList;
    hadPhysicsList = new HadronPhysicsLHEP_BERT_HP("hadron");
  } else if (name == "LHEP_BIC") {
    if (!hadPhysicsList) delete hadPhysicsList;
    hadPhysicsList = new HadronPhysicsLHEP_BIC("hadron");
  } else if (name == "LHEP_BIC_HP") {
    if (!hadPhysicsList) delete hadPhysicsList;
    hadPhysicsList = new HadronPhysicsLHEP_BIC_HP("hadron");
  } else if (name == "LowEnergy_EM") {
    if (!hadPhysicsList) delete emPhysicsList;
    emPhysicsList = new exrdmPhysListEmLowEnergy("lowe-em");
  } else if (name == "Standard_EM") {
    if (!hadPhysicsList) delete emPhysicsList;
    emPhysicsList = new exrdmPhysListEmStandard("standard-em");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmPhysicsList::SetCuts()
{

  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  G4cout << "world cuts are set" << G4endl;

  if( !TargetCuts ) SetTargetCut(cutForElectron);
  G4Region* region = (G4RegionStore::GetInstance())->GetRegion("Target");
  region->SetProductionCuts(TargetCuts);
  G4cout << "Target cuts are set" << G4endl;

  if( !DetectorCuts ) SetDetectorCut(cutForElectron);
  region = (G4RegionStore::GetInstance())->GetRegion("Detector");
  region->SetProductionCuts(DetectorCuts);
  G4cout << "Detector cuts are set" << G4endl;

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmPhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmPhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmPhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmPhysicsList::SetTargetCut(G4double cut)
{
  if( !TargetCuts ) TargetCuts = new G4ProductionCuts();

  TargetCuts->SetProductionCut(cut, idxG4GammaCut);
  TargetCuts->SetProductionCut(cut, idxG4ElectronCut);
  TargetCuts->SetProductionCut(cut, idxG4PositronCut);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmPhysicsList::SetDetectorCut(G4double cut)
{
  if( !DetectorCuts ) DetectorCuts = new G4ProductionCuts();

  DetectorCuts->SetProductionCut(cut, idxG4GammaCut);
  DetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
  DetectorCuts->SetProductionCut(cut, idxG4PositronCut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
