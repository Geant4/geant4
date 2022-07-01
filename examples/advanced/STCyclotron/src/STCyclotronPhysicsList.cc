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
// Author: F. Poignant, floriane.poignant@gmail.com
//

#include "STCyclotronPhysicsList.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4Deuteron.hh"

//Physic Lists (contained inside the Geant 4 distribution)
#include "G4EmStandardPhysics_option3.hh"
#include "G4DecayPhysics.hh"
#include "G4Decay.hh"
#include "G4StepLimiter.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4EmExtraPhysics.hh"
#include "G4EmParameters.hh"
#include "G4NuclideTable.hh"

#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"
#include "G4RadioactiveDecayPhysics.hh"
//#include "QGSP_BIC_HP.hh"
#include "G4PhysListFactory.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclideTable.hh"

STCyclotronPhysicsList::STCyclotronPhysicsList(STCyclotronDetectorConstruction* det)
  : G4VModularPhysicsList()
{
  // Mandatory for G4NuclideTable
  // Half-life threshold must be set small or many short-lived isomers 
  // will not be assigned life times (default to 0) 
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*second);
  G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);

  //---
  fDetector = det;
  
  G4LossTableManager::Instance();
 
  defaultCutValue = 0.1*mm;
  fCutForGamma     = defaultCutValue;
  fCutForElectron  = defaultCutValue;
  fCutForPositron  = defaultCutValue;

  fThickness_foil = defaultCutValue;
  fThickness_target = defaultCutValue;

  fCutTargetProton  = fThickness_target;
  fCutTargetElectron = fThickness_target;
  fCutTargetPositron = fThickness_target;
  fCutTargetGamma = fThickness_target;

  fCutFoilProton = fThickness_foil;
  fCutFoilElectron = fThickness_foil;
  fCutFoilPositron = fThickness_foil;
  fCutFoilGamma = fThickness_foil;
  
  //EM physics
  fEmPhysicsList = new G4EmStandardPhysics_option3(0);
  fEmName = G4String("emstandard_opt3");

  //Decay physics and all particles
  fDecPhysicsList = new G4DecayPhysics(0);
  fRaddecayList = new G4RadioactiveDecayPhysics(0);
  
  //Hadron physics
  fHadPhysicsList = new G4HadronPhysicsQGSP_BIC_AllHP(0);
  //fHadPhysicsList = new G4HadronPhysicsQGSP_BIC(0);


  
  

}

STCyclotronPhysicsList::~STCyclotronPhysicsList()
{

  delete fEmPhysicsList ;
  delete fDecPhysicsList;
  delete fRaddecayList;
  delete fHadPhysicsList;

  
}

void STCyclotronPhysicsList::ConstructParticle()
{

  G4Proton::ProtonDefinition();
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4Neutron::NeutronDefinition();
  G4Deuteron::DeuteronDefinition();

  fDecPhysicsList->ConstructParticle();
  

}

void STCyclotronPhysicsList::ConstructProcess()
{
  // Define transportation process
  

  AddTransportation();
  
  //electromagnetic physics list
  fEmPhysicsList->ConstructProcess();
  //em_config.AddModels();
  
  //decay physics list
  fDecPhysicsList->ConstructProcess();
  fRaddecayList->ConstructProcess();

  //hadronic physics lists
  fHadPhysicsList->ConstructProcess();

  //Get the value of the fThickness of foil and target
  fThickness_foil = fDetector->GetFoilThickness()*mm;
  fThickness_target = fDetector->GetTargetThickness()*mm;

  //Update the cuts with the 1/2 of the thickness of the foil/target
  //SetCuts();

  SetCutTarget(0.01,fThickness_target/2.,fThickness_target/2.,fThickness_target/2.);
  
  

}

void STCyclotronPhysicsList::SetCuts()
{

  if (verboseLevel >0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(fCutForGamma, "gamma");
  SetCutValue(fCutForElectron, "e-");
  SetCutValue(fCutForPositron, "e+");

  // Set cuts for detector
  SetCutFoil(0.1,fThickness_foil/2.,fThickness_foil/2.,fThickness_foil/2.);
  SetCutTarget(0.01,fThickness_target/2.,fThickness_target/2.,fThickness_target/2.);
  if (verboseLevel>0) DumpCutValuesTable();
}

void STCyclotronPhysicsList::SetCutForGamma(G4double cut)
{
  fCutForGamma = cut;
  SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

void STCyclotronPhysicsList::SetCutForElectron(G4double cut)
{
  fCutForElectron = cut;
  SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

void STCyclotronPhysicsList::SetCutForPositron(G4double cut)
{
  fCutForPositron = cut;
  SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

void STCyclotronPhysicsList::SetCutTarget(G4double cutProton, G4double cutElectron, G4double cutPositron, G4double cutGamma){
  
  fCutTargetProton = cutProton*mm;
  fCutTargetElectron = cutElectron*mm;
  fCutTargetPositron = cutPositron*mm;
  fCutTargetGamma = cutGamma*mm;

  G4String regionNameTarget = "Target";
  G4Region* regionTarget = G4RegionStore::GetInstance()->GetRegion(regionNameTarget);

  G4ProductionCuts* cutsTarget = new G4ProductionCuts ;
  cutsTarget -> SetProductionCut(fCutTargetGamma,"gamma");
  cutsTarget -> SetProductionCut(fCutTargetElectron,"e-");
  cutsTarget -> SetProductionCut(fCutTargetPositron,"e+");
  cutsTarget -> SetProductionCut(fCutTargetProton,"proton");
  
  regionTarget -> SetProductionCuts(cutsTarget);
    
}

void STCyclotronPhysicsList::SetCutFoil(G4double cutProton, G4double cutElectron, G4double cutPositron, G4double cutGamma){
  
  fCutFoilProton = cutProton*mm;
  fCutFoilElectron = cutElectron*mm;
  fCutFoilPositron = cutPositron*mm;
  fCutFoilGamma = cutGamma*mm;

  G4RegionStore::GetInstance()->GetRegion("Foil");
  G4String regionNameFoil = "Foil";
  G4Region* regionFoil = G4RegionStore::GetInstance()->GetRegion(regionNameFoil);
  
  G4ProductionCuts* cutsFoil = new G4ProductionCuts ;
  cutsFoil -> SetProductionCut(fCutFoilGamma,"gamma");
  cutsFoil -> SetProductionCut(fCutFoilElectron,"e-");
  cutsFoil -> SetProductionCut(fCutFoilPositron,"e+");
  cutsFoil -> SetProductionCut(fCutFoilProton,"proton");

  regionFoil -> SetProductionCuts(cutsFoil);

}

