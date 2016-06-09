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
// $Id: PhysicsList.cc,v 1.11 2006/06/29 17:24:21 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
/////////////////////////////////////////////////////////////////////////
//
// PhysicsList
//
// Created: 31.04.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
// 

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics72.hh"
#include "G4EmStandardPhysics71.hh"
#include "G4EmStandardPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4QStoppingPhysics.hh"
#include "G4LHEPStoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmProcessOptions.hh"

#include "HadronPhysicsQGSP_BIC.hh"
#include "HadronPhysicsLHEP.hh"
#include "HadronPhysicsLHEP_BIC.hh"
#include "HadronPhysicsLHEP_BIC_HP.hh"
#include "HadronPhysicsLHEP_BERT.hh"
#include "HadronPhysicsLHEP_BERT_HP.hh"
#include "HadronPhysicsQGSP_BERT_HP.hh"
#include "HadronPhysicsQGSP_BERT.hh"
#include "HadronPhysicsQGSC.hh"
#include "HadronPhysicsQGSP.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronInelasticQLHEP.hh"
#include "G4IonPhysics.hh"
#include "G4HadronProcessStore.hh"

#include "G4LossTableManager.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  pMessenger = new PhysicsListMessenger(this);

  // Particles
  particleList = new G4DecayPhysics("decays");

  // EM physics
  emPhysicsList = new G4EmStandardPhysics("standard");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PhysicsList::~PhysicsList()
{
  delete pMessenger;
  delete particleList;
  delete emPhysicsList;
  for(size_t i=0; i<hadronPhys.size(); i++) {
    delete hadronPhys[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::ConstructParticle()
{
  particleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  emPhysicsList->ConstructProcess();
  particleList->ConstructProcess();
  for(size_t i=0; i<hadronPhys.size(); i++) {
    hadronPhys[i]->ConstructProcess();
  }

  G4HadronProcessStore::Instance()->Dump(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0)
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;

  if (name == "LHEP") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("LElastic",
						     verboseLevel,false));
    hadronPhys.push_back( new HadronPhysicsLHEP());
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else if (name == "em_fast") {

    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics72("standard");

  } else if (name == "em_71") {

    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics71("standard");

  } else if (name == "LHEP_BERT") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("elastic",
						     verboseLevel,false));
    hadronPhys.push_back( new HadronPhysicsLHEP_BERT());
    hadronPhys.push_back( new G4QStoppingPhysics("stopping"));
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else if (name == "LHEP_BERT_HP") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("elastic",
						     verboseLevel,false));
    hadronPhys.push_back( new HadronPhysicsLHEP_BERT_HP());
    hadronPhys.push_back( new G4QStoppingPhysics("stopping"));
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else if (name == "LHEP_BIC") {

    SetStandardList(false);
    hadronPhys.push_back( new HadronPhysicsLHEP_BIC());

  } else if (name == "LHEP_BIC_HP") {

    SetStandardList(true);
    hadronPhys.push_back( new HadronPhysicsLHEP_BIC_HP());

  } else if (name == "QGSC") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("elastic",
						     verboseLevel,false));
    hadronPhys.push_back( new HadronPhysicsQGSC());
    hadronPhys.push_back( new G4QStoppingPhysics("stopping"));
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else if (name == "QGSP") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("elastic",
						     verboseLevel,false));
    hadronPhys.push_back( new HadronPhysicsQGSP());
    hadronPhys.push_back( new G4QStoppingPhysics("stopping"));
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else if (name == "LHEP_EMV") {

    AddPhysicsList("em_71");
    AddPhysicsList("LHEP");

  } else if (name == "QGSP_EMV") {

    AddPhysicsList("em_71");
    AddPhysicsList("QGSP");

  } else if (name == "QGSP_EMX") {

    AddPhysicsList("em_fast");
    AddPhysicsList("QGSP");

  } else if (name == "QGSP_BERT") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("elastic",
						     verboseLevel,false));
    hadronPhys.push_back( new HadronPhysicsQGSP_BERT());
    hadronPhys.push_back( new G4QStoppingPhysics("stopping"));
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else if (name == "QGSP_BERT_HP") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("elastic",
						     verboseLevel,true));
    hadronPhys.push_back( new HadronPhysicsQGSP_BERT_HP());
    hadronPhys.push_back( new G4QStoppingPhysics("stopping"));
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else if (name == "QGSP_BIC") {

    SetStandardList(false);
    hadronPhys.push_back( new HadronPhysicsQGSP_BIC());

  } else if (name == "QBBC") {

    SetStandardList(false);
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel));

  } else if (name == "QBBC_BERT") {

    SetStandardList(false);
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,true,false,false));

  } else if (name == "QBBC_HP") {

    SetStandardList(true);
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,false,true));

  } else if (name == "QBBC_CHIPS") {

    SetStandardList(false);
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,true,false));

  } else if (name == "QBBC_FTF") {

    SetStandardList(false);
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    true,false,false,false));

  } else if (name == "QQGSP") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("elastic",
						     verboseLevel,false));
    hadronPhys.push_back( new G4HadronInelasticQLHEP("QGSP",verboseLevel,
						 true,false,false,false));
    hadronPhys.push_back( new G4QStoppingPhysics("stopping"));
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else if (name == "QQGSP_BERT") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("elastic",
						     verboseLevel,false));
    hadronPhys.push_back( new G4HadronInelasticQLHEP("QGSP_BERT",verboseLevel,
						 true,true,false,false));
    hadronPhys.push_back( new G4QStoppingPhysics("stopping"));
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetStandardList(G4bool flagHP)
{
  hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
  hadronPhys.push_back( new G4HadronElasticPhysics("elastic",verboseLevel,
						   flagHP));
  hadronPhys.push_back( new G4QStoppingPhysics("stopping"));
  hadronPhys.push_back( new G4IonBinaryCascadePhysics("binary_ion"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetCuts()
{

  if (verboseLevel >0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

void PhysicsList::List()
{
  G4cout << "### PhysicsLists available: LHEP LHEP_BERT LHEP_BERT_HP LHEP_BIC LHEP_BIC_HP" << G4endl; 
  G4cout << "                            QGSC QGSP QGSP_BERT QGSP_BERT_HP QGSP_BIC QBBC" << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

