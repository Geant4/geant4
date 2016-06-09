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
// $Id: PhysicsList.cc,v 1.24 2007/12/12 12:00:24 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
/////////////////////////////////////////////////////////////////////////
//
// PhysicsList
//
// Created: 31.04.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
// 26.04.2007 Physics according to 8.3 Physics List (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
// 

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronDElasticPhysics.hh"
#include "G4HadronQElasticPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4QStoppingPhysics.hh"
#include "G4LHEPStoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmProcessOptions.hh"

#include "HadronPhysicsFTFC.hh"
#include "HadronPhysicsFTFP.hh"
#include "HadronPhysicsFTFP_BERT.hh"
#include "HadronPhysicsLHEP.hh"
#include "HadronPhysicsLHEP_BERT.hh"
#include "HadronPhysicsLHEP_EMV.hh"
#include "HadronPhysicsLHEP_PRECO_HP.hh"
#include "G4HadronInelasticQBBC.hh"
#include "HadronPhysicsQGSC.hh"
#include "HadronPhysicsQGSC_BERT.hh"
#include "HadronPhysicsQGSC_EFLOW.hh"
#include "HadronPhysicsQGSP.hh"
#include "HadronPhysicsQGSP_BERT.hh"
#include "HadronPhysicsQGSP_BERT_HP.hh"
#include "HadronPhysicsQGSP_BIC.hh"
#include "HadronPhysicsQGSP_BIC_HP.hh"

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
  dump            = false;

  pMessenger = new PhysicsListMessenger(this);

  // Particles
  particleList = new G4DecayPhysics("decays");

  // EM physics
  emPhysicsList = new G4EmStandardPhysics();
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

  if(dump) G4HadronProcessStore::Instance()->Dump(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0)
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;

  if (name == "emstandard_opt2") {

    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option2();

  } else if (name == "emstandard_opt1") {

    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option1();

  } else if (name == "FTFC") {

    SetBuilderList1();
    hadronPhys.push_back( new HadronPhysicsFTFC("hadron",true));
    dump = true;

  } else if (name == "FTFP") {

    SetBuilderList1();
    hadronPhys.push_back( new HadronPhysicsFTFP("hadron",true));
    dump = true;

  } else if (name == "FTFP_BERT") {

    SetBuilderList1();
    hadronPhys.push_back( new HadronPhysicsFTFP_BERT("hadron",true));
    dump = true;

  } else if (name == "FTFP_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("FTFP");

  } else if (name == "LHEP") {

    SetBuilderList2();
    hadronPhys.push_back( new HadronPhysicsLHEP("hadron"));
    dump = true;

  } else if (name == "LHEP_BERT") {

    SetBuilderList3();
    hadronPhys.push_back( new HadronPhysicsLHEP_BERT("hadron"));
    dump = true;

  } else if (name == "LHEP_EMV") {

    AddPhysicsList("emstandard_opt1");
    SetBuilderList3();
    hadronPhys.push_back( new HadronPhysicsLHEP_EMV("hadron"));
    dump = true;

  } else if (name == "LHEP_PRECO_HP") {

    SetBuilderList3(true);
    hadronPhys.push_back( new HadronPhysicsLHEP_PRECO_HP("hadron"));
    dump = true;

  } else if (name == "QBBC") {

    SetBuilderList0();
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,false,false,true));

  } else if (name == "QBBC_DEL") {

    SetBuilderList5();
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,false,false,true));

  } else if (name == "QBBC_HEL") {

    SetBuilderList6();
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,false,false,true));

  } else if (name == "QBBC_HP") {

    SetBuilderList0(true);
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,false,true,true));
  } else if (name == "QGSC") {

    SetBuilderList4();
    hadronPhys.push_back( new HadronPhysicsQGSC("hadron",true));
    dump = true;

  } else if (name == "QGSC_BERT") {

    SetBuilderList4();
    hadronPhys.push_back( new HadronPhysicsQGSC_BERT("hadron",true));
    dump = true;

  } else if (name == "QGSC_EFLOW") {

    SetBuilderList4();
    hadronPhys.push_back( new HadronPhysicsQGSC_EFLOW("hadron",true));
    dump = true;

  } else if (name == "QGSC_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("QGSC");

  } else if (name == "QGSP") {

    SetBuilderList1();
    hadronPhys.push_back( new HadronPhysicsQGSP("hadron",true));
    dump = true;

  } else if (name == "QGSP_BERT") {

    SetBuilderList1();
    hadronPhys.push_back( new HadronPhysicsQGSP_BERT("hadron",true));
    dump = true;

  } else if (name == "QGSP_BERT_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("QGSP_BERT");

  } else if (name == "QGSP_BERT_HP") {

    SetBuilderList1(true);
    hadronPhys.push_back( new HadronPhysicsQGSP_BERT_HP("hadron",true));

  } else if (name == "QGSP_BIC") {

    SetBuilderList0();
    hadronPhys.push_back( new HadronPhysicsQGSP_BIC("hadron",true));
    dump = true;

  } else if (name == "QGSP_BIC_HP") {

    SetBuilderList0(true);
    hadronPhys.push_back( new HadronPhysicsQGSP_BIC_HP("hadron",true));
    dump = true;

  } else if (name == "QGSP_DIF") {

    SetBuilderList1();
    HadronPhysicsQGSP * hp = new HadronPhysicsQGSP("hadron");
    hp->SetQuasiElastic(true);
    hp->SetProjectileDiffraction(true);
    hadronPhys.push_back(hp);
    dump = true;

  } else if (name == "QGSP_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("QGSP");
    dump = true;

  } else if (name == "QGSP_EMX") {

    AddPhysicsList("emstandard_opt2");
    AddPhysicsList("QGSP");
    dump = true;

  } else if (name == "QGSP_NQE") {

    SetBuilderList1();
    hadronPhys.push_back( new HadronPhysicsQGSP("hadron",false));
    dump = true;

  } else if (name == "QGSP_QEL") {

    SetBuilderList4();
    hadronPhys.push_back( new HadronPhysicsQGSP("hadron",true));
    dump = true;

  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList0(G4bool flagHP)
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronElasticPhysics("elastic",verboseLevel,
						    flagHP));
  hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonBinaryCascadePhysics("ionBIC"));
  hadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList1(G4bool flagHP)
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronElasticPhysics("elastic",verboseLevel,
						    flagHP));
  hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonPhysics("ion"));
  hadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList2(G4bool flagHP)
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronElasticPhysics("LElastic",verboseLevel,
						    flagHP));
  hadronPhys.push_back( new G4IonPhysics("ion"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList3(G4bool flagHP)
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronElasticPhysics("LElastic",verboseLevel,
						    flagHP));
  hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonPhysics("ion"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList4(G4bool)
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronQElasticPhysics("elastic",verboseLevel));
  hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonPhysics("ion"));
  hadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList5(G4bool flagHP)
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronDElasticPhysics(verboseLevel,flagHP));
  hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonBinaryCascadePhysics("ionBIC"));
  hadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList6(G4bool flagHP)
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronHElasticPhysics(verboseLevel,flagHP));
  hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonBinaryCascadePhysics("ionBIC"));
  hadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
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
  G4cout << "### PhysicsLists available: FTFC FTFP FTFP_BERT FTFP_EMV LHEP LHEP_BERT LHEP_EMV "
	 << G4endl;
  G4cout << "                            LHEP_PRECO_HP QBBC QBBC_DEL QBBC_HEL QBBC_HP QGSC "
	 << G4endl; 
  G4cout << "                            QGSC_BERT QGSC_EFLOW QGSC_EMV QGSP QGSP_BERT QGSP_BER_EMV "
	 << G4endl; 
  G4cout << "                            QGSP_BERT_HP QGSP_BIC QGSP_BIC_HP QGSP_DIF " 
	 << G4endl; 
  G4cout << "                            QGSP_EMV QGSP_EMX QGSP_NQE QGSP_QEL  "
	 << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

