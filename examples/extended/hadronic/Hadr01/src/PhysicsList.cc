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
// $Id: PhysicsList.cc,v 1.37 2010/11/19 20:12:32 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-04 $
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
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronElasticPhysicsLHEP.hh"
#include "G4HadronQElasticPhysics.hh"
#include "G4ChargeExchangePhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4NeutronCrossSectionXS.hh"
#include "G4QStoppingPhysics.hh"
#include "G4LHEPStoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmProcessOptions.hh"

#include "HadronPhysicsFTFP_BERT.hh"
#include "HadronPhysicsFTF_BIC.hh"
#include "HadronPhysicsLHEP.hh"
#include "HadronPhysicsLHEP_EMV.hh"
#include "G4HadronInelasticQBBC.hh"
#include "HadronPhysicsQGSC_BERT.hh"
#include "HadronPhysicsQGSP.hh"
#include "HadronPhysicsQGSP_BERT.hh"
#include "HadronPhysicsQGSP_BERT_HP.hh"
#include "HadronPhysicsQGSP_BIC.hh"
#include "HadronPhysicsQGSP_BIC_HP.hh"
#include "HadronPhysicsQGSP_FTFP_BERT.hh"
#include "HadronPhysicsQGS_BIC.hh"

#include "G4IonPhysics.hh"

#include "G4LossTableManager.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 0.7*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  cutForProton    = defaultCutValue;
  verboseLevel    = 1;

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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }
  if (name == "emstandard_opt2") {

    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option2();

  } else if (name == "emstandard_opt3") {

    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option3();

  } else if (name == "emstandard_opt1") {

    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option1();

  } else if (name == "emstandard_opt0") {

    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics();

  } else if (name == "FTFP_BERT_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("FTFP_BERT");

  } else if (name == "FTFP_BERT_EMX") {

    AddPhysicsList("emstandard_opt2");
    AddPhysicsList("FTFP_BERT");

  } else if (name == "FTFP_BERT") {

    SetBuilderList1();
    hadronPhys.push_back( new HadronPhysicsFTFP_BERT());

  } else if (name == "FTF_BIC") {

    SetBuilderList0();
    hadronPhys.push_back( new HadronPhysicsFTF_BIC());
    hadronPhys.push_back( new G4NeutronCrossSectionXS(verboseLevel));

  } else if (name == "LHEP") {

    SetBuilderList2();
    hadronPhys.push_back( new HadronPhysicsLHEP());

  } else if (name == "LHEP_EMV") {

    AddPhysicsList("emstandard_opt1");
    SetBuilderList2(true);
    hadronPhys.push_back( new HadronPhysicsLHEP_EMV());

  } else if (name == "QBBC") {

    AddPhysicsList("emstandard_opt2");
    SetBuilderList3();
    hadronPhys.push_back( new G4HadronInelasticQBBC());

  } else if (name == "QGSC_BERT") {

    SetBuilderList4();
    hadronPhys.push_back( new HadronPhysicsQGSC_BERT());

  } else if (name == "QGSP") {

    SetBuilderList1();
    hadronPhys.push_back( new HadronPhysicsQGSP());

  } else if (name == "QGSP_BERT") {

    SetBuilderList1();
    hadronPhys.push_back( new HadronPhysicsQGSP_BERT());

  } else if (name == "QGSP_FTFP_BERT") {

    SetBuilderList1();
    hadronPhys.push_back( new HadronPhysicsQGSP_FTFP_BERT());

  } else if (name == "QGSP_BERT_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("QGSP_BERT");

  } else if (name == "QGSP_BERT_EMX") {

    AddPhysicsList("emstandard_opt2");
    AddPhysicsList("QGSP_BERT");

  } else if (name == "QGSP_BERT_HP") {

    SetBuilderList1(true);
    hadronPhys.push_back( new HadronPhysicsQGSP_BERT_HP());

  } else if (name == "QGSP_BIC") {

    SetBuilderList0();
    hadronPhys.push_back( new HadronPhysicsQGSP_BIC());

  } else if (name == "QGSP_BIC_EMY") {

    AddPhysicsList("emstandard_opt3");
    SetBuilderList0();
    hadronPhys.push_back( new HadronPhysicsQGSP_BIC());

  } else if (name == "QGS_BIC") {

    SetBuilderList0();
    hadronPhys.push_back( new HadronPhysicsQGS_BIC());
    hadronPhys.push_back( new G4NeutronCrossSectionXS(verboseLevel));

  } else if (name == "QGSP_BIC_HP") {

    SetBuilderList0(true);
    hadronPhys.push_back( new HadronPhysicsQGSP_BIC_HP());

  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList0(G4bool flagHP)
{
  hadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  if(flagHP) {
    hadronPhys.push_back( new G4HadronElasticPhysicsHP(verboseLevel) );
  } else {
    hadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel) );
  }
  hadronPhys.push_back( new G4QStoppingPhysics(verboseLevel));
  hadronPhys.push_back( new G4IonBinaryCascadePhysics(verboseLevel));
  hadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList1(G4bool flagHP)
{
  hadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  if(flagHP) {
    hadronPhys.push_back( new G4HadronElasticPhysicsHP(verboseLevel) );
  } else {
    hadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel) );
  }
  hadronPhys.push_back( new G4QStoppingPhysics(verboseLevel));
  hadronPhys.push_back( new G4IonPhysics(verboseLevel));
  hadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList2(G4bool addStopping)
{
  hadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  hadronPhys.push_back( new G4HadronElasticPhysicsLHEP(verboseLevel));
  if(addStopping) { hadronPhys.push_back( new G4QStoppingPhysics(verboseLevel)); }
  hadronPhys.push_back( new G4IonPhysics(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList3()
{
  hadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  RegisterPhysics( new G4HadronElasticPhysicsXS(verboseLevel) );
  hadronPhys.push_back( new G4QStoppingPhysics(verboseLevel));
  hadronPhys.push_back( new G4IonBinaryCascadePhysics(verboseLevel));
  hadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList4()
{
  hadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  hadronPhys.push_back( new G4HadronQElasticPhysics(verboseLevel));
  hadronPhys.push_back( new G4QStoppingPhysics(verboseLevel));
  hadronPhys.push_back( new G4IonPhysics(verboseLevel));
  hadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
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
  SetCutValue(cutForProton, "proton");

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForProton(G4double cut)
{
  cutForProton = cut;
  SetParticleCuts(cutForProton, G4Proton::Proton());
}

void PhysicsList::List()
{
  G4cout << "### PhysicsLists available: FTFP_BERT FTFP_BERT_EMV FTFP_BERT_EMX FTF_BIC"
	 << G4endl;
  G4cout << "                            LHEP LHEP_EMV QBBC QGS_BIC QGSP"
	 << G4endl; 
  G4cout << "                            QGSC_BERT QGSP_BERT QGSP_BERT_EMV QGSP_BIC_EMY"
	 << G4endl; 
  G4cout << "                            QGSP_BERT_EMX QGSP_BERT_HP QGSP_BIC QGSP_BIC_HP" 
	 << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

