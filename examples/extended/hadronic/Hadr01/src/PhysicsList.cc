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
/// \file hadronic/Hadr01/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
//
// $Id$
//
/////////////////////////////////////////////////////////////////////////
//
// PhysicsList
//
// Created: 31.04.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
// 26.04.2007 Physics according to 8.3 Physics List (V.Ivanchenko)
// 16.10.2012 Renamed used classes (A.Ribon)
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
#include "G4EmStandardPhysics_option4.hh"
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
#include "G4StoppingPhysics.hh"
#include "G4LHEPStoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonLHEPPhysics.hh"
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

#include "G4LossTableManager.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 0.7*mm;
  fCutForGamma     = defaultCutValue;
  fCutForElectron  = defaultCutValue;
  fCutForPositron  = defaultCutValue;
  fCutForProton    = defaultCutValue;
  verboseLevel    = 1;

  fMessenger = new PhysicsListMessenger(this);

  // Particles
  fParticleList = new G4DecayPhysics("decays");

  // EM physics
  fEmPhysicsList = new G4EmStandardPhysics();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PhysicsList::~PhysicsList()
{
  delete fMessenger;
  delete fParticleList;
  delete fEmPhysicsList;
  for(size_t i=0; i<fHadronPhys.size(); i++) {
    delete fHadronPhys[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::ConstructParticle()
{
  fParticleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  fEmPhysicsList->ConstructProcess();
  fParticleList->ConstructProcess();
  for(size_t i=0; i<fHadronPhys.size(); i++) {
    fHadronPhys[i]->ConstructProcess();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }
  if (name == "emstandard_opt0") {

    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics();

  } else if (name == "emstandard_opt1") {

    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option1();

  } else if (name == "emstandard_opt2") {

    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option2();

  } else if (name == "emstandard_opt3") {

    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3();

  } else if (name == "emstandard_opt4") {

    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4();

  } else if (name == "FTFP_BERT_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("FTFP_BERT");

  } else if (name == "FTFP_BERT_EMX") {

    AddPhysicsList("emstandard_opt2");
    AddPhysicsList("FTFP_BERT");

  } else if (name == "FTFP_BERT") {

    SetBuilderList1();
    fHadronPhys.push_back( new HadronPhysicsFTFP_BERT());

  } else if (name == "FTF_BIC") {

    SetBuilderList0();
    fHadronPhys.push_back( new HadronPhysicsFTF_BIC());
    fHadronPhys.push_back( new G4NeutronCrossSectionXS(verboseLevel));

  } else if (name == "LHEP") {

    SetBuilderList2();
    fHadronPhys.push_back( new HadronPhysicsLHEP());

  } else if (name == "LHEP_EMV") {

    AddPhysicsList("emstandard_opt1");
    SetBuilderList2(true);
    fHadronPhys.push_back( new HadronPhysicsLHEP_EMV());

  } else if (name == "QBBC") {

    AddPhysicsList("emstandard_opt0");
    SetBuilderList3();
    fHadronPhys.push_back( new G4HadronInelasticQBBC());

  } else if (name == "QGSC_BERT") {

    SetBuilderList4();
    fHadronPhys.push_back( new HadronPhysicsQGSC_BERT());

  } else if (name == "QGSP") {

    SetBuilderList1();
    fHadronPhys.push_back( new HadronPhysicsQGSP());

  } else if (name == "QGSP_BERT") {

    SetBuilderList1();
    fHadronPhys.push_back( new HadronPhysicsQGSP_BERT());

  } else if (name == "QGSP_FTFP_BERT") {

    SetBuilderList1();
    fHadronPhys.push_back( new HadronPhysicsQGSP_FTFP_BERT());

  } else if (name == "QGSP_BERT_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("QGSP_BERT");

  } else if (name == "QGSP_BERT_EMX") {

    AddPhysicsList("emstandard_opt2");
    AddPhysicsList("QGSP_BERT");

  } else if (name == "QGSP_BERT_HP") {

    SetBuilderList1(true);
    fHadronPhys.push_back( new HadronPhysicsQGSP_BERT_HP());

  } else if (name == "QGSP_BIC") {

    SetBuilderList0();
    fHadronPhys.push_back( new HadronPhysicsQGSP_BIC());

  } else if (name == "QGSP_BIC_EMY") {

    AddPhysicsList("emstandard_opt3");
    SetBuilderList0();
    fHadronPhys.push_back( new HadronPhysicsQGSP_BIC());

  } else if (name == "QGS_BIC") {

    SetBuilderList0();
    fHadronPhys.push_back( new HadronPhysicsQGS_BIC());
    fHadronPhys.push_back( new G4NeutronCrossSectionXS(verboseLevel));

  } else if (name == "QGSP_BIC_HP") {

    SetBuilderList0(true);
    fHadronPhys.push_back( new HadronPhysicsQGSP_BIC_HP());

  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList0(G4bool flagHP)
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  if(flagHP) {
    fHadronPhys.push_back( new G4HadronElasticPhysicsHP(verboseLevel) );
  } else {
    fHadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel) );
  }
  fHadronPhys.push_back( new G4StoppingPhysics(verboseLevel));
  fHadronPhys.push_back( new G4IonBinaryCascadePhysics(verboseLevel));
  fHadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList1(G4bool flagHP)
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  if(flagHP) {
    fHadronPhys.push_back( new G4HadronElasticPhysicsHP(verboseLevel) );
  } else {
    fHadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel) );
  }
  fHadronPhys.push_back( new G4StoppingPhysics(verboseLevel));
  fHadronPhys.push_back( new G4IonLHEPPhysics(verboseLevel));
  fHadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList2(G4bool addStopping)
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  fHadronPhys.push_back( new G4HadronElasticPhysicsLHEP(verboseLevel));
  if(addStopping) { fHadronPhys.push_back( new G4StoppingPhysics(verboseLevel)); }
  fHadronPhys.push_back( new G4IonLHEPPhysics(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList3()
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  fHadronPhys.push_back( new G4HadronElasticPhysicsXS(verboseLevel) );
  fHadronPhys.push_back( new G4StoppingPhysics(verboseLevel));
  fHadronPhys.push_back( new G4IonPhysics(verboseLevel));
  fHadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList4()
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  fHadronPhys.push_back( new G4HadronQElasticPhysics(verboseLevel));
  fHadronPhys.push_back( new G4QStoppingPhysics(verboseLevel));  //16-Oct-2012 A.R. Leave CHIPS stopping
  fHadronPhys.push_back( new G4IonLHEPPhysics(verboseLevel));
  fHadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
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
  SetCutValue(fCutForGamma, "gamma");
  SetCutValue(fCutForElectron, "e-");
  SetCutValue(fCutForPositron, "e+");
  SetCutValue(fCutForProton, "proton");

  if (verboseLevel>0) { DumpCutValuesTable(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetCutForGamma(G4double cut)
{
  fCutForGamma = cut;
  SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
  fCutForElectron = cut;
  SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut)
{
  fCutForPositron = cut;
  SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForProton(G4double cut)
{
  fCutForProton = cut;
  SetParticleCuts(fCutForProton, G4Proton::Proton());
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

