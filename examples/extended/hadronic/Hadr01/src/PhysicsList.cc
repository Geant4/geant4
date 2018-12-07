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
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4ChargeExchangePhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4NeutronCrossSectionXS.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonPhysics.hh"
#include "G4IonPhysicsXS.hh"
#include "G4IonElasticPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmParameters.hh"
#include "G4PhysListFactoryMessenger.hh"

#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsFTFP_BERT_TRV.hh"
#include "G4HadronPhysicsFTF_BIC.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_FTFP_BERT.hh"
#include "G4HadronPhysicsQGS_BIC.hh"
#include "G4RadioactiveDecayPhysics.hh"

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
  SetDefaultCutValue(0.7*CLHEP::mm);

  verboseLevel = 1;

  fMessenger = new PhysicsListMessenger(this);
  fFactMessenger = new G4PhysListFactoryMessenger(this);

  // Particles
  fParticleList = new G4DecayPhysics(verboseLevel);

  // EM physics
  fEmPhysicsList = new G4EmStandardPhysics(verboseLevel);
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
    fEmPhysicsList = new G4EmStandardPhysics(verboseLevel);

  } else if (name == "emstandard_opt1") {

    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option1(verboseLevel);

  } else if (name == "emstandard_opt2") {

    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option2(verboseLevel);

  } else if (name == "emstandard_opt3") {

    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3(verboseLevel);

  } else if (name == "emstandard_opt4") {

    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4(verboseLevel);

  } else if (name == "emstandardGS") {

    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsGS(verboseLevel);

  } else if (name == "FTFP_BERT_EMV") {

    AddPhysicsList("FTFP_BERT");
    AddPhysicsList("emstandard_opt1");

  } else if (name == "FTFP_BERT_EMX") {

    AddPhysicsList("FTFP_BERT");
    AddPhysicsList("emstandard_opt2");

  } else if (name == "FTFP_BERT_EMY") {

    AddPhysicsList("FTFP_BERT");
    AddPhysicsList("emstandard_opt3");

  } else if (name == "FTFP_BERT_EMZ") {

    AddPhysicsList("FTFP_BERT");
    AddPhysicsList("emstandard_opt4");

  } else if (name == "FTFP_BERT") {

    SetBuilderList0(false);
    fHadronPhys.push_back( new G4HadronPhysicsFTFP_BERT(verboseLevel));

  } else if (name == "FTFP_BERT_TRV") {

    AddPhysicsList("emstandardGS");
    G4EmParameters::Instance()->SetMscStepLimitType( fUseSafety );

    SetBuilderList1(false);
    fHadronPhys.push_back( new G4HadronPhysicsFTFP_BERT_TRV(verboseLevel));

  } else if (name == "FTF_BIC") {

    SetBuilderList0(false);
    fHadronPhys.push_back( new G4HadronPhysicsFTF_BIC(verboseLevel));

  } else if (name == "QBBC") {

    AddPhysicsList("emstandard_opt0");
    SetBuilderList2();
    fHadronPhys.push_back( new G4HadronInelasticQBBC(verboseLevel));

  } else if (name == "QGSP_BERT") {

    SetBuilderList0(false);
    fHadronPhys.push_back( new G4HadronPhysicsQGSP_BERT(verboseLevel));

  } else if (name == "QGSP_FTFP_BERT") {

    SetBuilderList0(false);
    fHadronPhys.push_back( new G4HadronPhysicsQGSP_FTFP_BERT(verboseLevel));

  } else if (name == "QGSP_FTFP_BERT_EMV") {

    AddPhysicsList("QGSP_FTFP_BERT");
    AddPhysicsList("emstandard_opt1");

  } else if (name == "QGSP_BERT_EMV") {

    AddPhysicsList("QGSP_BERT");
    AddPhysicsList("emstandard_opt1");

  } else if (name == "QGSP_BERT_EMX") {

    AddPhysicsList("QGSP_BERT");
    AddPhysicsList("emstandard_opt2");

  } else if (name == "QGSP_BERT_HP") {

    SetBuilderList0(true);
    fHadronPhys.push_back( new G4HadronPhysicsQGSP_BERT_HP(verboseLevel));

  } else if (name == "QGSP_BIC") {

    SetBuilderList0(false);
    fHadronPhys.push_back( new G4HadronPhysicsQGSP_BIC(verboseLevel));
    fHadronPhys.push_back( new G4IonElasticPhysics(verboseLevel));

  } else if (name == "QGSP_BIC_EMY") {

    AddPhysicsList("QGSP_BIC");
    AddPhysicsList("emstandard_opt3");

  } else if (name == "QGS_BIC") {

    SetBuilderList0(false);
    fHadronPhys.push_back( new G4HadronPhysicsQGS_BIC(verboseLevel));

  } else if (name == "QGSP_BIC_HP") {

    SetBuilderList0(true);
    fHadronPhys.push_back( new G4HadronPhysicsQGSP_BIC_HP(verboseLevel));

  } else if (name == "RadioactiveDecay") {

    fHadronPhys.push_back( new G4RadioactiveDecayPhysics(verboseLevel));

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
  fHadronPhys.push_back( new G4IonPhysics(verboseLevel));
  fHadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList1(G4bool flagHP)
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  if(flagHP) {
    fHadronPhys.push_back( new G4HadronElasticPhysicsHP(verboseLevel) );
  } else {
    fHadronPhys.push_back( new G4HadronHElasticPhysics(verboseLevel) );
  }
  fHadronPhys.push_back( new G4StoppingPhysics(verboseLevel));
  fHadronPhys.push_back( new G4IonPhysics(verboseLevel));
  fHadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList2()
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  fHadronPhys.push_back( new G4HadronElasticPhysicsXS(verboseLevel) );
  fHadronPhys.push_back( new G4StoppingPhysics(verboseLevel));
  fHadronPhys.push_back( new G4IonPhysicsXS(verboseLevel));
  fHadronPhys.push_back( new G4IonElasticPhysics(verboseLevel));
  fHadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::List()
{
  G4cout << "### PhysicsLists available: FTFP_BERT FTFP_BERT_EMV "
         << "FTFP_BERT_EMX FTFP_BERT_EMZ FTFP_BERT_TRV"
         << G4endl;
  G4cout << "                            FTF_BIC QBBC QGSP_BERT "
         << "QGSP_BERT_EMV QGSP_BERT_EMX"
         << G4endl; 
  G4cout << "                            QGSP_BERT_HP QGSP_FTFP_BERT "
         << "QGSP_FTFP_BERT_EMV"
         << G4endl; 
  G4cout << "                            QGS_BIC QGSP_BIC QGSP_BIC_EMY "
         << "QGSP_BIC_HP" 
         << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

