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
// $Id: PhysicsList.cc,v 1.21 2007/06/01 16:07:31 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-00 $
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
#include "G4HadronQElasticPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4QStoppingPhysics.hh"
#include "G4LHEPStoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmProcessOptions.hh"

#include "HadronPhysicsFTFP.hh"
#include "HadronPhysicsFTFC.hh"
#include "HadronPhysicsLHEP.hh"
#include "HadronPhysicsLHEP_BERT.hh"
#include "HadronPhysicsLHEP_EMV.hh"
#include "HadronPhysicsLHEP_PRECO_HP.hh"
#include "G4HadronInelasticQBBC.hh"
#include "HadronPhysicsQGSC.hh"
#include "HadronPhysicsQGSC_EFLOW.hh"
#include "HadronPhysicsQGSP.hh"
#include "HadronPhysicsQGSP_BERT.hh"
#include "HadronPhysicsQGSP_BERT_HP.hh"
#include "HadronPhysicsQGSP_BERT_TRV.hh"
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

    hadronPhys.push_back( new HadronPhysicsFTFC("hadron",true));
    SetStandardList(false, false);
    dump = true;

  } else if (name == "FTFP") {

    hadronPhys.push_back( new HadronPhysicsFTFP("hadron",true));
    SetStandardList(false, false);
    dump = true;

  } else if (name == "FTFP_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("FTFP");
    dump = true;

  } else if (name == "LHEP") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("LElastic",
						     verboseLevel,false));
    hadronPhys.push_back( new HadronPhysicsLHEP("hadron"));
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else if (name == "LHEP_BERT") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("LElastic",
						     verboseLevel,false));
    hadronPhys.push_back( new HadronPhysicsLHEP_BERT("hadron"));
    hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else if (name == "LHEP_EMV") {

    AddPhysicsList("emstandard_opt1");
    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("LElastic",
						     verboseLevel,false));
    hadronPhys.push_back( new HadronPhysicsLHEP_EMV("hadron"));
    hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else if (name == "LHEP_PRECO_HP") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronElasticPhysics("LElastic",
						     verboseLevel,false));
    hadronPhys.push_back( new HadronPhysicsLHEP_PRECO_HP("hadron"));
    hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
    hadronPhys.push_back( new G4IonPhysics("ion"));

  } else if (name == "QBBC") {

    SetStandardList(false,true);
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,false,false,true));

  } else if (name == "QBBCG") {

    SetStandardList(false, false);
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,false,false,false));

  } else if (name == "QBEC") {

    SetStandardList(false,true);
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,true,false,false,true));

  } else if (name == "QBBC_HP") {

    SetStandardList(true,true);
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,false,true,true));

  } else if (name == "QBEC_HP") {

    SetStandardList(true,true);
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,true,false,true,true));

  } else if (name == "QGSC") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronQElasticPhysics("QElastic",verboseLevel));
    hadronPhys.push_back( new HadronPhysicsQGSC("hadron",true));
    hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel,false));
    hadronPhys.push_back( new G4IonPhysics("ion"));
    hadronPhys.push_back( new G4NeutronTrackingCut());
    dump = true;

  } else if (name == "QGSC_EFLOW") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronQElasticPhysics("QElastic",verboseLevel));
    hadronPhys.push_back( new HadronPhysicsQGSC_EFLOW("hadron",true));
    hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel,false));
    hadronPhys.push_back( new G4IonPhysics("ion"));
    hadronPhys.push_back( new G4NeutronTrackingCut());
    dump = true;

  } else if (name == "QGSC_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("QGSC");
    dump = true;

  } else if (name == "QGSP") {

    SetStandardList(false, false);
    hadronPhys.push_back( new HadronPhysicsQGSP("hadron",true));
    dump = true;

  } else if (name == "QGSP_BERT") {

    SetStandardList(false, false);
    hadronPhys.push_back( new HadronPhysicsQGSP_BERT("hadron",true));
    dump = true;

  } else if (name == "QGSP_BERT_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("QGSP_BERT");
    dump = true;

  } else if (name == "QGSP_BERT_HP") {

    SetStandardList(true, false);
    hadronPhys.push_back( new HadronPhysicsQGSP_BERT_HP("hadron",true));
    dump = true;

  } else if (name == "QGSP_BERT_NQE") {

    SetStandardList(false, false);
    hadronPhys.push_back( new HadronPhysicsQGSP_BERT("hadron",false));
    dump = true;

  } else if (name == "QGSP_BERT_TRV") {

    SetStandardList(false, false);
    hadronPhys.push_back( new HadronPhysicsQGSP_BERT_TRV("hadron",true));
    dump = true;

  } else if (name == "QGSP_BIC") {

    SetStandardList(false, false);
    hadronPhys.push_back( new HadronPhysicsQGSP_BIC("hadron",true));
    dump = true;

  } else if (name == "QGSP_BIC_HP") {

    SetStandardList(true, false);
    hadronPhys.push_back( new HadronPhysicsQGSP_BIC_HP("hadron",true));
    dump = true;

  } else if (name == "QGSP_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("QGSP");
    dump = true;

  } else if (name == "QGSP_EMV_NQE") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("QGSP_NQE");
    dump = true;

  } else if (name == "QGSP_EMX") {

    AddPhysicsList("emstandard_opt2");
    AddPhysicsList("QGSP");
    dump = true;

  } else if (name == "QGSP_NQE") {

    SetStandardList(false, false);
    hadronPhys.push_back( new HadronPhysicsQGSP("hadron",false));
    dump = true;

  } else if (name == "QGSP_QEL") {

    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronQElasticPhysics("QElastic",
						      verboseLevel));
    hadronPhys.push_back( new HadronPhysicsQGSP("hadron",true));
    hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
    hadronPhys.push_back( new G4IonPhysics("ion"));
    hadronPhys.push_back( new G4NeutronTrackingCut());
    dump = true;

  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetStandardList(G4bool flagHP, G4bool glauber)
{
  hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
  hadronPhys.push_back( new G4HadronElasticPhysics("elastic",verboseLevel,
						    flagHP,glauber));
  hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonBinaryCascadePhysics("binary_ion"));
  hadronPhys.push_back( new G4NeutronTrackingCut());
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
  G4cout << "### PhysicsLists available: FTFC FTFP FTFP_EMV LHEP LHEP_BERT LHEP_EMV"
	 << G4endl;
  G4cout << "                            LHEP_PRECO_HP QBBC QBBCG QBEC QBBC_HP QBEC_HP"
	 << G4endl; 
  G4cout << "                            QGSC QGSC_EFLOW QGSC_EMV QGSP QGSP_BERT QGSP_BER_EMV"
	 << G4endl; 
  G4cout << "                            QGSP_BERT_HP QGSP_BERT_NQE QGSP_BERT_TRV"
	 << G4endl; 
  G4cout << "                            QGSP_BIC QGSP_BIC_HP QGSP_EMV QGSP_EMV_NQE" 
	 << G4endl; 
  G4cout << "                            QGSC_EMX QGSP_NQE QGSP_QEL"
	 << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

