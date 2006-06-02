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
//
// $Id: PhysicsList.cc,v 1.2 2006-06-02 19:06:40 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4QStoppingPhysics.hh"
#include "G4LHEPStoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4LHEPIonPhysics.hh"
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
#include "G4HadronProcessStore.hh"

#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  stepMaxProcess  = 0;

  pMessenger = new PhysicsListMessenger(this);

  SetVerboseLevel(1);

  // Particles
  particleList = new G4DecayPhysics("decays");

  // EM physics
  emName = G4String("standard");
  emPhysicsList = new G4EmStandardPhysics(emName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete pMessenger;
  delete particleList;
  delete emPhysicsList;
  for(size_t i=0; i<hadronPhys.size(); i++) {
    delete hadronPhys[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  particleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  emPhysicsList->ConstructProcess();
  particleList->ConstructProcess();
  for(size_t i=0; i<hadronPhys.size(); i++) {
    hadronPhys[i]->ConstructProcess();
  }

  // Define energy interval for loss processes
  G4EmProcessOptions emOptions;
  emOptions.SetMinEnergy(0.1*keV);
  emOptions.SetMaxEnergy(100.*GeV);
  emOptions.SetDEDXBinning(180);
  emOptions.SetLambdaBinning(180);
  G4HadronProcessStore::Instance()->Dump(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0)
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;

  if (name == "elastic" || name == "elasticB") {

    hadronPhys.push_back( new G4HadronElasticPhysics(name));

  } else if (name == "binary") {

    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",1));

  } else if (name == "binary_ion") {

    hadronPhys.push_back( new G4IonBinaryCascadePhysics(name));

  } else if (name == "gamma_nuc") {

    hadronPhys.push_back( new G4EmExtraPhysics(name));

  } else if (name == "stopping") {

    hadronPhys.push_back( new G4QStoppingPhysics(name));

  } else if (name == "LHEP") {

    hadronPhys.push_back( new HadronPhysicsLHEP());

  } else if (name == "LHEP_BERT") {

    hadronPhys.push_back( new HadronPhysicsLHEP_BERT());

  } else if (name == "LHEP_BERT_HP") {

    hadronPhys.push_back( new HadronPhysicsLHEP_BERT_HP());

  } else if (name == "LHEP_BIC") {

    hadronPhys.push_back( new HadronPhysicsLHEP_BIC());

  } else if (name == "LHEP_BIC_HP") {

    hadronPhys.push_back( new HadronPhysicsLHEP_BIC_HP());

  } else if (name == "QGSC") {

    hadronPhys.push_back( new HadronPhysicsQGSC());

  } else if (name == "QGSP") {

    hadronPhys.push_back( new HadronPhysicsQGSP());

  } else if (name == "QGSP_BERT") {

    hadronPhys.push_back( new HadronPhysicsQGSP_BERT());

  } else if (name == "QGSP_BERT_HP") {

    hadronPhys.push_back( new HadronPhysicsQGSP_BERT_HP());

  } else if (name == "QGSP_BIC") {

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

    hadronPhys.push_back( new G4HadronElasticPhysics("LElastic",
						     verboseLevel,false));
    hadronPhys.push_back( new G4LHEPIonPhysics("lhep_ion"));
    hadronPhys.push_back( new G4LHEPStoppingPhysics("stopping"));
    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronInelasticQLHEP("QGSP",verboseLevel,
						 true,false,false,false));

  } else if (name == "QQGSP_BERT") {

    hadronPhys.push_back( new G4HadronElasticPhysics("LElastic",
						     verboseLevel,false));
    hadronPhys.push_back( new G4LHEPIonPhysics("lhep_ion"));
    hadronPhys.push_back( new G4LHEPStoppingPhysics("stopping"));
    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronInelasticQLHEP("QGSP_BERT",verboseLevel,
						 true,true,false,false));

  } else if (name == "QQGSP_EL") {

    hadronPhys.push_back( new G4HadronElasticPhysics("elastic",
						     verboseLevel,false));
    hadronPhys.push_back( new G4LHEPIonPhysics("lhep_ion"));
    hadronPhys.push_back( new G4LHEPStoppingPhysics("stopping"));
    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronInelasticQLHEP("QGSP",verboseLevel,
						 true,false,false,false));
  } else if (name == "QQGSP_NEL") {

    hadronPhys.push_back( new G4LHEPIonPhysics("lhep_ion"));
    hadronPhys.push_back( new G4LHEPStoppingPhysics("stopping"));
    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
    hadronPhys.push_back( new G4HadronInelasticQLHEP("QGSP",verboseLevel,
						 true,false,false,false));

  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetStandardList(G4bool flagHP)
{
  hadronPhys.push_back( new G4HadronElasticPhysics("elastic",verboseLevel,flagHP));
  hadronPhys.push_back( new G4IonBinaryCascadePhysics("binary_ion"));
  hadronPhys.push_back( new G4QStoppingPhysics("stopping"));
  hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

