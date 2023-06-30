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
/// \file XrayTESdetPhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "G4VUserPhysicsList.hh"
#include "XrayTESdetPhysicsList.hh"
#include "XrayTESdetPhysicsListMessenger.hh"

#include "G4RadioactiveDecayPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsSS.hh"

#include "G4EmStandardPhysics_SpacePhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4ChargeExchangePhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4NeutronCrossSectionXS.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonPhysics.hh"
#include "G4EmExtraPhysics.hh"

#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsFTF_BIC.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_FTFP_BERT.hh"
#include "G4HadronPhysicsQGS_BIC.hh"

#include "G4LossTableManager.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

XrayTESdetPhysicsList::XrayTESdetPhysicsList() : fMessenger(nullptr), fEmPhysicsList(nullptr), fParticleList(nullptr)
{
  fMessenger = new XrayTESdetPhysicsListMessenger(this);
  G4LossTableManager::Instance();
  verboseLevel = 1;

  // Particles
  G4cout << "1 - Defining DecayPhysics" << G4endl;
  fParticleList = new G4DecayPhysics("decays");

  // Radioactive decay
  G4cout << "2 - Defining RadioactiveDecayPhysics" << G4endl;
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  // EM physics
  G4cout << "3 - Defining Standard em" << G4endl;
  fEmName = G4String("local");
  //fEmPhysicsList = new PhysListEmStandard(fEmName);
  fEmPhysicsList = new G4EmStandardPhysics();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

XrayTESdetPhysicsList::~XrayTESdetPhysicsList()
{
  delete fParticleList;
  delete fEmPhysicsList;
  delete fMessenger;
  for(size_t i=0; i<fHadronPhys.size(); i++)
  {
    delete fHadronPhys[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void XrayTESdetPhysicsList::ConstructParticle()
{
  fParticleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void XrayTESdetPhysicsList::ConstructProcess()
{
  AddTransportation();
  fEmPhysicsList->ConstructProcess();
  fParticleList->ConstructProcess();
  for(size_t i=0; i<fHadronPhys.size(); i++)
  {
    fHadronPhys[i]->ConstructProcess();
  }

  // Em options
  G4cout << "4 - Defining em options" << G4endl;
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDeexActiveRegion("InnerRegion", true, true, true);

  param->SetAuger(true);
  param->SetAugerCascade(false);
  param->SetFluo(true);
  param->SetPixe(true);
  param->SetDeexcitationIgnoreCut(false);

  param->SetMuHadLateralDisplacement(false);
  param->SetBremsstrahlungTh(10*TeV);
  param->Dump();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void XrayTESdetPhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0)
  {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == fEmName)
  {
    return;
  }

  if (name == "SpacePhysics")
  {
    delete fEmPhysicsList;
    G4cout << "X - Defining SpacePhysics, no hadronic part" << G4endl;
    fEmName = name;
    fEmPhysicsList = new G4EmStandardPhysics_SpacePhysics();
  }
  else if (name == "SpacePhysics_QBBC")
  {
    delete fEmPhysicsList;
    G4cout << "X - Defining SpacePhysics" << G4endl;
    fEmName = name;
    fEmPhysicsList = new G4EmStandardPhysics_SpacePhysics();
    fSetBuilderList1(false);
    fHadronPhys.push_back(new G4HadronInelasticQBBC(verboseLevel));
  }
  else if (name == "opt4_QGSP")
  {
    delete fEmPhysicsList;
    fEmName = name;
    G4cout << "X - Defining opt4+QGSP no HP physics" << G4endl;
    fEmPhysicsList = new G4EmStandardPhysics_option4();
    fSetBuilderList1(false);
    fHadronPhys.push_back( new G4HadronPhysicsQGSP_BERT_HP());
  }
  else if (name == "opt4_QBBC")
  {
    delete fEmPhysicsList;
    fEmName = name;
    G4cout << "X - Defining opt4+QBBC physics" << G4endl;
    fEmPhysicsList = new G4EmStandardPhysics_option4();
    fSetBuilderList1(false);
    fHadronPhys.push_back(new G4HadronInelasticQBBC(verboseLevel));
  }
  else if (name == "emstandard_opt4")
  {
    delete fEmPhysicsList;
    G4cout << "X - Defining opt4" << G4endl;
    fEmName = name;
    fEmPhysicsList  = new G4EmStandardPhysics_option4();
  }
  else if (name == "emstandardSS")
  {
    delete fEmPhysicsList;
    G4cout << "X - Defining SS" << G4endl;
    fEmName = name;
    fEmPhysicsList  = new G4EmStandardPhysicsSS();
  }
  else
  {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void XrayTESdetPhysicsList::fSetBuilderList1(G4bool flagHP)
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  if(flagHP)
  {
    fHadronPhys.push_back( new G4HadronElasticPhysicsHP(verboseLevel));
  }
  else
  {
    fHadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel));
  }
  fHadronPhys.push_back(new G4StoppingPhysics(verboseLevel));
  fHadronPhys.push_back(new G4IonPhysics(verboseLevel));
  fHadronPhys.push_back(new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
