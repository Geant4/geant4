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
//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
//
// ClassName:   test31PhysicsList
//  
// Description: test31 PhysicsList for Geant4 tests 
//
// Authors:   07.04.01  V.Ivanchenko 
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//

#include "test31PhysicsList.hh"
#include "test31DetectorConstruction.hh"
#include "test31PhysicsListMessenger.hh"
#include "test31LowEPhysicsList.hh"
#include "test31StanPhysicsList.hh"
#include "test31StanMAPhysicsList.hh"
#include "test31HadronPhysicsList.hh"
#include "test31HadronPhysicsList1.hh"
#include "test31StepCut.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4IonC12.hh"
#include "G4IonAr40.hh"
#include "G4IonFe56.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4EnergyLossTables.hh"
#include "G4ios.hh"
#include "g4std/iomanip"                
#include "G4Decay.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31PhysicsList::test31PhysicsList(const test31DetectorConstruction* p):
  pDet(p),
  physicsIsDefined(false)
{
  InitializeMe(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::InitializeMe()
{
  verbose = pDet->GetVerbose();
 
  // default cuts
  cutForGamma     = 1.0*cm;
  cutForElectron  = 1.0*mm;
  cutForProton    = 1.0*mm;
  maxChargedStep  = DBL_MAX; 
  lowEnergyLimit  = 10.0*eV;
  highEnergyLimit  = 100.0*GeV;

  theMessenger = new test31PhysicsListMessenger(this);

  emPhysics = G4String("");
  hadronPhysics = G4String("");
  decayPhysics = G4String("none");
  SetEMPhysicsList(G4String("LowEnergy"));  
  SetHadronPhysicsList(G4String("none"));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31PhysicsList::~test31PhysicsList()
{
  delete theMessenger; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructMyBosons();
  ConstructMyLeptons();
  ConstructMyBarions();
  ConstructMyMesons();
  ConstructMyIons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::ConstructMyBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  
  // gamma
  G4Gamma::GammaDefinition();
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::ConstructMyLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::ConstructMyMesons()
{
 //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::ConstructMyBarions()
{
//  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
  G4Lambda::LambdaDefinition();
  G4SigmaZero::SigmaZeroDefinition();
  G4SigmaPlus::SigmaPlusDefinition();
  G4SigmaMinus::SigmaMinusDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::ConstructMyIons()
{
//  Ions
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4Alpha::AlphaDefinition();
  G4IonC12::IonC12Definition();
  G4IonAr40::IonAr40Definition();
  G4IonFe56::IonFe56Definition();
  G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::ConstructProcess()
{
  verbose = pDet->GetVerbose();
  AddTransportation();
  if(theEMList)  theEMList->ConstructEM();
  if(theHadList) theHadList->ConstructHad();
  ConstructDecay();
  physicsIsDefined = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::ConstructDecay()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  test31StepCut* theStepCut = new test31StepCut(G4String("user cut"),this);

  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (theDecayProcess->IsApplicable(*particle) && decayPhysics != "none") { 
      pmanager ->AddProcess(theDecayProcess);

      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
    if (particle->GetPDGCharge()) { 
      pmanager ->AddProcess(theStepCut);
      pmanager ->SetProcessOrdering(theStepCut, idxPostStep);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::SetCuts()
{

  //special for low energy physics
  G4Gamma   ::SetEnergyRange(lowEnergyLimit,highEnergyLimit);
  G4Electron::SetEnergyRange(lowEnergyLimit,highEnergyLimit);
  G4Positron::SetEnergyRange(lowEnergyLimit,highEnergyLimit);
   
  if (verbose > 0){
    G4cout << "test31PhysicsList::SetCuts: "
           << "CutLength = " << maxChargedStep/mm << " mm" << G4endl;
  }  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 

   G4cout << "test31PhysicsList: Set cuts for all particles " << G4endl; 

   SetCutValue(cutForGamma,"gamma");

   SetCutValue(cutForElectron,"e-");
   SetCutValue(cutForElectron,"e+");

   SetCutValue(cutForProton,"mu-");
   SetCutValue(cutForProton,"mu+");

  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton 

  SetCutValue(cutForProton, "proton");

  SetCutValue(cutForProton, "anti_proton");

  SetCutValueForOthers(cutForProton);

  if (verbose > 1) {
    G4cout << "test31PhysicsList: Dump the table" << G4endl;
  }
  DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::SetGammaCut(G4double val)
{
  if(physicsIsDefined) ResetCuts();
  cutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::SetElectronCut(G4double val)
{
  if(physicsIsDefined) ResetCuts();
  cutForElectron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::SetProtonCut(G4double val)
{
  if(physicsIsDefined) ResetCuts();
  cutForProton = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::SetElectronCutByEnergy(G4double energy)
{
  G4Material* currMat = pDet->GetAbsorberMaterial();
  G4ParticleDefinition* part = G4Electron::ElectronDefinition();

  // Get range from electron energy and set it as a cut
  SetElectronCut(G4EnergyLossTables::GetRange(part,energy,currMat));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::SetLowEnergyLimit(G4double energy)
{
  if(physicsIsDefined) ResetCuts();
  lowEnergyLimit = energy;
  G4cout << "test31PhysicsList: lowEnergyLimit = " 
         << lowEnergyLimit/eV << " eV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::SetHighEnergyLimit(G4double energy)
{
  if(physicsIsDefined) ResetCuts();
  highEnergyLimit = energy;
  G4cout << "test31PhysicsList: highEnergyLimit = " 
         << highEnergyLimit/GeV << " GeV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::SetMaxStep(G4double step)
{
  maxChargedStep = step;
  G4cout << "test31PhysicsList: MaxChargedStep = " 
         << maxChargedStep/mm << " mm" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::SetEMPhysicsList(const G4String& name)
{
  verbose = pDet->GetVerbose();
  // Name is not changed
  if(name == emPhysics) return;

  // Define hadronic process class

  if(name == "Standard") {
    theEMList = new test31StanPhysicsList();
    theEMList->SetVerbose(verbose);

  } else if (name == "MAstandard") {
    theEMList = new test31StanMAPhysicsList();
    theEMList->SetVerbose(verbose);

  } else if (name == "LowEnergy") {
    theEMList = new test31LowEPhysicsList();
    theEMList->SetVerbose(verbose);

  } else if (name == "none") {
    theEMList = 0;

  } else {
    G4cout << "HsPhysicsList: There are no EMPhysicsList called <"
           << name << ">, so old one is used" << G4endl;
    return;
  }

  emPhysics = name;   
  G4cout << "HsPhysicsList: <" << name 
         << "> EMPhysicsList is set" << G4endl;
  if(physicsIsDefined) ResetCuts(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31PhysicsList::SetHadronPhysicsList(const G4String& name)
{
  verbose = pDet->GetVerbose();
  // Name is not changed
  if(name == hadronPhysics) return;

  // Define hadronic process class
  G4String stand = G4String("Standard");
  G4String hado  = G4String("LowEnergy");
  G4String frag  = G4String("IonFragmentation");
  G4String none  = G4String("none");

  if(name == "LEparametrised") {
    theHadList = new test31HadronPhysicsList();
    theHadList->SetVerbose(verbose);

  } else if (name == "CHIPS") {
    theHadList = new test31HadronPhysicsList1();
    theHadList->SetVerbose(verbose);

  } else if (name == "none") {
    theHadList = 0;

  } else {
    G4cout << "HsPhysicsList: There are no HadronicPhysicsList called <"
           << name << ">, so old one is used" << G4endl;
    return;
  }

  hadronPhysics = name;   
  G4cout << "HsPhysicsList: <" << name 
         << "> HadronicPhysicsList is set" << G4endl;
  if(physicsIsDefined) ResetCuts(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






