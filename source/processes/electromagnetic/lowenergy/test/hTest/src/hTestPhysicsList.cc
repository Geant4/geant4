//---------------------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//---------------------------------------------------------------------------
//
// ClassName:   hTestPhysicsList
//  
// Description: hTest PhysicsList for Geant4 tests 
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

#include "hTestPhysicsList.hh"
#include "hTestDetectorConstruction.hh"
#include "hTestPhysicsListMessenger.hh"
#include "hTestLowEPhysicsList.hh"
#include "hTestStanPhysicsList.hh"
#include "hTestHadronPhysicsList.hh"
#include "hTestHadronPhysicsList1.hh"

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

hTestPhysicsList::hTestPhysicsList(hTestDetectorConstruction* p):
  pDet(p),
  physicsIsDefined(false)
{
  InitializeMe(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::InitializeMe()
{
  verbose = pDet->GetVerbose();
 
  // default cuts
  cutForGamma     = 1.0*cm;
  cutForElectron  = 1.0*mm;
  cutForProton    = 1.0*mm;
  maxChargedStep  = DBL_MAX; 
  lowEnergyLimit  = 250.0*eV;
  highEnergyLimit  = 10.0*GeV;

  theMessenger = new hTestPhysicsListMessenger(this);

  emPhysics = G4String("");
  hadronPhysics = G4String("");
  SetEMPhysicsList(G4String("LowEnergy"));  
  SetHadronPhysicsList(G4String("none"));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPhysicsList::~hTestPhysicsList()
{
  delete theMessenger; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructBarions();
  ConstructMesons();
  ConstructIons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  
  // gamma
  G4Gamma::GammaDefinition();
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructLeptons()
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

void hTestPhysicsList::ConstructMesons()
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

void hTestPhysicsList::ConstructBarions()
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

void hTestPhysicsList::ConstructIons()
{
//  Ions
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4Alpha::AlphaDefinition();
  G4IonC12::IonC12Definition();
  G4IonAr40::IonAr40Definition();
  G4IonFe56::IonFe56Definition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructProcess()
{
  verbose = pDet->GetVerbose();
  AddTransportation();
  if(theEMList)  theEMList->ConstructEM();
  if(theHadList) theHadList->ConstructHad();
  ConstructDecay();
  physicsIsDefined = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructDecay()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);

      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetCuts()
{

  //special for low energy physics
  G4Gamma   ::SetEnergyRange(lowEnergyLimit,highEnergyLimit);
  G4Electron::SetEnergyRange(lowEnergyLimit,highEnergyLimit);
  G4Positron::SetEnergyRange(lowEnergyLimit,highEnergyLimit);
   
  if (verbose > 0){
    G4cout << "hTestPhysicsList::SetCuts: "
           << "CutLength = " << maxChargedStep/mm << " mm" << G4endl;
  }  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 

   G4cout << "hTestPhysicsList: Set cuts for all particles " << G4endl; 

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
    G4cout << "hTestPhysicsList: Dump the table" << G4endl;
  }
  DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetGammaCut(G4double val)
{
  if(physicsIsDefined) ResetCuts();
  cutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetElectronCut(G4double val)
{
  if(physicsIsDefined) ResetCuts();
  cutForElectron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetProtonCut(G4double val)
{
  if(physicsIsDefined) ResetCuts();
  cutForProton = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetElectronCutByEnergy(G4double energy)
{
  G4ParticleTable* thehTestParticleTable =  G4ParticleTable::GetParticleTable();
  G4Material* currMat = pDet->GetAbsorberMaterial();
  G4ParticleDefinition* part = G4Electron::ElectronDefinition();

  // Get range from electron energy and set it as a cut
  SetElectronCut(G4EnergyLossTables::GetRange(part,energy,currMat));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetLowEnergyLimit(G4double energy)
{
  if(physicsIsDefined) ResetCuts();
  lowEnergyLimit = energy;
  G4cout << "hTestPhysicsList: lowEnergyLimit = " 
         << lowEnergyLimit/eV << " eV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetHighEnergyLimit(G4double energy)
{
  if(physicsIsDefined) ResetCuts();
  highEnergyLimit = energy;
  G4cout << "hTestPhysicsList: highEnergyLimit = " 
         << highEnergyLimit/GeV << " GeV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetMaxStep(G4double step)
{
  maxChargedStep = step;
  G4cout << "hTestPhysicsList: MaxChargedStep = " 
         << maxChargedStep/mm << " mm" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetEMPhysicsList(const G4String& name)
{
  verbose = pDet->GetVerbose();
  // Name is not changed
  if(name == emPhysics) return;

  // Define hadronic process class
  G4String stand = G4String("Standard");
  G4String lowe  = G4String("LowEnergy");
  G4String none  = G4String("none");

  if(name == stand) {
    theEMList = new hTestStanPhysicsList();
    theEMList->SetVerbose(verbose);

  } else if (name == lowe) {
    theEMList = new hTestLowEPhysicsList();
    theEMList->SetVerbose(verbose);

  } else if (name == none) {
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

void hTestPhysicsList::SetHadronPhysicsList(const G4String& name)
{
  verbose = pDet->GetVerbose();
  // Name is not changed
  if(name == hadronPhysics) return;

  // Define hadronic process class
  G4String stand = G4String("Standard");
  G4String hado  = G4String("LowEnergy");
  G4String frag  = G4String("IonFragmentation");
  G4String none  = G4String("none");

  if(name == stand) {
    theHadList = new hTestHadronPhysicsList();
    theHadList->SetVerbose(verbose);

  } else if (name == hado) {
    theHadList = new hTestHadronPhysicsList1();
    theHadList->SetVerbose(verbose);

    //not implemented yet
  } else if (name == frag) {
    theHadList = 0;

  } else if (name == none) {
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






