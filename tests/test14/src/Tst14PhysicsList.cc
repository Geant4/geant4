// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14PhysicsList.cc,v 1.1 1999-05-29 14:12:11 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
   
#include "G4Timer.hh"

#include "Tst14PhysicsList.hh"
#include "Tst14DetectorConstruction.hh"
#include "Tst14PhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4EnergyLossTables.hh"
#include "G4ios.hh"
#include <iomanip.h>                

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst14PhysicsList::Tst14PhysicsList(Tst14DetectorConstruction* p)
  :  G4VUserPhysicsList(), theLowEnergyBremstrahlung(0),
     theLowEnergyPhotoElectric(0), theLowEnergyCompton(0),
     theLowEnergyGammaConversion(0), theLowEnergyRayleigh(0),
     //
     theeminusIonisation(NULL), theeminusBremsstrahlung(NULL),
     theeplusIonisation(NULL), theeplusBremsstrahlung(NULL),
     theeminusStepCut(NULL), theeplusStepCut(NULL),
     MaxChargedStep(DBL_MAX)
{
  pDet = p;

  defaultCutValue = 1.000*mm;
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  cutForProton = defaultCutValue;

  SetVerboseLevel(1);
  physicsListMessenger = new Tst14PhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst14PhysicsList::~Tst14PhysicsList()
{
  delete physicsListMessenger; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBarions();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PhysicsList::ConstructBosons()
{
  // gamma
  G4Gamma::GammaDefinition();
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PhysicsList::ConstructLeptons()
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

void Tst14PhysicsList::ConstructMesons()
{
 //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PhysicsList::ConstructBarions()
{
//  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyRayleigh.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyBremsstrahlung.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"

#include "Tst14StepCut.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma

      // Construct processes for gamma
      theLowEnergyPhotoElectric = new G4LowEnergyPhotoElectric();      
      theLowEnergyCompton = new G4LowEnergyCompton();
      theLowEnergyRayleigh = new G4LowEnergyRayleigh();
      theLowEnergyGammaConversion = new G4LowEnergyGammaConversion();
      
      pmanager->AddDiscreteProcess(theLowEnergyPhotoElectric);
      pmanager->AddDiscreteProcess(theLowEnergyCompton);
      pmanager->AddDiscreteProcess(theLowEnergyRayleigh);
      pmanager->AddDiscreteProcess(theLowEnergyGammaConversion);
      
    } else if (particleName == "e-") {
      // Construct processes for electron <---------------------------
      theeminusIonisation = new G4eIonisation();
      theeminusBremsstrahlung = new G4eBremsstrahlung();
      theLowEnergyBremstrahlung = new G4LowEnergyBremsstrahlung();
      theeminusStepCut = new Tst14StepCut();

      pmanager->AddProcess(theLowEnergyBremstrahlung);
      pmanager->SetProcessOrdering(theLowEnergyBremstrahlung,idxAlongStep,1);
      pmanager->SetProcessOrdering(theLowEnergyBremstrahlung,idxPostStep,1);

      pmanager->AddProcess(theeminusIonisation,-1,2,2);
      pmanager->AddProcess(theeminusBremsstrahlung,-1,-1,3);      
      pmanager->AddProcess(theeminusStepCut,-1,-1,4);
        theeminusStepCut->SetMaxStep(MaxChargedStep) ;

    } else if (particleName == "e+") {
    // Construct processes for positron <------------------------------
      theeplusIonisation = new G4eIonisation();
      theeplusBremsstrahlung = new G4eBremsstrahlung();
      theeplusStepCut = new Tst14StepCut();
      
      pmanager->AddProcess(theeplusIonisation,-1,2,2);
      pmanager->AddProcess(theeplusBremsstrahlung,-1,-1,3);
      pmanager->AddProcess(theeplusStepCut,-1,-1,5);
      theeplusStepCut->SetMaxStep(MaxChargedStep);
  
    } 
  }
}


#include "G4Decay.hh"

void Tst14PhysicsList::ConstructGeneral(){

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

void Tst14PhysicsList::SetCuts()
{
  G4Timer theTimer ;
  theTimer.Start() ;
  if (verboseLevel >0){
    G4cout << "Tst14PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << endl;
  }  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 

  SetCutValue(cutForGamma,"gamma");
  
  SetCutValue(cutForElectron,"e-");
  SetCutValue(cutForElectron,"e+");
  
  SetCutValue(defaultCutValue,"mu-");
  SetCutValue(defaultCutValue,"mu+");
  
  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton 
  
  SetCutValue(defaultCutValue, "proton");
  SetCutValue(defaultCutValue, "anti_proton");
  
  SetCutValueForOthers(defaultCutValue);
  
  if (verboseLevel>1) {
    DumpCutValuesTable();
  }

  theTimer.Stop();
  G4cout.precision(6);
  G4cout << endl ;
  G4cout << "total time(SetCuts)=" << theTimer.GetUserElapsed() << " s " <<endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  //  G4Gamma::SetEnergyRange(2.5e-4*MeV,1e5*MeV);
  cutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PhysicsList::SetElectronCut(G4double val)
{
  ResetCuts();
  // G4Electron::SetEnergyRange(2.5e-4*MeV,1e5*MeV);
  cutForElectron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PhysicsList::SetProtonCut(G4double val)
{
  ResetCuts();
  cutForProton = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PhysicsList::SetCutsByEnergy(G4double val)
{

  G4Material* currMat = pDet->GetAbsorberMaterial();

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  G4ParticleDefinition* part;
  G4double cut;

  part = theParticleTable->FindParticle("e-");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "e-");

  part = theParticleTable->FindParticle("e+");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "e+");

  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton
  part = theParticleTable->FindParticle("proton");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "proton");

  part = theParticleTable->FindParticle("anti_proton");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "anti_proton");

  SetCutValueForOthers(cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PhysicsList::GetRange(G4double val)
{
  G4Material* currMat = pDet->GetAbsorberMaterial();

  G4ParticleDefinition* part;
  G4double cut;
  part = theParticleTable->FindParticle("e-");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  G4cout << "material : " << currMat->GetName() << endl;
  G4cout << "particle : " << part->GetParticleName() << endl;
  G4cout << "energy   : " << val / keV << " (keV)" << endl;
  G4cout << "range    : " << cut / mm << " (mm)" << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst14PhysicsList::SetMaxStep(G4double step)
{
  MaxChargedStep = step ;
  G4cout << " MaxChargedStep=" << MaxChargedStep << endl;
  G4cout << endl;
}



















