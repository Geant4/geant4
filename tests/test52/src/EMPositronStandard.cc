
#include "EMPositronStandard.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4StepLimiter.hh"
#include "G4EmProcessOptions.hh"


EMPositronStandard::EMPositronStandard(const G4String& name): 
   G4VPhysicsConstructor(name) {
 
  G4cout << "ELECTROMAGNETIC PROCESS(ES): G4MultipleScattering (positron)" 
         << G4endl
         << "                             G4eIonisation (positron)" 
         << G4endl
         << "                             G4eBremsstrahlung (positron)" 
         << G4endl
         << "                             G4eplusAnnihilation (positron)" 
         << G4endl
         << "APPLIED MODEL(S): -" 
         << G4endl;
 
  facRange = 0.02;
}


EMPositronStandard::~EMPositronStandard() { 

}


void EMPositronStandard::ConstructProcess() {

  // ****************
  // *** Positron ***
  // ****************

  // G4EmProcessOptions* positronEmProcessOptions = new G4EmProcessOptions();
  // positronEmProcessOptions -> SetDEDXBinning(480);

  G4eMultipleScattering* positrMultipScatProcess = new G4eMultipleScattering();
  G4eIonisation* positrIonisationProcess = new G4eIonisation();
  G4eBremsstrahlung* positrBremsstrProcess = new G4eBremsstrahlung();
  G4eplusAnnihilation* positrAnnihilationProcess = new G4eplusAnnihilation();

  G4StepLimiter* positrStepLimiter = new G4StepLimiter();

  positrMultipScatProcess -> SetRangeFactor(facRange);

  G4ParticleDefinition* particle = G4Positron::Positron(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddProcess(positrMultipScatProcess, -1, 1, 1);
  processManager -> AddProcess(positrIonisationProcess, -1, 2, 2);
  processManager -> AddProcess(positrBremsstrProcess, -1, -1, 3);
  processManager -> AddProcess(positrAnnihilationProcess, 0, -1, 4);
  processManager -> AddProcess(positrStepLimiter, -1, -1, 5);
}
