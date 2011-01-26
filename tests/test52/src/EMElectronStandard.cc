
#include "EMElectronStandard.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4StepLimiter.hh"
//#include "G4EmProcessOptions.hh"


EMElectronStandard::EMElectronStandard(const G4String& name): 
   G4VPhysicsConstructor(name) {
 
  G4cout << "ELECTROMAGNETIC PROCESS(ES): G4MultipleScattering (electron)" 
         << G4endl
         << "                             G4eIonisation (electron)" 
         << G4endl
         << "                             G4eBremsstrahlung (electron)" 
         << G4endl
         << "APPLIED MODEL(S): -" 
         << G4endl;

  facRange = 0.02;
}


EMElectronStandard::~EMElectronStandard() { 

}


void EMElectronStandard::ConstructProcess() {
  
  // ****************
  // *** Electron ***
  // ****************

  // G4EmProcessOptions* elecEmProcessOptions = new G4EmProcessOptions();
  // elecEmProcessOptions -> SetDEDXBinning(480);

  G4eMultipleScattering* elecMultipScatProcess = new G4eMultipleScattering();
  G4eIonisation* elecIonisationProcess = new G4eIonisation();
  G4eBremsstrahlung* elecBremsstrProcess = new G4eBremsstrahlung();

  G4StepLimiter* elecStepLimiter = new G4StepLimiter();

  elecMultipScatProcess -> SetRangeFactor(facRange);

  G4ParticleDefinition* particle = G4Electron::Electron(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddProcess(elecMultipScatProcess, -1, 1, 1);
  processManager -> AddProcess(elecIonisationProcess, -1, 2, 2);
  processManager -> AddProcess(elecBremsstrProcess, -1, -1, 3);
  processManager -> AddProcess(elecStepLimiter, -1, -1, 4);
}
