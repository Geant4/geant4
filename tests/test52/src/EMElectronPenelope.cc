
#include "EMElectronPenelope.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4PenelopeIonisation.hh"
#include "G4PenelopeBremsstrahlung.hh"
#include "G4StepLimiter.hh"


EMElectronPenelope::EMElectronPenelope(const G4String& name): 
   G4VPhysicsConstructor(name) {
 
  G4cout << "ELECTROMAGNETIC PROCESS(ES): G4MultipleScattering (electron)" 
         << G4endl
         << "                             G4PenelopeIonisation (electron)" 
         << G4endl
         << "                             G4PenelopeBremsstrahlung (electron)" 
         << G4endl
         << "APPLIED MODEL(S): -" 
         << G4endl;

  facRange = 0.02;  
}


EMElectronPenelope::~EMElectronPenelope() { 

}


void EMElectronPenelope::ConstructProcess() {

  // ****************
  // *** Electron ***
  // ****************
  
  G4MultipleScattering* elecMultipScatProcess = new G4MultipleScattering();
  G4PenelopeIonisation* elecIonisationProcess = new G4PenelopeIonisation();
  G4PenelopeBremsstrahlung* elecBremsstrProcess = new G4PenelopeBremsstrahlung();

  G4StepLimiter* elecStepLimiter = new G4StepLimiter();

  elecMultipScatProcess -> SetRangeFactor(facRange);

  G4ParticleDefinition* particle = G4Electron::Electron(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddProcess(elecMultipScatProcess, -1, 1, 1);
  processManager -> AddProcess(elecIonisationProcess, -1, 2, 2);
  processManager -> AddProcess(elecBremsstrProcess, -1, -1, 3);
  processManager -> AddProcess(elecStepLimiter, -1, -1, 4);
}
