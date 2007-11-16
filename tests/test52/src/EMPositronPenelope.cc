
#include "EMPositronPenelope.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4PenelopeIonisation.hh"
#include "G4PenelopeBremsstrahlung.hh"
#include "G4PenelopeAnnihilation.hh"
#include "G4StepLimiter.hh"


EMPositronPenelope::EMPositronPenelope(const G4String& name): 
  G4VPhysicsConstructor(name) { 

  G4cout << "ELECTROMAGNETIC PROCESS(ES): G4MultipleScattering (positron)" 
         << G4endl
         << "                             G4PenelopeIonisation (positron)" 
         << G4endl
         << "                             G4PenelopeBremsstrahlung (positron)" 
         << G4endl
         << "                             G4PenelopeAnnihilation (positron)" 
         << G4endl
         << "APPLIED MODEL(S): -" 
         << G4endl;

  facRange = 0.02;
}


EMPositronPenelope::~EMPositronPenelope() { 

}


void EMPositronPenelope::ConstructProcess() {

  // ****************
  // *** Positron ***
  // ****************

  G4MultipleScattering* positrMultipScatProcess = new G4MultipleScattering();
  G4PenelopeIonisation* positrIonisationProcess = new G4PenelopeIonisation();
  G4PenelopeBremsstrahlung* positrBremsstrProcess = 
                                               new G4PenelopeBremsstrahlung();
  G4PenelopeAnnihilation* positrAnnihilationProcess = 
                                               new G4PenelopeAnnihilation();

  G4StepLimiter* positrStepLimiter = new G4StepLimiter();

  positrMultipScatProcess -> SetRangeFactor(facRange);

  G4ParticleDefinition* particle = G4Positron::Positron(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddProcess(positrMultipScatProcess, -1, 1, 1);
  processManager -> AddProcess(positrIonisationProcess, -1, 2, 2);
  processManager -> AddProcess(positrBremsstrProcess, -1, -1, 3);
  processManager -> AddProcess(positrAnnihilationProcess, 0, -1, 4);
  processManager -> AddProcess(positrStepLimiter, -1, -1,  5);
}
