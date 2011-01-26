
#include "EMElectronPenelope.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4eMultipleScattering.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4eIonisation.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4PenelopeBremsstrahlungModel.hh"
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
  
  G4eMultipleScattering* elecMultipScatProcess = new G4eMultipleScattering();

  G4eIonisation* elecIonisationProcess = new G4eIonisation();
  G4PenelopeIonisationModel* theIoniPenelope = 
    new G4PenelopeIonisationModel();
  elecIonisationProcess->AddEmModel(0,theIoniPenelope,new G4UniversalFluctuation());

  G4eBremsstrahlung* elecBremsstrProcess = new G4eBremsstrahlung();
  G4PenelopeBremsstrahlungModel* theBremPenelope = new 
    G4PenelopeBremsstrahlungModel();	
  elecBremsstrProcess->AddEmModel(0,theBremPenelope);

  G4StepLimiter* elecStepLimiter = new G4StepLimiter();

  elecMultipScatProcess -> SetRangeFactor(facRange);

  G4ParticleDefinition* particle = G4Electron::Electron(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddProcess(elecMultipScatProcess, -1, 1, 1);
  processManager -> AddProcess(elecIonisationProcess, -1, 2, 2);
  processManager -> AddProcess(elecBremsstrProcess, -1, -1, 3);
  processManager -> AddProcess(elecStepLimiter, -1, -1, 4);
}
