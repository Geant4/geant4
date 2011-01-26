
#include "EMElectronEEDL.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4eMultipleScattering.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4StepLimiter.hh"


EMElectronEEDL::EMElectronEEDL(const G4String& name): 
    G4VPhysicsConstructor(name) {
 
  G4cout << "ELECTROMAGNETIC PROCESS(ES): G4eMultipleScattering (electron)" 
         << G4endl
         << "                             G4LivermoreIonisationModel (electron)" 
         << G4endl
         << "                             G4LivermoreBremsstrahlungModel (electron)" 
         << G4endl
         << "APPLIED MODEL(S): -" 
         << G4endl;

  facRange = 0.02;
}


EMElectronEEDL::~EMElectronEEDL() { 

}


void EMElectronEEDL::ConstructProcess() {

  // ****************
  // *** Electron ***
  // ****************

  G4eMultipleScattering* elecMultipScatProcess = new G4eMultipleScattering();
  elecMultipScatProcess->AddEmModel(0, new G4GoudsmitSaundersonMscModel());

  G4eIonisation* elecIonisationProcess = new G4eIonisation();
  G4LivermoreIonisationModel* theIoniLivermore = new
    G4LivermoreIonisationModel();
  elecIonisationProcess->AddEmModel(0, theIoniLivermore, new G4UniversalFluctuation() );
  elecIonisationProcess->SetStepFunction(0.2, 100*um); //     

  G4eBremsstrahlung* elecBremsstrProcess = new G4eBremsstrahlung();
  G4LivermoreBremsstrahlungModel* theBremLivermore = new
    G4LivermoreBremsstrahlungModel();
  elecBremsstrProcess->AddEmModel(0, theBremLivermore);

  G4StepLimiter* elecStepLimiter = new G4StepLimiter();

  elecMultipScatProcess -> SetRangeFactor(facRange);

  G4ParticleDefinition* particle = G4Electron::Electron(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddProcess(elecMultipScatProcess, -1, 1, 1);
  processManager -> AddProcess(elecIonisationProcess, -1, 2, 2);
  processManager -> AddProcess(elecBremsstrProcess, -1, -1, 3);
  processManager -> AddProcess(elecStepLimiter, -1, -1, 4);
}
