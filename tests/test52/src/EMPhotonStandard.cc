
#include "EMPhotonStandard.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4StepLimiter.hh"
// #include "G4EmProcessOptions.hh"


EMPhotonStandard::EMPhotonStandard(const G4String& name): 
   G4VPhysicsConstructor(name) {
 
  G4cout << "ELECTROMAGNETIC PROCESS(ES): G4PhotoElectricEffect (photon)" 
         << G4endl
         << "                             G4ComptonScattering (photon)" 
         << G4endl
         << "                             G4GammaConversion (photon)" 
         << G4endl
         << "APPLIED MODEL(S): -" 
         << G4endl;
}


EMPhotonStandard::~EMPhotonStandard() { 

}


void EMPhotonStandard::ConstructProcess() {

  // **************
  // *** Photon ***
  // **************

  // G4EmProcessOptions* photonEmProcessOptions = new G4EmProcessOptions();
  // photonEmProcessOptions -> SetDEDXBinning(480);
   
  G4PhotoElectricEffect* photPhotoElectricProcess = 
                                            new G4PhotoElectricEffect();
  G4ComptonScattering* photComptonProcess = new G4ComptonScattering;
  G4GammaConversion* photGammaConvProcess = new G4GammaConversion;

  G4StepLimiter* photStepLimiter = new G4StepLimiter();

  G4ParticleDefinition* particle = G4Gamma::Gamma(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(photPhotoElectricProcess);
  processManager -> AddDiscreteProcess(photComptonProcess);
  processManager -> AddDiscreteProcess(photGammaConvProcess);
  processManager -> AddDiscreteProcess(photStepLimiter);
}
