
#include "EMPhotonEPDL.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"
#include "G4StepLimiter.hh"


EMPhotonEPDL::EMPhotonEPDL(const G4String& name): 
   G4VPhysicsConstructor(name) {
 
  G4cout << "ELECTROMAGNETIC PROCESS(ES): G4LowEnergyPhotoElectric (photon)" 
         << G4endl
         << "                             G4LowEnergyCompton (photon)" 
         << G4endl
         << "                             G4LowEnergyGammaConversion (photon)" 
         << G4endl
         << "                             G4LowEnergyRayleigh (photon)" 
         << G4endl
         << "APPLIED MODEL(S): -" 
         << G4endl;
}


EMPhotonEPDL::~EMPhotonEPDL() { 

}


void EMPhotonEPDL::ConstructProcess() {

  // **************
  // *** Photon ***
  // **************
  
  G4LowEnergyPhotoElectric* photPhotoElectricProcess = 
                                   new G4LowEnergyPhotoElectric();
  G4LowEnergyCompton* photComptonProcess = new G4LowEnergyCompton;
  G4LowEnergyGammaConversion* photGammaConvProcess = 
                                   new G4LowEnergyGammaConversion;
  G4LowEnergyRayleigh* photRayleighProcess = new G4LowEnergyRayleigh;

  G4StepLimiter* photStepLimiter = new G4StepLimiter();

  G4ParticleDefinition* particle = G4Gamma::Gamma(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(photPhotoElectricProcess);
  processManager -> AddDiscreteProcess(photComptonProcess);
  processManager -> AddDiscreteProcess(photGammaConvProcess);
  processManager -> AddDiscreteProcess(photRayleighProcess);
  processManager -> AddDiscreteProcess(photStepLimiter);

}
