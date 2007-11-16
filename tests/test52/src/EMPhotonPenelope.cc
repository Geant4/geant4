
#include "EMPhotonPenelope.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4PenelopeCompton.hh"
#include "G4PenelopeGammaConversion.hh"
#include "G4PenelopePhotoElectric.hh"
#include "G4PenelopeRayleigh.hh"
#include "G4StepLimiter.hh"


EMPhotonPenelope::EMPhotonPenelope(const G4String& name): 
   G4VPhysicsConstructor(name) { 

  G4cout << "ELECTROMAGNETIC PROCESS(ES): G4PenelopePhotoElectric (photon)" 
         << G4endl
         << "                             G4PenelopeCompton (photon)" 
         << G4endl
         << "                             G4PenelopeGammaConversion (photon)" 
         << G4endl
         << "                             G4PenelopeRayleigh (photon)" 
         << G4endl
         << "APPLIED MODEL(S): -" 
         << G4endl;
}


EMPhotonPenelope::~EMPhotonPenelope() { 

}


void EMPhotonPenelope::ConstructProcess() {

  // **************
  // *** Photon ***
  // **************
  
  G4PenelopePhotoElectric* photPhotoElectricProcess = 
                                          new G4PenelopePhotoElectric();
  G4PenelopeCompton* photComptonProcess = new G4PenelopeCompton;
  G4PenelopeGammaConversion* photGammaConvProcess = 
                                          new G4PenelopeGammaConversion;
  G4PenelopeRayleigh* photRayleighProcess = new G4PenelopeRayleigh;

  G4StepLimiter* photStepLimiter = new G4StepLimiter();

  G4ParticleDefinition* particle = G4Gamma::Gamma(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(photPhotoElectricProcess);
  processManager -> AddDiscreteProcess(photComptonProcess);
  processManager -> AddDiscreteProcess(photGammaConvProcess);
  processManager -> AddDiscreteProcess(photRayleighProcess);
  processManager -> AddDiscreteProcess(photStepLimiter);
}


