
#include "EMPhotonPenelope.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4PenelopePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4PenelopeComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4RayleighScattering.hh"
#include "G4PenelopeRayleighModel.hh"
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
  G4ParticleDefinition* particle = G4Gamma::Gamma(); 
  G4ProcessManager* manager = particle -> GetProcessManager();

  //Photo-electric effect
  G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
  G4PenelopePhotoElectricModel* thePEPenelopeModel = new 
    G4PenelopePhotoElectricModel();
  thePhotoElectricEffect->AddEmModel(0,thePEPenelopeModel);
  manager->AddDiscreteProcess(thePhotoElectricEffect);

  //Compton scattering
  G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
  G4PenelopeComptonModel* theComptonPenelopeModel = 
    new G4PenelopeComptonModel();
  theComptonScattering->AddEmModel(0,theComptonPenelopeModel);
  manager->AddDiscreteProcess(theComptonScattering);
	  
  //Gamma conversion
  G4GammaConversion* theGammaConversion = new G4GammaConversion();
  G4PenelopeGammaConversionModel* theGCPenelopeModel = 
    new G4PenelopeGammaConversionModel();
  theGammaConversion->AddEmModel(0,theGCPenelopeModel);
  manager->AddDiscreteProcess(theGammaConversion);
	  
  //Rayleigh scattering
  G4RayleighScattering* theRayleigh = new G4RayleighScattering();
  G4PenelopeRayleighModel* theRayleighPenelopeModel = 
    new G4PenelopeRayleighModel();  
  theRayleigh->AddEmModel(0,theRayleighPenelopeModel);
  manager->AddDiscreteProcess(theRayleigh);
}


