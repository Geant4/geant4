//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "EMPhotonEPDL.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4StepLimiter.hh"


EMPhotonEPDL::EMPhotonEPDL(const G4String& name): 
   G4VPhysicsConstructor(name) {
 
  G4cout << "ELECTROMAGNETIC PROCESS(ES): G4LivermorePhotoElectric (photon)" 
         << G4endl
         << "                             G4LivermoreCompton (photon)" 
         << G4endl
         << "                             G4LivermoreGammaConversion (photon)" 
         << G4endl
         << "                             G4LivermoreRayleigh (photon)" 
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
  G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
  G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel = 
    new G4LivermorePhotoElectricModel();
  thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);

  G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
  G4LivermoreComptonModel* theLivermoreComptonModel = 
    new G4LivermoreComptonModel();
  theComptonScattering->AddEmModel(0, theLivermoreComptonModel);


  G4GammaConversion* theGammaConversion = new G4GammaConversion();
  G4LivermoreGammaConversionModel* theLivermoreGammaConversionModel = 
    new G4LivermoreGammaConversionModel();  
  theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);

  G4RayleighScattering* theRayleigh = new G4RayleighScattering();
  G4LivermoreRayleighModel* theRayleighModel = new G4LivermoreRayleighModel();
  theRayleigh->AddEmModel(0, theRayleighModel);

  G4StepLimiter* photStepLimiter = new G4StepLimiter();

  G4ParticleDefinition* particle = G4Gamma::Gamma(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(thePhotoElectricEffect);
  processManager -> AddDiscreteProcess(theComptonScattering);
  processManager -> AddDiscreteProcess(theGammaConversion);
  processManager -> AddDiscreteProcess(theRayleigh);
  processManager -> AddDiscreteProcess(photStepLimiter);

}
