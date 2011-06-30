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
