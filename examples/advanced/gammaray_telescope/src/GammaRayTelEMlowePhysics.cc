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
//
//
// 

#include "GammaRayTelEMlowePhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   

// gamma

#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"


// e-/e+

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4eplusAnnihilation.hh"


GammaRayTelEMlowePhysics::GammaRayTelEMlowePhysics(const G4String& name)
               :  G4VPhysicsConstructor(name)
{
}

GammaRayTelEMlowePhysics::~GammaRayTelEMlowePhysics()
{
}

void GammaRayTelEMlowePhysics::ConstructParticle()
{
}


#include "G4ProcessManager.hh"


void GammaRayTelEMlowePhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;
  
  // Gamma Physics
  pManager = G4Gamma::Gamma()->GetProcessManager();

  G4RayleighScattering* theRayleigh = new G4RayleighScattering();
  theRayleigh->SetEmModel(new G4LivermoreRayleighModel()); 
  G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
  thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
  G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
  theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
  G4GammaConversion* theGammaConversion = new G4GammaConversion();
  theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel());

  pManager->AddDiscreteProcess(theRayleigh);
  pManager->AddDiscreteProcess(thePhotoElectricEffect);
  pManager->AddDiscreteProcess(theComptonScattering);
  pManager->AddDiscreteProcess(theGammaConversion);

  
  // Electron Physics

  pManager = G4Electron::Electron()->GetProcessManager();


  G4eMultipleScattering* msc = new G4eMultipleScattering();

  G4eIonisation* eIonisation = new G4eIonisation();
  eIonisation->SetEmModel(new G4LivermoreIonisationModel());

  G4eBremsstrahlung* eBremsstrahlung = new G4eBremsstrahlung();
  eBremsstrahlung->SetEmModel(new G4LivermoreBremsstrahlungModel());

  pManager->AddProcess(msc,-1, 1, 1);
  pManager->AddProcess(eIonisation,-1, 2, 2);
  pManager->AddProcess(eBremsstrahlung, -1,-3, 3);


  //Positron Physics

  pManager = G4Positron::Positron()->GetProcessManager();

  G4eMultipleScattering* msc_p = new G4eMultipleScattering();
  pManager->AddProcess(msc_p,-1, 1, 1);
  
  G4eIonisation* eIonisation_p = new G4eIonisation();
  pManager->AddProcess(eIonisation_p, -1, 2, 2);
  pManager->AddProcess(new G4eBremsstrahlung(), -1,-1, 3);

  pManager->AddProcess(new G4eplusAnnihilation(),0,-1, 4);

}



