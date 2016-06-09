//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "LXeEMPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   


LXeEMPhysics::LXeEMPhysics(const G4String& name)
               :  G4VPhysicsConstructor(name)
{
}

LXeEMPhysics::~LXeEMPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

void LXeEMPhysics::ConstructParticle()
{
  // gamma
  G4Gamma::GammaDefinition();
 
  // electron
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
}


#include "G4ProcessManager.hh"


void LXeEMPhysics::ConstructProcess()
{
  thePhotoEffect = new G4PhotoElectricEffect();
  theComptonEffect = new G4ComptonScattering();
  thePairProduction = new G4GammaConversion();
  
    // Electron physics
  theElectronMultipleScattering = new G4MultipleScattering();
  theElectronIonisation = new G4eIonisation();
  theElectronBremsStrahlung = new G4eBremsstrahlung();
  
    //Positron physics
  thePositronMultipleScattering = new G4MultipleScattering();
  thePositronIonisation = new G4eIonisation(); 
  thePositronBremsStrahlung = new G4eBremsstrahlung();  
  theAnnihilation = new G4eplusAnnihilation();


  G4ProcessManager * pManager = 0;
  
  // Gamma Physics
  pManager = G4Gamma::Gamma()->GetProcessManager();
  pManager->AddDiscreteProcess(thePhotoEffect);
  pManager->AddDiscreteProcess(theComptonEffect);
  pManager->AddDiscreteProcess(thePairProduction);

  // Electron Physics
  pManager = G4Electron::Electron()->GetProcessManager();

  pManager->AddProcess(theElectronMultipleScattering, -1, 1, 1);
  pManager->AddProcess(theElectronIonisation,         -1, 2, 2);
  pManager->AddProcess(theElectronBremsStrahlung,     -1, 3, 3);  


  //Positron Physics
  pManager = G4Positron::Positron()->GetProcessManager();
 
  pManager->AddProcess(thePositronMultipleScattering, -1, 1, 1);
  pManager->AddProcess(thePositronIonisation,         -1, 2, 2);
  pManager->AddProcess(thePositronBremsStrahlung,     -1, 3, 3);  
  pManager->AddProcess(theAnnihilation,                0,-1, 4);  

}



