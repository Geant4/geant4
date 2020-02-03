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

#include "GammaRayTelEMstdPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   

// gamma

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"


// e-
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"

// e+
#include "G4eIonisation.hh" 
#include "G4eBremsstrahlung.hh" 
#include "G4eplusAnnihilation.hh"


GammaRayTelEMstdPhysics::GammaRayTelEMstdPhysics(const G4String& name)
               :  G4VPhysicsConstructor(name)
{
}

GammaRayTelEMstdPhysics::~GammaRayTelEMstdPhysics()
{
}

void GammaRayTelEMstdPhysics::ConstructParticle()
{
}


#include "G4ProcessManager.hh"


void GammaRayTelEMstdPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;
  
  // Gamma Physics
  pManager = G4Gamma::Gamma()->GetProcessManager();

  pManager->AddDiscreteProcess(new G4PhotoElectricEffect);
  pManager->AddDiscreteProcess(new G4ComptonScattering);
  pManager->AddDiscreteProcess(new G4GammaConversion);

     
  // Electron Physics

  pManager = G4Electron::Electron()->GetProcessManager();

  pManager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
  pManager->AddProcess(new G4eIonisation,         -1, 2, 2);
  pManager->AddProcess(new G4eBremsstrahlung,     -1, 3, 3);
	    
  // Positron Physics

  pManager = G4Positron::Positron()->GetProcessManager();

  pManager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
  pManager->AddProcess(new G4eIonisation,         -1, 2, 2);
  pManager->AddProcess(new G4eBremsstrahlung,     -1, 3, 3);
  pManager->AddProcess(new G4eplusAnnihilation,    0,-1, 4);

}


