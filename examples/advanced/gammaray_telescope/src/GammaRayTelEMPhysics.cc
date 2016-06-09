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
// $Id: GammaRayTelEMPhysics.cc,v 1.4 2006/06/29 15:56:36 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// 

#include "GammaRayTelEMPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   


GammaRayTelEMPhysics::GammaRayTelEMPhysics(const G4String& name)
               :  G4VPhysicsConstructor(name)
{
}

GammaRayTelEMPhysics::~GammaRayTelEMPhysics()
{
}

void GammaRayTelEMPhysics::ConstructParticle()
{
}


#include "G4ProcessManager.hh"


void GammaRayTelEMPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;
  
  // Gamma Physics
  pManager = G4Gamma::Gamma()->GetProcessManager();

  // std

  pManager->AddDiscreteProcess(&thePhotoEffect);
  pManager->AddDiscreteProcess(&theComptonEffect);
  pManager->AddDiscreteProcess(&thePairProduction);

  // lowe

  pManager->AddDiscreteProcess(&theLowEnPhoto);
  pManager->AddDiscreteProcess(&theLowEnCompton);
  pManager->AddDiscreteProcess(&theLowEnPair);
  pManager->AddDiscreteProcess(&theLowEnRayleigh);
  
  // Electron Physics
  pManager = G4Electron::Electron()->GetProcessManager();

   // add processes

  pManager->AddDiscreteProcess(&theElectronBremsStrahlung);  
  pManager->AddProcess(&theElectronIonisation, ordInActive,2, 2);

  pManager->AddDiscreteProcess(&theLowEnBremss);
  pManager->AddProcess(&theLowEnIon, ordInActive,2, 2);


  pManager->AddProcess(&theElectronMultipleScattering);
  pManager->SetProcessOrdering(&theElectronMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&theElectronMultipleScattering, idxPostStep,  1);



  //Positron Physics
  pManager = G4Positron::Positron()->GetProcessManager();
  // add processes
  pManager->AddDiscreteProcess(&thePositronBremsStrahlung);

  pManager->AddDiscreteProcess(&theAnnihilation);

  pManager->AddRestProcess(&theAnnihilation);

  pManager->AddProcess(&thePositronIonisation, ordInActive,2, 2);

  pManager->AddProcess(&thePositronMultipleScattering);
  pManager->SetProcessOrdering(&thePositronMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&thePositronMultipleScattering, idxPostStep,  1);

}



