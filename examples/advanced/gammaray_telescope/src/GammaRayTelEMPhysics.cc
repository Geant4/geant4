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
//
// $Id: GammaRayTelEMPhysics.cc,v 1.2 2003/06/16 16:46:24 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

void GammaRayTelEMPhysics::ConstructParticle()
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



