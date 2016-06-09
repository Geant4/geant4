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
#include "LXeOpticalPhysics.hh"
#include "LXeOpPhysMessenger.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"

LXeOpticalPhysics::LXeOpticalPhysics(const G4String& name)
  :  G4VPhysicsConstructor(name)
{
  //  new LXeOpPhysMessenger(this);
}

LXeOpticalPhysics::~LXeOpticalPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4OpticalPhoton.hh"

void LXeOpticalPhysics::ConstructParticle()
{
  G4OpticalPhoton::OpticalPhotonDefinition();
}


#include "G4ProcessManager.hh"

void LXeOpticalPhysics::ConstructProcess()
{
  theScintProcess = new G4Scintillation();
  theCerenkovProcess=new G4Cerenkov();
  theAbsorptionProcess=new G4OpAbsorption();
  theRayleighScattering=new G4OpRayleigh();
  theBoundaryProcess=new G4OpBoundaryProcess();
  theWLSProcess=new G4OpWLS();

  theWLSProcess->UseTimeProfile("delta");
//theWLSProcess->UseTimeProfile("exponential");

  G4ProcessManager * pManager = 0;
  
  pManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
  pManager->AddDiscreteProcess(theAbsorptionProcess);
  pManager->AddDiscreteProcess(theRayleighScattering);
  theBoundaryProcess->SetModel(unified);
  pManager->AddDiscreteProcess(theBoundaryProcess);
  pManager->AddDiscreteProcess(theWLSProcess);
  
  theScintProcess->SetScintillationYieldFactor(1.);
  theScintProcess->SetScintillationExcitationRatio(0.0);
  theScintProcess->SetTrackSecondariesFirst(true);

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    pManager = particle->GetProcessManager();
    if(theCerenkovProcess->IsApplicable(*particle)){
      pManager->AddProcess(theCerenkovProcess);
      pManager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
    }
    if(theScintProcess->IsApplicable(*particle)){
      pManager->AddProcess(theScintProcess);
      pManager->SetProcessOrderingToLast(theScintProcess,idxAtRest);
      pManager->SetProcessOrderingToLast(theScintProcess,idxPostStep);
    }
  
  }
}

void LXeOpticalPhysics::SetScintYieldFactor(G4double yf){
  if(theScintProcess)
    theScintProcess->SetScintillationYieldFactor(yf);
}


