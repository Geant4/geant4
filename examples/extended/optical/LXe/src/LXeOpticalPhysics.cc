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
      pManager->AddContinuousProcess(theCerenkovProcess);
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

