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
#include "IonPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// Nuclei
#include "G4IonConstructor.hh"

#include "globals.hh"
#include "G4ios.hh"

IonPhysics::IonPhysics(const G4String& name)
                 :  G4VPhysicsConstructor(name), wasActivated(false)
{
}

IonPhysics::~IonPhysics()
{
  if(wasActivated)
  {
  G4ProcessManager * pManager = 0;
  
  pManager = G4GenericIon::GenericIon()->GetProcessManager();
  //if(pManager) pManager->RemoveProcess(&theIonElasticProcess);
  //if(pManager) pManager->RemoveProcess(&fIonIonisation);
  // if(pManager) pManager->RemoveProcess(&fIonMultipleScattering);

  pManager = G4Deuteron::Deuteron()->GetProcessManager();
  if(pManager) pManager->RemoveProcess(&theDElasticProcess);
  if(pManager) pManager->RemoveProcess(&fDeuteronProcess);
  if(pManager) pManager->RemoveProcess(&fDeuteronIonisation);
  if(pManager) pManager->RemoveProcess(&fDeuteronMultipleScattering);
 
  pManager = G4Triton::Triton()->GetProcessManager();
  if(pManager) pManager->RemoveProcess(&theTElasticProcess);
  if(pManager) pManager->RemoveProcess(&fTritonProcess);
  if(pManager) pManager->RemoveProcess(&fTritonIonisation);
  if(pManager) pManager->RemoveProcess(&fTritonMultipleScattering);
 
  pManager = G4Alpha::Alpha()->GetProcessManager();
  if(pManager) pManager->RemoveProcess(&theAElasticProcess);
  if(pManager) pManager->RemoveProcess(&fAlphaProcess);
  if(pManager) pManager->RemoveProcess(&fAlphaIonisation);
  if(pManager) pManager->RemoveProcess(&fAlphaMultipleScattering);
 
  pManager = G4He3::He3()->GetProcessManager();
  //if(pManager) pManager->RemoveProcess(&theHe3ElasticProcess);
  //if(pManager) pManager->RemoveProcess(&fHe3Ionisation);
  //if(pManager) pManager->RemoveProcess(&fHe3MultipleScattering);
  }
}

void IonPhysics::ConstructParticle()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void IonPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;
  
  // Elastic Process
  theElasticModel = new G4LElastic();
  theIonElasticProcess.RegisterMe(theElasticModel);
  theDElasticProcess.RegisterMe(theElasticModel);
  theTElasticProcess.RegisterMe(theElasticModel);
  theAElasticProcess.RegisterMe(theElasticModel);
  theHe3ElasticProcess.RegisterMe(theElasticModel);

  // Generic Ion
  pManager = G4GenericIon::GenericIon()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theIonElasticProcess);

  pManager->AddProcess(&fIonIonisation, ordInActive, 2, 2);

  // pManager->AddProcess(&fIonMultipleScattering);
  // pManager->SetProcessOrdering(&fIonMultipleScattering, idxAlongStep,  1);
  // pManager->SetProcessOrdering(&fIonMultipleScattering, idxPostStep,  1);

  // Deuteron 
  pManager = G4Deuteron::Deuteron()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theDElasticProcess);

  fDeuteronModel = new G4LEDeuteronInelastic();
  fDeuteronProcess.RegisterMe(fDeuteronModel);
  pManager->AddDiscreteProcess(&fDeuteronProcess);

  pManager->AddProcess(&fDeuteronIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fDeuteronMultipleScattering);
  pManager->SetProcessOrdering(&fDeuteronMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fDeuteronMultipleScattering, idxPostStep,  1);
 
  // Triton 
  pManager = G4Triton::Triton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theTElasticProcess);

  fTritonModel = new G4LETritonInelastic();
  fTritonProcess.RegisterMe(fTritonModel);
  pManager->AddDiscreteProcess(&fTritonProcess);

  pManager->AddProcess(&fTritonIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fTritonMultipleScattering);
  pManager->SetProcessOrdering(&fTritonMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fTritonMultipleScattering, idxPostStep,  1);
 
  // Alpha 
  pManager = G4Alpha::Alpha()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theAElasticProcess);

  fAlphaModel = new G4LEAlphaInelastic();
  fAlphaProcess.RegisterMe(fAlphaModel);
  pManager->AddDiscreteProcess(&fAlphaProcess);

  pManager->AddProcess(&fAlphaIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fAlphaMultipleScattering);
  pManager->SetProcessOrdering(&fAlphaMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fAlphaMultipleScattering, idxPostStep,  1);
 
  // He3
  pManager = G4He3::He3()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theHe3ElasticProcess);

  pManager->AddProcess(&fHe3Ionisation, ordInActive, 2, 2);

  pManager->AddProcess(&fHe3MultipleScattering);
  pManager->SetProcessOrdering(&fHe3MultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fHe3MultipleScattering, idxPostStep,  1);
   
  wasActivated = true;
}



// 2002 by J.P. Wellisch
