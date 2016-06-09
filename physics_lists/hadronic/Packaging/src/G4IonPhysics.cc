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
// $Id: G4IonPhysics.cc,v 1.2 2005/12/02 12:40:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonPhysics
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4IonPhysics.hh"
#include "G4HadronElasticProcess.hh"
#include "G4LElastic.hh"

#include "G4DeuteronInelasticProcess.hh"
#include "G4LEDeuteronInelastic.hh"

#include "G4TritonInelasticProcess.hh"
#include "G4LETritonInelastic.hh"

#include "G4AlphaInelasticProcess.hh"
#include "G4LEAlphaInelastic.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// Nuclei
#include "G4IonConstructor.hh"

G4IonPhysics::G4IonPhysics(const G4String& name)
                  :  G4VPhysicsConstructor(name), wasActivated(false)
{
}

G4IonPhysics::~G4IonPhysics()
{
  if(wasActivated) {

    delete theElasticModel;
    delete theIonElasticProcess;
    delete theDElasticProcess;
    delete fDeuteronProcess;
    delete fDeuteronModel;
    delete theTElasticProcess;
    delete fTritonProcess;
    delete fTritonModel;
    delete theAElasticProcess;
    delete fAlphaProcess;
    delete fAlphaModel;
    delete theHe3ElasticProcess;

   }
 }

void G4IonPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;

  // Elastic Process
  theElasticModel = new G4LElastic();

  theIonElasticProcess = new G4HadronElasticProcess();
  theDElasticProcess   = new G4HadronElasticProcess();
  theTElasticProcess   = new G4HadronElasticProcess();
  theAElasticProcess   = new G4HadronElasticProcess();
  theHe3ElasticProcess = new G4HadronElasticProcess();

  theIonElasticProcess->RegisterMe(theElasticModel);
  theDElasticProcess->RegisterMe(theElasticModel);
  theTElasticProcess->RegisterMe(theElasticModel);
  theAElasticProcess->RegisterMe(theElasticModel);
  theHe3ElasticProcess->RegisterMe(theElasticModel);

  // Generic Ion
  pManager = G4GenericIon::GenericIon()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theIonElasticProcess);

  // Deuteron
  pManager = G4Deuteron::Deuteron()->GetProcessManager();
  // add process
  fDeuteronModel = new G4LEDeuteronInelastic();
  fDeuteronProcess = new G4DeuteronInelasticProcess();
  fDeuteronProcess->RegisterMe(fDeuteronModel);
  pManager->AddDiscreteProcess(theDElasticProcess);
  pManager->AddDiscreteProcess(fDeuteronProcess);

  // Triton
  pManager = G4Triton::Triton()->GetProcessManager();
  // add process
  fTritonModel = new G4LETritonInelastic();
  fTritonProcess = new G4TritonInelasticProcess();
  fTritonProcess->RegisterMe(fTritonModel);
  pManager->AddDiscreteProcess(theTElasticProcess);
  pManager->AddDiscreteProcess(fTritonProcess);

  // Alpha
  pManager = G4Alpha::Alpha()->GetProcessManager();
  // add process
  fAlphaModel = new G4LEAlphaInelastic();
  fAlphaProcess = new G4AlphaInelasticProcess();
  fAlphaProcess->RegisterMe(fAlphaModel);
  pManager->AddDiscreteProcess(theAElasticProcess);
  pManager->AddDiscreteProcess(fAlphaProcess);

  // He3
  pManager = G4He3::He3()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(theHe3ElasticProcess);

  wasActivated = true;
}

 void G4IonPhysics::ConstructParticle()
 {
   //  Construct light ions
   G4IonConstructor pConstructor;
   pConstructor.ConstructParticle();  
 }



 // 2002 by J.P. Wellisch
