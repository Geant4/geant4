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
// $Id: G4IonPhysics.cc,v 1.4 2006-06-16 17:41:31 ribon Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

    delete fDeuteronProcess;
    delete fDeuteronModel;
    delete fTritonProcess;
    delete fTritonModel;
    delete fAlphaProcess;
    delete fAlphaModel;

   }
 }

void G4IonPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;

 // Deuteron
  pManager = G4Deuteron::Deuteron()->GetProcessManager();
  // add process
  fDeuteronModel = new G4LEDeuteronInelastic();
  fDeuteronProcess = new G4DeuteronInelasticProcess();
  fDeuteronProcess->RegisterMe(fDeuteronModel);
  pManager->AddDiscreteProcess(fDeuteronProcess);

  // Triton
  pManager = G4Triton::Triton()->GetProcessManager();
  // add process
  fTritonModel = new G4LETritonInelastic();
  fTritonProcess = new G4TritonInelasticProcess();
  fTritonProcess->RegisterMe(fTritonModel);
  pManager->AddDiscreteProcess(fTritonProcess);

  // Alpha
  pManager = G4Alpha::Alpha()->GetProcessManager();
  // add process
  fAlphaModel = new G4LEAlphaInelastic();
  fAlphaProcess = new G4AlphaInelasticProcess();
  fAlphaProcess->RegisterMe(fAlphaModel);
  pManager->AddDiscreteProcess(fAlphaProcess);

  wasActivated = true;
}

 void G4IonPhysics::ConstructParticle()
 {
   //  Construct light ions
   G4IonConstructor pConstructor;
   pConstructor.ConstructParticle();  
 }



 // 2002 by J.P. Wellisch
