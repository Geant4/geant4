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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonLHEPPhysics
//
// Author:      A.Ribon 16-Oct-2012
//              Copied from the original G4IonPhysics and renamed.
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4IonLHEPPhysics.hh"

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
#include "G4BuilderType.hh"

#include "G4HadronicDeprecate.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4IonLHEPPhysics);

G4IonLHEPPhysics::G4IonLHEPPhysics(G4int)
                  :  G4VPhysicsConstructor("IonPhysics")
		   , wasActivated(false)
{
  G4HadronicDeprecate("G4IonLHEPPhysics");
  SetPhysicsType(bIons);
  fDeuteronProcess = 0;
  fDeuteronModel = 0;
  fTritonProcess = 0;
  fTritonModel = 0;
  fAlphaProcess = 0;
  fAlphaModel = 0;
}

G4IonLHEPPhysics::G4IonLHEPPhysics(const G4String& name)
                  :  G4VPhysicsConstructor(name), wasActivated(false)
{
  G4HadronicDeprecate("G4IonLHEPPhysics");
  SetPhysicsType(bIons);
  fDeuteronProcess = 0;
  fDeuteronModel = 0;
  fTritonProcess = 0;
  fTritonModel = 0;
  fAlphaProcess = 0;
  fAlphaModel = 0;
}

G4IonLHEPPhysics::~G4IonLHEPPhysics()
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

void G4IonLHEPPhysics::ConstructProcess()
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

 void G4IonLHEPPhysics::ConstructParticle()
 {
   //  Construct light ions
   G4IonConstructor pConstructor;
   pConstructor.ConstructParticle();  
 }
