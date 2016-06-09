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
// $Id: ExN04IonPhysics.cc,v 1.4 2006/06/29 17:37:13 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 

#include "ExN04IonPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   


ExN04IonPhysics::ExN04IonPhysics(const G4String& name)
                 :  G4VPhysicsConstructor(name)
{
}

ExN04IonPhysics::~ExN04IonPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

// Nuclei
#include "G4IonConstructor.hh"

void ExN04IonPhysics::ConstructParticle()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}


#include "G4ProcessManager.hh"


void ExN04IonPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;
  
  // Elastic Process
  theElasticModel = new G4LElastic();
  theElasticProcess.RegisterMe(theElasticModel);

  // Generic Ion
  pManager = G4GenericIon::GenericIon()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  pManager->AddProcess(&fIonMultipleScattering, -1, 1, 1);
  pManager->AddProcess(&fIonIonisation,         -1, 2, 2);

  // Deuteron 
  pManager = G4Deuteron::Deuteron()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  fDeuteronModel = new G4LEDeuteronInelastic();
  fDeuteronProcess.RegisterMe(fDeuteronModel);
  pManager->AddDiscreteProcess(&fDeuteronProcess);

  pManager->AddProcess(&fDeuteronMultipleScattering, -1, 1, 1);
  pManager->AddProcess(&fDeuteronIonisation,        -1, 2, 2);

  // Triton 
  pManager = G4Triton::Triton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  fTritonModel = new G4LETritonInelastic();
  fTritonProcess.RegisterMe(fTritonModel);
  pManager->AddDiscreteProcess(&fTritonProcess);

  pManager->AddProcess(&fTritonMultipleScattering, -1, 1, 1);
  pManager->AddProcess(&fTritonIonisation,        -1, 2, 2);
 
  // Alpha 
  pManager = G4Alpha::Alpha()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  fAlphaModel = new G4LEAlphaInelastic();
  fAlphaProcess.RegisterMe(fAlphaModel);
  pManager->AddDiscreteProcess(&fAlphaProcess);

  pManager->AddProcess(&fAlphaMultipleScattering, -1, 1, 1);
  pManager->AddProcess(&fAlphaIonisation,        -1, 2, 2);

  // He3
  pManager = G4He3::He3()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  pManager->AddProcess(&fHe3MultipleScattering, -1, 1, 1);
  pManager->AddProcess(&fHe3Ionisation,        -1, 2, 2);

   
}



