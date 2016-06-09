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
// ClassName:   G4QIonPhysics
//
// Author:  03/06/2010, M. Kossov
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4QIonPhysics.hh"
#include "G4BuilderType.hh"

G4QIonPhysics::G4QIonPhysics(const G4String& name)
 : G4VPhysicsConstructor(name)
    , fQAAInelasticProcess(0),fQAAElasticProcess(0)
    , wasActivated(false)
{
  SetPhysicsType(bIons);
}

void G4QIonPhysics::ConstructProcess()
{
  G4ProcessManager* pManager = 0;
  fQAAInelasticProcess = new G4QLowEnergy();
  fQAAElasticProcess   = new G4QIonIonElastic();

 // Deuteron
  pManager = G4Deuteron::Deuteron()->GetProcessManager();
  pManager->AddDiscreteProcess(fQAAInelasticProcess);
  pManager->AddDiscreteProcess(fQAAElasticProcess);

  // Triton
  pManager = G4Triton::Triton()->GetProcessManager();
  pManager->AddDiscreteProcess(fQAAInelasticProcess);
  pManager->AddDiscreteProcess(fQAAElasticProcess);

  // He3
  pManager = G4He3::He3()->GetProcessManager();
  pManager->AddDiscreteProcess(fQAAInelasticProcess);
  pManager->AddDiscreteProcess(fQAAElasticProcess);

  // Alpha
  pManager = G4Alpha::Alpha()->GetProcessManager();
  pManager->AddDiscreteProcess(fQAAInelasticProcess);
  pManager->AddDiscreteProcess(fQAAElasticProcess);

  // GenericIon
  pManager = G4GenericIon::GenericIon()->GetProcessManager();
  pManager->AddDiscreteProcess(fQAAInelasticProcess);
  pManager->AddDiscreteProcess(fQAAElasticProcess);

  wasActivated = true;
}

 void G4QIonPhysics::ConstructParticle()
 {
   //  Construct light ions
   G4IonConstructor pConstructor;
   pConstructor.ConstructParticle();  
 }
