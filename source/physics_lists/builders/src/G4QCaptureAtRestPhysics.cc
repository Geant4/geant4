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
// ClassName:   G4QCaptureAtRestPhysics
//
// Author: 16 Nov 2009 M. Kosov
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4QCaptureAtRestPhysics.hh"

G4QCaptureAtRestPhysics::G4QCaptureAtRestPhysics(G4int ver)
  :  G4VPhysicsConstructor("CHIPS nuclear_capture"), captureProcess(0)
   , verbose(ver), wasActivated(false)
{
  if(verbose > 1) G4cout << "###> G4QCaptureAtRestPhysics is initialized" << G4endl;
}

G4QCaptureAtRestPhysics::G4QCaptureAtRestPhysics(const G4String& name, G4int ver)
  :  G4VPhysicsConstructor(name), captureProcess(0), verbose(ver), wasActivated(false)
{
  if(verbose > 1) G4cout << "###> G4QCaptureAtRestPhysics is initialized" << G4endl;
}

G4QCaptureAtRestPhysics::~G4QCaptureAtRestPhysics()
{
  if(wasActivated) delete captureProcess;
}

void G4QCaptureAtRestPhysics::ConstructParticle()
{
// G4cout << "G4QCaptureAtRestPhysics::ConstructParticle" << G4endl;
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

}

void G4QCaptureAtRestPhysics::ConstructProcess()
{
  if(verbose>1) G4cout<<"###> G4QCaptureAtRestPhysics::ConstructProcess: "<<wasActivated
                      <<G4endl;
  if(wasActivated) return;
  wasActivated = true;
  captureProcess = new G4QCaptureAtRest();

  // Add Stopping Process
  G4ParticleDefinition* particle=0;
  G4ProcessManager* pmanager=0;

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    particle = theParticleIterator->value();
    pmanager = particle->GetProcessManager();
    if(particle->GetPDGCharge() < 0. && particle != G4Electron::Electron() &&
       !(particle->IsShortLived()) && captureProcess->IsApplicable(*particle) ) 
    { 
      pmanager->AddRestProcess(captureProcess);
      if(verbose > 1) G4cout << "###> G4QCaptureAtRestPhysics is added for " 
                             <<particle->GetParticleName() << G4endl;
    }
  }
}


