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
// ClassName:   G4QEmExtraPhysics
//
// Author: 23 May 2007 V. Ivanchenko
//
// Modified: 19 Nov 2009 M.Kosov: G4QInelastic instead of G4QCollision
//
//----------------------------------------------------------------------------
//

#include "G4QEmExtraPhysics.hh"
#include "G4QInelastic.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4LeptonConstructor.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4HadronicDeprecate.hh"


G4QEmExtraPhysics::G4QEmExtraPhysics(G4int ver)
  :  G4VPhysicsConstructor("QEmExtra"), hProcess(0), verbose(ver), wasActivated(false) 
{
  G4HadronicDeprecate("G4QEmExtraPhysics");
  if(verbose > 1) G4cout << "### G4QEmExtraPhysics" << G4endl;
}

G4QEmExtraPhysics::~G4QEmExtraPhysics()
{
  delete hProcess;
}

void G4QEmExtraPhysics::ConstructParticle()
{
// G4cout << "G4QEmExtraPhysics::ConstructParticle" << G4endl;
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();
}

void G4QEmExtraPhysics::ConstructProcess()
{
  if(verbose > 1) G4cout << "### G4QEmExtraPhysics::ConstructProcess " 
			 << wasActivated << G4endl;
  if(wasActivated) return;
  wasActivated = true;

  hProcess = new G4QInelastic();

  G4ParticleDefinition* particle = G4Gamma::Gamma();
  G4ProcessManager* pmanager = particle->GetProcessManager();
  pmanager->AddDiscreteProcess(hProcess);

  particle = G4Electron::Electron();
  pmanager = particle->GetProcessManager();
  pmanager->AddDiscreteProcess(hProcess);
  
  particle = G4Positron::Positron();
  pmanager = particle->GetProcessManager();
  pmanager->AddDiscreteProcess(hProcess);

  particle = G4MuonMinus::MuonMinus();
  pmanager = particle->GetProcessManager();
  pmanager->AddDiscreteProcess(hProcess);

  particle = G4MuonPlus::MuonPlus();
  pmanager = particle->GetProcessManager();
  pmanager->AddDiscreteProcess(hProcess);
  
}


