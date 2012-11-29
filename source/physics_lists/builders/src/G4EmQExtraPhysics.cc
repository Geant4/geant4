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
// ClassName:   G4EmQExtraPhysics
//
// Author: 16-Oct-2012 A. Ribon
//         Copied from the original G4EmExtraPhysics and renamed
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4EmQExtraPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4SynchrotronRadiation.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4ProcessManager.hh"
#include "G4BuilderType.hh"
#include "G4HadronicDeprecate.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmQExtraPhysics);

G4EmQExtraPhysics::G4EmQExtraPhysics(G4int ver): 
  G4VPhysicsConstructor("G4GammaLeptoNuclearPhys"), wasBuilt(false), gnActivated(false), 
  munActivated(false), synActivated(false), synchOn(false), gammNucOn(true), muNucOn(false), 
  theElectronSynch(0), thePositronSynch(0), theGNPhysics(0), muNucProcess(0), muNucModel(0),
  verbose(ver)
{
  G4HadronicDeprecate("G4EmQExtraPhysics");
  theMessenger = new G4EmQMessenger(this);
  SetPhysicsType(bEmExtra);
}

G4EmQExtraPhysics::G4EmQExtraPhysics(const G4String&): 
  G4VPhysicsConstructor("G4GammaLeptoNuclearPhys"), wasBuilt(false), gnActivated(false), 
  munActivated(false), synActivated(false), synchOn(false), gammNucOn(true), muNucOn(false), 
  theElectronSynch(0), thePositronSynch(0), theGNPhysics(0), muNucProcess(0), muNucModel(0),
  verbose(1)
{
  G4HadronicDeprecate("G4EmQExtraPhysics");
  theMessenger = new G4EmQMessenger(this);
  SetPhysicsType(bEmExtra);
}

G4EmQExtraPhysics::~G4EmQExtraPhysics()
{
  delete theMessenger;
  delete theElectronSynch;
  delete thePositronSynch;
  delete theGNPhysics;
  delete muNucProcess;
  delete muNucModel;
}

void G4EmQExtraPhysics::Synch(G4String & newState)
{
  if(newState == "on" || newState == "ON") {
    synchOn = true;
    if(wasBuilt) BuildSynch();
  } else synchOn = false;
}

void G4EmQExtraPhysics::GammaNuclear(G4String & newState)
{
  if(newState == "on" || newState == "ON") {
    gammNucOn = true;
    if(wasBuilt) BuildGammaNuclear();
  } else  gammNucOn = false;
}

void G4EmQExtraPhysics::MuonNuclear(G4String & newState)
{
  if(newState == "on" || newState == "ON") {
    muNucOn = true;
    if(wasBuilt) BuildMuonNuclear();
  } else muNucOn = false;
}

void G4EmQExtraPhysics::ConstructParticle()
{
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();
}

void G4EmQExtraPhysics::ConstructProcess()
{
  if(wasBuilt) return;
  wasBuilt = true;

  if (synchOn)   BuildSynch();
  if (gammNucOn) BuildGammaNuclear();
  if (muNucOn)   BuildMuonNuclear();
}

void G4EmQExtraPhysics::BuildMuonNuclear()
{
  if(munActivated) return;
  munActivated = true;
  G4ProcessManager * pManager = 0;

  muNucProcess = new G4MuonNuclearProcess();
  muNucModel = new G4MuonVDNuclearModel();
  muNucProcess->RegisterMe(muNucModel);

  pManager = G4MuonPlus::MuonPlus()->GetProcessManager();
  pManager->AddDiscreteProcess(muNucProcess);

  pManager = G4MuonMinus::MuonMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(muNucProcess);
}

void G4EmQExtraPhysics::BuildGammaNuclear()
{
  if(gnActivated) return;
  gnActivated = true;

  theGNPhysics = new G4ElectroNuclearBuilder();
  theGNPhysics->Build();
}

void G4EmQExtraPhysics::BuildSynch()
{
  if(synActivated) return;
  synActivated = true;
  G4ProcessManager * pManager = 0;

  pManager = G4Electron::Electron()->GetProcessManager();
  theElectronSynch = new G4SynchrotronRadiation();
  pManager->AddDiscreteProcess(theElectronSynch);

  pManager = G4Positron::Positron()->GetProcessManager();
  thePositronSynch = new G4SynchrotronRadiation();
  pManager->AddDiscreteProcess(thePositronSynch);
}
