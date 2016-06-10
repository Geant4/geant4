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
// $Id: G4EmExtraPhysics.cc 66704 2013-01-10 18:20:17Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmExtraPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard
// 19.06.2006 V.Ivanchenko add mu-nuclear process
// 16.10.2012 A.Ribon: renamed G4EmExtraBertiniPhysics as G4EmExtraPhysics
//
//----------------------------------------------------------------------------
//

#include "G4EmExtraPhysics.hh"

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

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmExtraPhysics);

G4ThreadLocal G4bool G4EmExtraPhysics::wasBuilt     = false;
G4ThreadLocal G4bool G4EmExtraPhysics::munActivated = false;
G4ThreadLocal G4bool G4EmExtraPhysics::gnActivated  = false;
G4ThreadLocal G4bool G4EmExtraPhysics::synActivated = false;
G4ThreadLocal G4bool G4EmExtraPhysics::synchOn      = false;
G4ThreadLocal G4bool G4EmExtraPhysics::gammNucOn    = true;
G4ThreadLocal G4bool G4EmExtraPhysics::muNucOn      = false;


G4EmExtraPhysics::G4EmExtraPhysics(G4int ver): 
  G4VPhysicsConstructor("G4GammaLeptoNuclearPhys"),
  verbose(ver)
{
  theMessenger = new G4EmMessenger(this);
  SetPhysicsType(bEmExtra);
  if(verbose > 1) G4cout << "### G4EmExtraPhysics" << G4endl;
}

G4EmExtraPhysics::G4EmExtraPhysics(const G4String&): 
  G4VPhysicsConstructor("G4GammaLeptoNuclearPhys"),
  verbose(1)
{
  theMessenger = new G4EmMessenger(this);
  SetPhysicsType(bEmExtra);
  if(verbose > 1) G4cout << "### G4EmExtraPhysics" << G4endl;
}

G4EmExtraPhysics::~G4EmExtraPhysics()
{
  delete theMessenger;
}

void G4EmExtraPhysics::Synch(G4String & newState)
{
  if(newState == "on" || newState == "ON") {
    synchOn = true;
    if(wasBuilt) BuildSynch();
  } else synchOn = false;
}

void G4EmExtraPhysics::GammaNuclear(G4String & newState)
{
  if(newState == "on" || newState == "ON") {
    gammNucOn = true;
    if(wasBuilt) BuildGammaNuclear();
  } else  gammNucOn = false;
}

void G4EmExtraPhysics::MuonNuclear(G4String & newState)
{
  if(newState == "on" || newState == "ON") {
    muNucOn = true;
    if(wasBuilt) BuildMuonNuclear();
  } else muNucOn = false;
}

void G4EmExtraPhysics::ConstructParticle()
{
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();
}

void G4EmExtraPhysics::ConstructProcess()
{
  if(wasBuilt) return;
  wasBuilt = true;

  if (synchOn)   BuildSynch();
  if (gammNucOn) BuildGammaNuclear();
  if (muNucOn)   BuildMuonNuclear();
}

void G4EmExtraPhysics::BuildMuonNuclear()
{
  if(munActivated) return;
  munActivated = true;
  G4ProcessManager * pManager = 0;

  G4MuonNuclearProcess* muNucProcess = new G4MuonNuclearProcess();
  G4MuonVDNuclearModel* muNucModel = new G4MuonVDNuclearModel();
  muNucProcess->RegisterMe(muNucModel);

  pManager = G4MuonPlus::MuonPlus()->GetProcessManager();
  pManager->AddDiscreteProcess(muNucProcess);

  pManager = G4MuonMinus::MuonMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(muNucProcess);
}

void G4EmExtraPhysics::BuildGammaNuclear()
{
  if(gnActivated) return;
  gnActivated = true;

  G4BertiniElectroNuclearBuilder* theGNPhysics = new G4BertiniElectroNuclearBuilder();
  theGNPhysics->Build();
}

void G4EmExtraPhysics::BuildSynch()
{
  if(synActivated) return;
  synActivated = true;
  G4ProcessManager * pManager = 0;

  pManager = G4Electron::Electron()->GetProcessManager();
  G4SynchrotronRadiation* theElectronSynch = new G4SynchrotronRadiation();
  pManager->AddDiscreteProcess(theElectronSynch);

  pManager = G4Positron::Positron()->GetProcessManager();
  G4SynchrotronRadiation* thePositronSynch = new G4SynchrotronRadiation();
  pManager->AddDiscreteProcess(thePositronSynch);
}
