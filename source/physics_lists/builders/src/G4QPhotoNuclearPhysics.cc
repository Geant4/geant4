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
// ClassName:   G4QPhotoNuclearPhysics
//
// Author: 2009 M. V. Kosov
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4QPhotoNuclearPhysics.hh"


G4QPhotoNuclearPhysics::G4QPhotoNuclearPhysics(G4int): 
  G4VPhysicsConstructor("CHIPS photo-nuclear"), wasBuilt(false), SynchRActivated(false),
  GamNucActivated(false), EleNucActivated(false), MuoNucActivated(false),
  TauNucActivated(false), synchrOn(true), synchrMinGam(227.), gamNucOn(true),
  eleNucOn(true), muoNucOn(true), tauNucOn(true), photoNucBias(1.)
{
  theMessenger = G4QMessenger::GetPointer();
  theMessenger->Add(this);
}

G4QPhotoNuclearPhysics::G4QPhotoNuclearPhysics(const G4String& name): 
  G4VPhysicsConstructor(name), wasBuilt(false), SynchRActivated(false),
  GamNucActivated(false), EleNucActivated(false), MuoNucActivated(false),
  TauNucActivated(false), synchrOn(true), synchrMinGam(227.), gamNucOn(true),
  eleNucOn(true), muoNucOn(true), tauNucOn(true), photoNucBias(1.)
{
  theMessenger = G4QMessenger::GetPointer();
  theMessenger->Add(this);
}

G4QPhotoNuclearPhysics::~G4QPhotoNuclearPhysics()
{
  if(wasBuilt)
  {
    delete inelastic;
    if(synchrOn) delete synchrad;
  }
}

void G4QPhotoNuclearPhysics::SetSynchRadOnOff(G4String& newSwitch)
{
  if(wasBuilt) G4cout<<"G4QPhotoNuclearPhysics:No, processes are already builded!"<<G4endl;
  else if(newSwitch == "on" || newSwitch == "ON" || newSwitch == "On") synchrOn = true;
  else synchrOn = false;
}

void G4QPhotoNuclearPhysics::SetGammaNuclearOnOff(G4String& newSwitch)
{
  if(wasBuilt) G4cout<<"G4QPhotoNuclearPhysics:No, processes are already builded!"<<G4endl;
  else if(newSwitch == "on" || newSwitch == "ON" || newSwitch == "On") gamNucOn = true;
  else gamNucOn = false;
}

void G4QPhotoNuclearPhysics::SetElPosNuclearOnOff(G4String& newSwitch)
{
  if(wasBuilt) G4cout<<"G4QPhotoNuclearPhysics:No, processes are already builded!"<<G4endl;
  else if(newSwitch == "on" || newSwitch == "ON" || newSwitch == "On") eleNucOn = true;
  else eleNucOn = false;
}

void G4QPhotoNuclearPhysics::SetMuonNuclearOnOff(G4String& newSwitch)
{
  if(wasBuilt) G4cout<<"G4QPhotoNuclearPhysics:No, processes are already builded!"<<G4endl;
  else if(newSwitch == "on" || newSwitch == "ON" || newSwitch == "On") muoNucOn = true;
  else muoNucOn = false;
}

void G4QPhotoNuclearPhysics::SetTauNuclearOnOff(G4String& newSwitch)
{
  if(wasBuilt) G4cout<<"G4QPhotoNuclearPhysics:No, processes are already builded!"<<G4endl;
  else if(newSwitch == "on" || newSwitch == "ON" || newSwitch == "On") tauNucOn = true;
  else tauNucOn = false;
}

void G4QPhotoNuclearPhysics::SetMinGammaSR(G4double newValue)
{
  if(wasBuilt) G4cout<<"G4QPhotoNuclearPhysics:No, processes are already builded!"<<G4endl;
  else synchrMinGam = newValue;
}

void G4QPhotoNuclearPhysics::SetPhotoNucBias(G4double newValue)
{
  if(wasBuilt) G4cout<<"G4QPhotoNuclearPhysics:No, processes are already builded!"<<G4endl;
  else photoNucBias = newValue;
}

void G4QPhotoNuclearPhysics::ConstructParticle()
{
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();
  G4TauPlus::TauPlus();
  G4TauMinus::TauMinus();
  if (synchrOn)
  {
    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();
  }
}

void G4QPhotoNuclearPhysics::ConstructProcess()
{
  if(wasBuilt) return;
  wasBuilt = true;

  inelastic = new G4QInelastic("photoNuclear");
  inelastic->SetPhotNucBias(photoNucBias);

  if (synchrOn)   BuildSynchRad();
  if (gamNucOn)   BuildGammaNuclear();
  if (eleNucOn)   BuildElectroNuclear();
  if (muoNucOn)   BuildMuonNuclear();
  if (tauNucOn)   BuildTauNuclear();
}

void G4QPhotoNuclearPhysics::BuildGammaNuclear()
{
  if(GamNucActivated) return;
  GamNucActivated = true;
  G4ProcessManager* pManager  = G4Gamma::Gamma()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);
}

void G4QPhotoNuclearPhysics::BuildElectroNuclear()
{
  if(EleNucActivated) return;
  EleNucActivated = true;
  G4ProcessManager * pManager = 0;

  pManager  = G4Electron::Electron()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);

  pManager  = G4Positron::Positron()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);
}

void G4QPhotoNuclearPhysics::BuildMuonNuclear()
{
  if(MuoNucActivated) return;
  MuoNucActivated = true;
  G4ProcessManager * pManager = 0;

  pManager  = G4MuonPlus::MuonPlus()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);

  pManager  = G4MuonMinus::MuonMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);
}

void G4QPhotoNuclearPhysics::BuildTauNuclear()
{
  if(TauNucActivated) return;
  TauNucActivated = true;
  G4ProcessManager * pManager = 0;

  pManager  = G4TauPlus::TauPlus()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);

  pManager  = G4TauMinus::TauMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);
}

// The CHIPS Synchrotron radiation process is working for all charged particles
void G4QPhotoNuclearPhysics::BuildSynchRad()
{
  if(SynchRActivated) return;
  SynchRActivated = true;
  synchrad = new G4QSynchRad();
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4double charge = particle->GetPDGCharge();
    if(charge != 0.0)
    {
      G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager->AddDiscreteProcess(synchrad);
    }
  }
}
