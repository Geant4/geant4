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
// ClassName:   G4QNeutrinoPhysics
//
// Author: 2009 M. V. Kosov
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4QNeutrinoPhysics.hh"

G4QNeutrinoPhysics::G4QNeutrinoPhysics(G4int): 
  G4VPhysicsConstructor("CHIPS neutrino-nuclear"), wasBuilt(false), nuEleActivated(false),
  nuMuoActivated(false), nuTauActivated(false), nuEleOn(false),
  nuMuoOn(false), nuTauOn(false), nuNucBias(1.), inelastic(0)
{
  theMessenger = G4QMessenger::GetPointer();
  theMessenger->Add(this);
}

G4QNeutrinoPhysics::G4QNeutrinoPhysics(const G4String& name): 
  G4VPhysicsConstructor(name), wasBuilt(false), nuEleActivated(false),
  nuMuoActivated(false), nuTauActivated(false), nuEleOn(false),
  nuMuoOn(false), nuTauOn(false), nuNucBias(1.), inelastic(0)
{
  theMessenger = G4QMessenger::GetPointer();
  theMessenger->Add(this);
}

G4QNeutrinoPhysics::~G4QNeutrinoPhysics()
{
  if(wasBuilt && inelastic) delete inelastic;
}

void G4QNeutrinoPhysics::ConstructParticle()
{
  G4NeutrinoE::NeutrinoE();
  G4AntiNeutrinoE::AntiNeutrinoE();
  G4NeutrinoMu::NeutrinoMu();
  G4AntiNeutrinoMu::AntiNeutrinoMu();
  G4NeutrinoTau::NeutrinoTau();
  G4AntiNeutrinoTau::AntiNeutrinoTau();
}

void G4QNeutrinoPhysics::SetNuElNuclearOnOff(G4String& newSwitch)
{
  if(wasBuilt) G4cout<<"G4QNeutrinoPhysics:No, processes are already builded!"<<G4endl;
  else if(newSwitch == "on" || newSwitch == "ON" || newSwitch == "On") nuEleOn = true;
  else nuEleOn = false;
}

void G4QNeutrinoPhysics::SetNuMuNuclearOnOff(G4String& newSwitch)
{
  if(wasBuilt) G4cout<<"G4QNeutrinoPhysics:No, processes are already builded!"<<G4endl;
  else if(newSwitch == "on" || newSwitch == "ON" || newSwitch == "On") nuMuoOn = true;
  else nuMuoOn = false;
}

void G4QNeutrinoPhysics::SetNuTauNuclearOnOff(G4String& newSwitch)
{
  if(wasBuilt) G4cout<<"G4QNeutrinoPhysics:No, processes are already builded!"<<G4endl;
  else if(newSwitch == "on" || newSwitch == "ON" || newSwitch == "On") nuTauOn = true;
  else nuTauOn = false;
}

void G4QNeutrinoPhysics::SetNuNuclearBias(G4double newValue)
{
  if(wasBuilt) G4cout<<"G4QNeutrinoPhysics:No, processes are already builded!"<<G4endl;
  else nuNucBias = newValue;
}

void G4QNeutrinoPhysics::ConstructProcess()
{
  if(wasBuilt) return;
  if(nuEleOn || nuMuoOn || nuTauOn)
  {
    wasBuilt = true;
    G4cout<<"Builded=>G4QNeutrinoPhysics: "<<nuEleOn<<", "<<nuMuoOn<<", "<<nuTauOn<<G4endl;
    inelastic = new G4QInelastic("neutrinoNuclear");
    inelastic->SetWeakNucBias(nuNucBias); // enough only once (static)
    if (nuEleOn)   BuildNuEleNuclear();
    if (nuMuoOn)   BuildNuMuoNuclear();
    if (nuTauOn)   BuildNuTauNuclear();
  }
}

void G4QNeutrinoPhysics::BuildNuEleNuclear()
{
  if(nuEleActivated) return;
  nuEleActivated = true;
  G4ProcessManager * pManager = 0;

  pManager  = G4NeutrinoE::NeutrinoE()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);

  pManager  = G4AntiNeutrinoE::AntiNeutrinoE()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);
}

void G4QNeutrinoPhysics::BuildNuMuoNuclear()
{
  if(nuMuoActivated) return;
  nuMuoActivated = true;
  G4ProcessManager * pManager = 0;

  pManager  = G4NeutrinoMu::NeutrinoMu()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);

  pManager  = G4AntiNeutrinoMu::AntiNeutrinoMu()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);
}

void G4QNeutrinoPhysics::BuildNuTauNuclear()
{
  if(nuTauActivated) return;
  nuTauActivated = true;
  G4ProcessManager * pManager = 0;

  pManager  = G4NeutrinoTau::NeutrinoTau()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);

  pManager  = G4AntiNeutrinoTau::AntiNeutrinoTau()->GetProcessManager();
  pManager->AddDiscreteProcess(inelastic);
}
