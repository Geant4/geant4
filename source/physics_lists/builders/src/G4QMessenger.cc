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
// ClassName:   G4QMessenger
//
// Author: 2009 M. V. Kossov
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4QMessenger.hh"
#include "G4QNeutrinoPhysics.hh"
#include "G4QPhotoNuclearPhysics.hh"

G4QMessenger::G4QMessenger()
{
  rootDir = new G4UIdirectory("/CHIPS_physics/");
  rootDir->SetGuidance("messenger of the CHIPS processes");
  weakDir=0;
  theWeak=0;
  theNuElN=0;
  theNuMuN=0;
  theNuTaN=0;
  biasNuNuc=0;
  photoDir=0;
  thePhoto=0;
  theSynchR=0;
  minGamSR=0;
  theGamN=0;
  theEleN=0;
  theMuoN=0;
  theTauN=0;
  biasPhotoN=0;
}

G4QMessenger::~G4QMessenger()
{
  delete rootDir;
  if(photoDir)
  {
    delete photoDir;
    delete theSynchR;
    delete minGamSR;
    delete theGamN;
    delete theEleN;
    delete theMuoN;
    delete theTauN;
    delete biasPhotoN;
  }
  else if(weakDir)
  {
    delete weakDir;
    delete theNuElN;
    delete theNuMuN;
    delete theNuTaN;
    delete biasNuNuc;
  }
}

// Returns Pointer to the G4QMessenger class
G4QMessenger* G4QMessenger::GetPointer()
{
  static G4QMessenger theMessenger; //**Static body of CHIPS Messenger**
  return &theMessenger;
}

void G4QMessenger::Add(G4QNeutrinoPhysics* weak)
{
  theWeak = weak;

  weakDir = new G4UIdirectory("/CHIPS_physics/neutrino/");
  weakDir->SetGuidance("weak (neutrino) processes");

  // commands for neutrino_el-nuclear physics.
  theNuElN = new G4UIcmdWithAString("/CHIPS_physics/neutrino/NuElNuclear",this);
  theNuElN->SetGuidance("Switching of nu_el-nuclear physics.");
  theNuElN->SetParameterName("status",false);
  theNuElN->SetCandidates("on off");
  theNuElN->SetDefaultValue("off");
  theNuElN->AvailableForStates(G4State_PreInit, G4State_Idle);

  // commands for neutrino_mu-nuclear physics.
  theNuMuN = new G4UIcmdWithAString("/CHIPS_physics/neutrino/NuMuNuclear",this);
  theNuMuN->SetGuidance("Switching of nu_mu-nuclear physics.");
  theNuMuN->SetParameterName("status",false);
  theNuMuN->SetCandidates("on off");
  theNuMuN->SetDefaultValue("off");
  theNuMuN->AvailableForStates(G4State_PreInit, G4State_Idle);

  // commands for neutrino_tau-nuclear physics.
  theNuTaN = new G4UIcmdWithAString("/CHIPS_physics/neutrino/NuTauNuclear",this);
  theNuTaN->SetGuidance("Switching of nu_tau-nuclear physics.");
  theNuTaN->SetParameterName("status",false);
  theNuTaN->SetCandidates("on off");
  theNuTaN->SetDefaultValue("off");
  theNuTaN->AvailableForStates(G4State_PreInit, G4State_Idle);

  // command for biasing of neutrino-nuclear reactions
  biasNuNuc = new G4UIcmdWithADouble("/CHIPS_physics/neutrino/NuNuc_Biasing", this);  
  biasNuNuc->SetGuidance("Set a biasing coefficient for neutrino-nuclear ractions");
  biasNuNuc->SetParameterName("NuNuc_Biasing", false);
  biasNuNuc->SetRange("NuNuc_Biasing > 0.");
  biasNuNuc->SetDefaultValue(1.);
  biasNuNuc->AvailableForStates(G4State_PreInit, G4State_Idle);
}

void G4QMessenger::Add(G4QPhotoNuclearPhysics* photo)
{
  thePhoto = photo;

  weakDir = new G4UIdirectory("/CHIPS_physics/photoNuclear/");
  weakDir->SetGuidance("weak (neutrino) processes");

  // use G4UIcmdWithADouble for weighting of processes

  // command for synchrotron radiation.
  theSynchR = new G4UIcmdWithAString("/CHIPS_physics/photoNuclear/SynchRadiation",this);
  theSynchR->SetGuidance("Switching on/off synchrotron radiation.");
  theSynchR->SetParameterName("status",false);
  theSynchR->SetCandidates("on off");
  theSynchR->SetDefaultValue("off");
  theSynchR->AvailableForStates(G4State_PreInit, G4State_Idle);

  minGamSR = new G4UIcmdWithADouble("/CHIPS_physics/photoNuclear/MinGamma_SynchRad", this);
  minGamSR->SetGuidance("Set a minimum gamma for Synchratron Radiation");
  minGamSR->SetParameterName("MinGamma_SynchRad", false);
  minGamSR->SetRange("MinGamma_SynchRad >> 1.");
  minGamSR->SetDefaultValue(227.);
  minGamSR->AvailableForStates(G4State_PreInit, G4State_Idle);

  // command for gamma-nuclear physics.
  theGamN = new G4UIcmdWithAString("/CHIPS_physics/photoNuclear/GammaNuclear",this);
  theGamN->SetGuidance("Switching of gamma-nuclear physics.");
  theGamN->SetParameterName("status",false);
  theGamN->SetCandidates("on off");
  theGamN->SetDefaultValue("on");
  theGamN->AvailableForStates(G4State_PreInit, G4State_Idle);

  // command for electro-nuclear physics.
  theEleN = new G4UIcmdWithAString("/CHIPS_physics/photoNuclear/ElectroNuclear",this);
  theEleN->SetGuidance("Switching of electron-nuclear physics.");
  theEleN->SetParameterName("status",false);
  theEleN->SetCandidates("on off");
  theEleN->SetDefaultValue("off");
  theEleN->AvailableForStates(G4State_PreInit, G4State_Idle);

  // command for muon-nuclear physics.
  theMuoN = new G4UIcmdWithAString("/CHIPS_physics/photoNuclear/MuonNuclear",this);
  theMuoN->SetGuidance("Switching of muon nuclear physics.");
  theMuoN->SetParameterName("status",false);
  theMuoN->SetCandidates("on off");
  theMuoN->SetDefaultValue("off");
  theMuoN->AvailableForStates(G4State_PreInit, G4State_Idle);

  // command for tau-nuclear physics.
  theTauN = new G4UIcmdWithAString("/CHIPS_physics/photoNuclear/TauNuclear",this);
  theTauN->SetGuidance("Switching of tau nuclear physics.");
  theTauN->SetParameterName("status",false);
  theTauN->SetCandidates("on off");
  theTauN->SetDefaultValue("off");
  theTauN->AvailableForStates(G4State_PreInit, G4State_Idle);


  biasPhotoN = new G4UIcmdWithADouble("/CHIPS_physics/photoNuclear/PhotoN_Biasing", this);
  biasPhotoN->SetGuidance("Set a biasing coefficient for photo-nuclear ractions");
  biasPhotoN->SetParameterName("PhotoN_Biasing", false);
  biasPhotoN->SetRange("PhotoN_Biasing > 0.");
  biasPhotoN->SetDefaultValue(1.);
  biasPhotoN->AvailableForStates(G4State_PreInit, G4State_Idle);
}

void G4QMessenger::SetNewValue(G4UIcommand* aComm, G4String aS)
{
  if(photoDir)
  {
    if     (aComm==theSynchR) thePhoto->SetSynchRadOnOff(aS);
    else if(aComm==minGamSR)  thePhoto->SetMinGammaSR(minGamSR->GetNewDoubleValue(aS));
    else if(aComm==theGamN)   thePhoto->SetGammaNuclearOnOff(aS);
    else if(aComm==theMuoN)   thePhoto->SetElPosNuclearOnOff(aS);
    else if(aComm==theMuoN)   thePhoto->SetMuonNuclearOnOff(aS);
    else if(aComm==theMuoN)   thePhoto->SetTauNuclearOnOff(aS);
    else if(aComm==biasPhotoN)thePhoto->SetPhotoNucBias(biasPhotoN->GetNewDoubleValue(aS));
  }
  else if(weakDir)
  {
    if     (aComm==theNuElN)  theWeak->SetNuElNuclearOnOff(aS);
    else if(aComm==theNuMuN)  theWeak->SetNuMuNuclearOnOff(aS);
    else if(aComm==theNuTaN)  theWeak->SetNuTauNuclearOnOff(aS);
    else if(aComm==biasNuNuc) theWeak->SetNuNuclearBias(biasNuNuc->GetNewDoubleValue(aS));
  }
}

G4String G4QMessenger::GetCurrentValue(G4UIcommand* aComm)
{
  if(photoDir)
  {
    if     (aComm==theSynchR) return thePhoto->GetSynchRadOnOff();
    else if(aComm==theGamN)   return thePhoto->GetGammaNuclearOnOff();
    else if(aComm==theEleN)   return thePhoto->GetElPosNuclearOnOff();
    else if(aComm==theMuoN)   return thePhoto->GetMuonNuclearOnOff();
    else if(aComm==theTauN)   return thePhoto->GetTauNuclearOnOff();
  }
  else if(weakDir)
  {
    if     (aComm==theNuElN)  return theWeak->GetNuElNuclearOnOff();
    else if(aComm==theNuMuN)  return theWeak->GetNuMuNuclearOnOff();
    else if(aComm==theNuTaN)  return theWeak->GetNuTauNuclearOnOff();
  }
  return "not_defined";
}
