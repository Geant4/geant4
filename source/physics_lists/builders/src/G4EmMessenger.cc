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
// ClassName:   G4EmMessenger
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 09.11.2005 V.Ivanchenko edit to provide a standard
// 19.06.2006 V.Ivanchenko add mu-nuclear process
//
//----------------------------------------------------------------------------
//

#include "G4EmMessenger.hh"
#include "G4EmExtraPhysics.hh"

G4EmMessenger::G4EmMessenger(G4EmExtraPhysics* ab)
{
  theB = ab;
  aDir1 = new G4UIdirectory("/physics_engine/");
  aDir1->SetGuidance("commands related to the physics simulation engine.");

  // general stuff.
  aDir2 = new G4UIdirectory("/physics_engine/tailor/");
  aDir2->SetGuidance("tailoring the processes");

  // command for synchrotron radiation.
  theSynch = new G4UIcmdWithAString("/physics_engine/tailor/SyncRadiation",this);
  theSynch->SetGuidance("Switching on/off synchrotron radiation.");
  theSynch->SetParameterName("status",false);
  theSynch->SetCandidates("on off");
  theSynch->SetDefaultValue("off");
  theSynch->AvailableForStates(G4State_PreInit,G4State_Idle);

  // command for gamma nuclear physics.
  theGN = new G4UIcmdWithAString("/physics_engine/tailor/GammaNuclear",this);
  theGN->SetGuidance("Switching on gamma nuclear physics.");
  theGN->SetParameterName("status",false);
  theGN->SetCandidates("on off");
  theGN->SetDefaultValue("on");
  theGN->AvailableForStates(G4State_PreInit,G4State_Idle);

  // command for muon nuclear physics.
  theMUN = new G4UIcmdWithAString("/physics_engine/tailor/MuonNuclear",this);
  theMUN->SetGuidance("Switching on muon nuclear physics.");
  theMUN->SetParameterName("status",false);
  theMUN->SetCandidates("on off");
  theMUN->SetDefaultValue("off");
  theMUN->AvailableForStates(G4State_PreInit,G4State_Idle);
}

G4EmMessenger::~G4EmMessenger()
{
  delete theSynch;
  delete theGN;
  delete theMUN;
  delete aDir1;
  delete aDir2;
}

void G4EmMessenger::SetNewValue(G4UIcommand* aComm, G4String aS)
{
  if(aComm==theSynch) theB->Synch(aS);
  if(aComm==theGN)    theB->GammaNuclear(aS);
  if(aComm==theMUN)   theB->MuonNuclear(aS);
}
