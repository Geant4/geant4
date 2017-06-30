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
// $Id: G4EmMessenger.cc 66704 2013-01-10 18:20:17Z gunter $
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


//A. Dotti (8Jun2013): This class does not need changes for MT
// Note that in general "physics" realated commands should not 
// be executed by threads, but this is a special case. Actually the command
// executes a building of processes if it was not build before, thus we need 
// all threads to process commands.
// The logic of thread-private objects is in G4EmExtraPhysics class

G4EmMessenger::G4EmMessenger(G4EmExtraPhysics* ab)
{
  theB = ab;
  aDir1 = new G4UIdirectory("/physics_lists/");
  aDir1->SetGuidance("commands related to the physics simulation engine.");

  // general stuff.
  aDir2 = new G4UIdirectory("/physics_lists/em/");
  aDir2->SetGuidance("tailoring the processes");

  // command for synchrotron radiation.
  theSynch = new G4UIcmdWithABool("/physics_lists/em/SyncRadiation",this);
  theSynch->SetGuidance("Switching on/off synchrotron radiation.");
  theSynch->AvailableForStates(G4State_PreInit);

  // command for synchrotron radiation.
  theSynchAll = new G4UIcmdWithABool("/physics_lists/em/SyncRadiationAll",this);
  theSynchAll->SetGuidance("Switching on/off synchrotron radiation for all charged.");
  theSynchAll->AvailableForStates(G4State_PreInit);

  // command for gamma nuclear physics.
  theGN = new G4UIcmdWithABool("/physics_lists/em/GammaNuclear",this);
  theGN->SetGuidance("Switching on gamma nuclear physics.");
  theGN->AvailableForStates(G4State_PreInit);

  theEN = new G4UIcmdWithABool("/physics_lists/em/ElectroNuclear",this);
  theEN->SetGuidance("Switching on e+- nuclear physics.");
  theEN->AvailableForStates(G4State_PreInit);

  // command for muon nuclear physics.
  theMUN = new G4UIcmdWithABool("/physics_lists/em/MuonNuclear",this);
  theMUN->SetGuidance("Switching on muon nuclear physics.");
  theMUN->AvailableForStates(G4State_PreInit);

  theGMM = new G4UIcmdWithABool("/physics_lists/em/GammaToMuons",this);
  theGMM->SetGuidance("Switching on gamma conversion to muon pair.");
  theGMM->AvailableForStates(G4State_PreInit);

  thePMM = new G4UIcmdWithABool("/physics_lists/em/PositronToMuons",this);
  thePMM->SetGuidance("Switching on positron conversion to muon pair.");
  thePMM->AvailableForStates(G4State_PreInit);

  thePH = new G4UIcmdWithABool("/physics_lists/em/PositronToHadrons",this);
  thePH->SetGuidance("Switching on positron conversion to hadrons.");
  thePH->AvailableForStates(G4State_PreInit);

  theGMM1 = new G4UIcmdWithADouble("/physics_lists/em/GammaToMuonsFactor",this);
  theGMM1->SetGuidance("Factor for gamma conversion to muon pair.");
  theGMM1->AvailableForStates(G4State_PreInit);

  thePMM1 = new G4UIcmdWithADouble("/physics_lists/em/PositronToMuonsFactor",this);
  thePMM1->SetGuidance("Factor for positron conversion to muon pair.");
  thePMM1->AvailableForStates(G4State_PreInit);

  thePH1 = new G4UIcmdWithADouble("/physics_lists/em/PositronToHadronsFactor",this);
  thePH1->SetGuidance("Factor for positron conversion to hadrons.");
  thePH1->AvailableForStates(G4State_PreInit);
}

G4EmMessenger::~G4EmMessenger()
{
  delete theSynch;
  delete theSynchAll;
  delete theGN;
  delete theEN;
  delete theMUN;
  delete theGMM;
  delete thePMM;
  delete thePH;
  delete theGMM1;
  delete thePMM1;
  delete thePH1;
  delete aDir1;
  delete aDir2;
}

void G4EmMessenger::SetNewValue(G4UIcommand* aComm, G4String aS)
{
  if(aComm==theSynch)    theB->Synch(theSynch->GetNewBoolValue(aS));
  if(aComm==theSynchAll) theB->SynchAll(theSynchAll->GetNewBoolValue(aS));
  if(aComm==theGN)       theB->GammaNuclear(theGN->GetNewBoolValue(aS));
  if(aComm==theEN)       theB->ElectroNuclear(theEN->GetNewBoolValue(aS));
  if(aComm==theMUN)      theB->MuonNuclear(theMUN->GetNewBoolValue(aS));
  if(aComm==theGMM)      theB->GammaToMuMu(theGMM->GetNewBoolValue(aS));
  if(aComm==thePMM)      theB->PositronToMuMu(thePMM->GetNewBoolValue(aS));
  if(aComm==thePH)       theB->PositronToHadrons(thePH->GetNewBoolValue(aS));
  if(aComm==theGMM1)     theB->GammaToMuMuFactor(theGMM1->GetNewDoubleValue(aS));
  if(aComm==thePMM1)     theB->PositronToMuMuFactor(thePMM1->GetNewDoubleValue(aS));
  if(aComm==thePH1)      theB->PositronToHadronsFactor(thePH1->GetNewDoubleValue(aS));
}
