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
  aDir1 = new G4UIdirectory("/physics_lists/", false);
  aDir1->SetGuidance("commands for physics list configuration.");

  // general stuff.
  aDir2 = new G4UIdirectory("/physics_lists/em/", false);
  aDir2->SetGuidance("Extra EM processes configuration.");

  // command for synchrotron radiation.
  theSynch = new G4UIcmdWithABool("/physics_lists/em/SyncRadiation",this);
  theSynch->SetGuidance("Switching on/off synchrotron radiation.");
  theSynch->AvailableForStates(G4State_PreInit);
  theSynch->SetToBeBroadcasted(false);

  // command for synchrotron radiation.
  theSynchAll = new G4UIcmdWithABool("/physics_lists/em/SyncRadiationAll",this);
  theSynchAll->SetGuidance("Switching on/off synchrotron radiation for all charged.");
  theSynchAll->AvailableForStates(G4State_PreInit);
  theSynchAll->SetToBeBroadcasted(false);

  // command for gamma nuclear physics.
  theGN = new G4UIcmdWithABool("/physics_lists/em/GammaNuclear",this);
  theGN->SetGuidance("Switching on gamma nuclear physics.");
  theGN->AvailableForStates(G4State_PreInit);
  theGN->SetToBeBroadcasted(false);

  // command for gamma nuclear physics.
  theXS = new G4UIcmdWithABool("/physics_lists/em/UseGammaNuclearXS",this);
  theXS->SetGuidance("Use XS gamma nuclear cross section.");
  theXS->AvailableForStates(G4State_PreInit);
  theXS->SetToBeBroadcasted(false);

  // command for lend gamma nuclear physics.
  theGLENDN = new G4UIcmdWithABool("/physics_lists/em/LENDGammaNuclear",this);
  theGLENDN->SetGuidance("Switching on LEND gamma nuclear physics.");
  theGLENDN->AvailableForStates(G4State_PreInit);
  theGLENDN->SetToBeBroadcasted(false);

  theEN = new G4UIcmdWithABool("/physics_lists/em/ElectroNuclear",this);
  theEN->SetGuidance("Switching on e+- nuclear physics.");
  theEN->AvailableForStates(G4State_PreInit);
  theEN->SetToBeBroadcasted(false);

  // command for muon nuclear physics.
  theMUN = new G4UIcmdWithABool("/physics_lists/em/MuonNuclear",this);
  theMUN->SetGuidance("Switching on muon nuclear physics.");
  theMUN->AvailableForStates(G4State_PreInit);
  theMUN->SetToBeBroadcasted(false);

  theGMM = new G4UIcmdWithABool("/physics_lists/em/GammaToMuons",this);
  theGMM->SetGuidance("Switching on gamma conversion to muon pair.");
  theGMM->AvailableForStates(G4State_PreInit);
  theGMM->SetToBeBroadcasted(false);

  theMMM = new G4UIcmdWithABool("/physics_lists/em/MuonToMuons",this);
  theMMM->SetGuidance("Switching on muon pair production by muons.");
  theMMM->AvailableForStates(G4State_PreInit);
  theMMM->SetToBeBroadcasted(false);

  thePMM = new G4UIcmdWithABool("/physics_lists/em/PositronToMuons",this);
  thePMM->SetGuidance("Switching on positron conversion to muon pair.");
  thePMM->AvailableForStates(G4State_PreInit);
  thePMM->SetToBeBroadcasted(false);

  thePH = new G4UIcmdWithABool("/physics_lists/em/PositronToHadrons",this);
  thePH->SetGuidance("Switching on positron conversion to hadrons.");
  thePH->AvailableForStates(G4State_PreInit);
  thePH->SetToBeBroadcasted(false);

  theGMM1 = new G4UIcmdWithADouble("/physics_lists/em/GammaToMuonsFactor",this);
  theGMM1->SetGuidance("Factor for gamma conversion to muon pair.");
  theGMM1->AvailableForStates(G4State_PreInit);
  theGMM1->SetToBeBroadcasted(false);

  thePMM1 = new G4UIcmdWithADouble("/physics_lists/em/PositronToMuonsFactor",this);
  thePMM1->SetGuidance("Factor for positron conversion to muon pair.");
  thePMM1->AvailableForStates(G4State_PreInit);
  thePMM1->SetToBeBroadcasted(false);

  thePH1 = new G4UIcmdWithADouble("/physics_lists/em/PositronToHadronsFactor",this);
  thePH1->SetGuidance("Factor for positron conversion to hadrons.");
  thePH1->AvailableForStates(G4State_PreInit);
  thePH1->SetToBeBroadcasted(false);

  theGNlowe = new G4UIcmdWithADoubleAndUnit("/physics_lists/em/GammaNuclearLEModelLimit",this);
  theGNlowe->SetGuidance("Upper energy limit for low-energy gamma-nuclear model");
  theGNlowe->SetParameterName("emin",true);
  theGNlowe->SetUnitCategory("Energy");
  theGNlowe->AvailableForStates(G4State_PreInit);
  theGNlowe->SetToBeBroadcasted(false);
}

G4EmMessenger::~G4EmMessenger()
{
  delete theSynch;
  delete theSynchAll;
  delete theGN;
  delete theGLENDN;
  delete theEN;
  delete theMUN;
  delete theGMM;
  delete theMMM;
  delete thePMM;
  delete thePH;

  delete theGMM1;
  delete thePMM1;
  delete thePH1;
  delete theGNlowe;
  delete theXS;

  delete aDir1;
  delete aDir2;
}

void G4EmMessenger::SetNewValue(G4UIcommand* aComm, G4String aS)
{
  if (aComm==theSynch)
    theB->Synch(theSynch->GetNewBoolValue(aS));
  else if (aComm==theSynchAll)
    theB->SynchAll(theSynchAll->GetNewBoolValue(aS));
  else if (aComm==theGN)
    theB->GammaNuclear(theGN->GetNewBoolValue(aS));
  else if (aComm==theGLENDN)
    theB->LENDGammaNuclear(theGLENDN->GetNewBoolValue(aS));
  else if (aComm==theEN)
    theB->ElectroNuclear(theEN->GetNewBoolValue(aS));
  else if (aComm==theMUN)
    theB->MuonNuclear(theMUN->GetNewBoolValue(aS));
  else if (aComm==theGMM)
    theB->GammaToMuMu(theGMM->GetNewBoolValue(aS));
  else if (aComm==theMMM)
    theB->MuonToMuMu(theMMM->GetNewBoolValue(aS));
  else if (aComm==thePMM)
    theB->PositronToMuMu(thePMM->GetNewBoolValue(aS));
  else if (aComm==thePH)
    theB->PositronToHadrons(thePH->GetNewBoolValue(aS));
  else if (aComm==theXS)
    theB->SetUseGammaNuclearXS(theXS->GetNewBoolValue(aS));
  else if (aComm==theGMM1)
    theB->GammaToMuMuFactor(theGMM1->GetNewDoubleValue(aS));
  else if (aComm==thePMM1)
    theB->PositronToMuMuFactor(thePMM1->GetNewDoubleValue(aS));
  else if (aComm==thePH1)
    theB->PositronToHadronsFactor(thePH1->GetNewDoubleValue(aS));
  else if (aComm==theGNlowe)
    theB->GammaNuclearLEModelLimit(theGNlowe->GetNewDoubleValue(aS));
}
