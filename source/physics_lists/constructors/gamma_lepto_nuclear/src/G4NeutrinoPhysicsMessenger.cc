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
// ClassName:   G4NeutrinoPhysicsMessenger
//
// Author: 2023 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4NeutrinoPhysicsMessenger.hh"
#include "G4NeutrinoPhysics.hh"

G4NeutrinoPhysicsMessenger::G4NeutrinoPhysicsMessenger(G4NeutrinoPhysics* ab)
  : theB(ab)
{
  // general stuff.
  aDir = new G4UIdirectory("/physics_lists/nu/", false);
  aDir->SetGuidance("tailoring the neutrino processes.");

  theNu = new G4UIcmdWithABool("/physics_lists/nu/NeutrinoActivation",this);
  theNu->SetGuidance("Activation of neutrino-nucleus processes");
  theNu->AvailableForStates(G4State_PreInit);
  theNu->SetToBeBroadcasted(false);

  theNuETX = new G4UIcmdWithABool("/physics_lists/nu/NuETotXscActivation",this);
  theNuETX->SetGuidance("Activation of neutrino-electron processes");
  theNuETX->AvailableForStates(G4State_PreInit);
  theNuETX->SetToBeBroadcasted(false);

  theNuEleCcBF = new G4UIcmdWithADouble("/physics_lists/nu/NuEleCcBias",this);
  theNuEleCcBF->SetGuidance("Neutrino-electron charge current bias factor");
  theNuEleCcBF->AvailableForStates(G4State_PreInit);
  theNuEleCcBF->SetToBeBroadcasted(false);

  theNuEleNcBF = new G4UIcmdWithADouble("/physics_lists/nu/NuEleNcBias",this);
  theNuEleNcBF->SetGuidance("Neutrino-electron neutral current bias factor");
  theNuEleNcBF->AvailableForStates(G4State_PreInit);
  theNuEleNcBF->SetToBeBroadcasted(false);

  theNuNucleusBF = new G4UIcmdWithADouble("/physics_lists/nu/NuNucleusBias",this);
  theNuNucleusBF->SetGuidance("Neutrino-nucleus cross section bias factor");
  theNuNucleusBF->AvailableForStates(G4State_PreInit);
  theNuNucleusBF->SetToBeBroadcasted(false);

  theNuOscDistanceBF = new G4UIcmdWithADouble("/physics_lists/nu/NuOscDistanceBias",this);
  theNuOscDistanceBF->SetGuidance("Neutrino-oscillation distance bias factor");
  theNuOscDistanceBF->AvailableForStates(G4State_PreInit);
  theNuOscDistanceBF->SetToBeBroadcasted(false);

  theNuDN = new G4UIcmdWithAString("/physics_lists/nu/NuDetectorName",this);  
  theNuDN->SetGuidance("Set neutrino detector name");
  theNuDN->AvailableForStates(G4State_PreInit);
  theNuDN->SetToBeBroadcasted(false);

  theNuODN = new G4UIcmdWithAString("/physics_lists/nu/NuOscDistanceName",this);  
  theNuODN->SetGuidance("Set neutrino oscillation distance region name");
  theNuODN->AvailableForStates(G4State_PreInit);
  theNuODN->SetToBeBroadcasted(false);
}

G4NeutrinoPhysicsMessenger::~G4NeutrinoPhysicsMessenger()
{
  delete theNu;
  delete theNuETX;

  delete theNuEleCcBF;
  delete theNuEleNcBF;
  delete theNuNucleusBF;
  delete theNuOscDistanceBF;

  delete theNuDN;
  delete theNuODN;

  delete aDir;
}

void G4NeutrinoPhysicsMessenger::SetNewValue(G4UIcommand* aComm, G4String aS)
{
  if (aComm==theNuETX)
    theB->NuETotXscActivated(theNuETX->GetNewBoolValue(aS));
  else if (aComm==theNuEleCcBF)
    theB->SetNuEleCcBias(theNuEleCcBF->GetNewDoubleValue(aS));
  else if (aComm==theNuEleNcBF)
    theB->SetNuEleNcBias(theNuEleNcBF->GetNewDoubleValue(aS));
  else if (aComm==theNuNucleusBF)
    theB->SetNuNucleusBias(theNuNucleusBF->GetNewDoubleValue(aS));
  else if (aComm==theNuOscDistanceBF)
    theB->SetNuOscDistanceBias(theNuOscDistanceBF->GetNewDoubleValue(aS));
  else if(aComm==theNuDN)
    theB->SetNuDetectorName(aS);
  else if(aComm==theNuODN)
    theB->SetNuOscDistanceName(aS);
}
