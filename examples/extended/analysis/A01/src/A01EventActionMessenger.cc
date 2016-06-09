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
/// \file analysis/A01/src/A01EventActionMessenger.cc
/// \brief Implementation of the A01EventActionMessenger class
//
// $Id$
// --------------------------------------------------------------
//

#include "A01EventActionMessenger.hh"
#include "A01EventAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"

A01EventActionMessenger::A01EventActionMessenger(A01EventAction * mpga)
:fTarget (mpga)
{
  fVerboseCmd = new G4UIcmdWithAnInteger("/mydet/verbose",this);
  fVerboseCmd->SetGuidance("Verbose level for each event.");
  fVerboseCmd->SetGuidance(" Event summary will be displayed for every 'level' events.");
  fVerboseCmd->SetParameterName("level",true);
  fVerboseCmd->SetRange("level>=0");
  fVerboseCmd->SetDefaultValue(1);
}

A01EventActionMessenger::~A01EventActionMessenger()
{
  delete fVerboseCmd;
}

void A01EventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==fVerboseCmd )
  { fTarget ->SetVerbose(fVerboseCmd->GetNewIntValue(newValue)); }
}

G4String A01EventActionMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==fVerboseCmd )
  { cv = fVerboseCmd->ConvertToString(fTarget ->GetVerbose()); }

  return cv;
}
