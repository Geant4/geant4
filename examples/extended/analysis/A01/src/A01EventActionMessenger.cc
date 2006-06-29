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
// $Id: A01EventActionMessenger.cc,v 1.6 2006-06-29 16:32:39 gunter Exp $
// --------------------------------------------------------------
//

#include "A01EventActionMessenger.hh"
#include "A01EventAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"

A01EventActionMessenger::A01EventActionMessenger(A01EventAction * mpga)
:target(mpga)
{
  verboseCmd = new G4UIcmdWithAnInteger("/mydet/verbose",this);
  verboseCmd->SetGuidance("Verbose level for each event.");
  verboseCmd->SetGuidance(" Event summary will be displayed for every 'level' events.");
  verboseCmd->SetParameterName("level",true);
  verboseCmd->SetRange("level>=0");
  verboseCmd->SetDefaultValue(1);
}

A01EventActionMessenger::~A01EventActionMessenger()
{
  delete verboseCmd;
}

void A01EventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==verboseCmd )
  { target->SetVerbose(verboseCmd->GetNewIntValue(newValue)); }
}

G4String A01EventActionMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==verboseCmd )
  { cv = verboseCmd->ConvertToString(target->GetVerbose()); }

  return cv;
}
