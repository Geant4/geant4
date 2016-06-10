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
// $Id: G4tgrMessenger.cc 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgrMessenger

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"

G4ThreadLocal G4int G4tgrMessenger::theVerboseLevel = 0;


// --------------------------------------------------------------------
G4tgrMessenger::G4tgrMessenger()
{
  tgDirectory = new G4UIdirectory("/geometry/textInput/");
  tgDirectory->SetGuidance("Geometry from text file control commands.");
  verboseCmd = new G4UIcmdWithAnInteger("/geometry/textInput/verbose",this);
  verboseCmd->SetGuidance("Set Verbose level of geometry text input category.");
  verboseCmd->SetGuidance(" 0 : Silent");
  verboseCmd->SetGuidance(" 1 : info verbosity");
  verboseCmd->SetGuidance(" 2 : debug verbosity");
  verboseCmd->SetParameterName("level",false);
  verboseCmd->SetRange("level>=0");
}


// --------------------------------------------------------------------
G4tgrMessenger::~G4tgrMessenger()
{
  delete verboseCmd;
  delete tgDirectory;
}


// --------------------------------------------------------------------
G4int G4tgrMessenger::GetVerboseLevel()
{
  return theVerboseLevel;
}


// --------------------------------------------------------------------
void G4tgrMessenger::SetVerboseLevel( G4int verb )
{
  theVerboseLevel = verb;
}


// --------------------------------------------------------------------
void G4tgrMessenger::SetNewValue(G4UIcommand * command, G4String newValues)
{
  if( command == verboseCmd )
  {
    G4tgrMessenger::SetVerboseLevel(verboseCmd->GetNewIntValue(newValues)); 
  }
}


// --------------------------------------------------------------------
G4String G4tgrMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command == verboseCmd )
  { 
    cv = verboseCmd->ConvertToString(G4tgrMessenger::GetVerboseLevel()); 
  }
  return cv;
}
