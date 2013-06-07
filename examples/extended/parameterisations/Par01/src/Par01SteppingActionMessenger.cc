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
// $Id$
//

#include "Par01SteppingActionMessenger.hh"

#include <sstream>

#include "Par01SteppingAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

Par01SteppingActionMessenger::Par01SteppingActionMessenger(Par01SteppingAction* SA)
:fSteppingAction(SA)
{
  fStepDirectory = new G4UIdirectory("/step/");
  fStepDirectory->SetGuidance("Step draw control command.");

  fDrawStepCmd = new G4UIcmdWithABool("/step/draw",this);
  fDrawStepCmd->SetGuidance("Draw each step on the fly. (default = true)");
  fDrawStepCmd->SetParameterName("draw", true);
  fDrawStepCmd->SetDefaultValue(true);
}

void Par01SteppingActionMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if( command->GetCommandName() == "draw" )
  {
    G4int vl;
    const char* t = newValues;
    std::istringstream is(t);
    is >> vl;
    fSteppingAction->SetDrawFlag(vl!=0);
  }
}

