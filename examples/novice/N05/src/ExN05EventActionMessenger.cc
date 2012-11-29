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

#include "ExN05EventActionMessenger.hh"

#include <sstream>

#include "ExN05EventAction.hh"
#include "G4UIcmdWithABool.hh"

ExN05EventActionMessenger::ExN05EventActionMessenger(ExN05EventAction* EA)
:EventAction(EA)
{
  // uses "/event" directory already
  // defined in G4EvManMessenger
  drawEventCmd = new G4UIcmdWithABool("/event/Draw",this);
  drawEventCmd->SetGuidance("Set drawFlag to Draw an event. (default = false)");
  drawEventCmd->SetParameterName("Draw", false);
  drawEventCmd->SetDefaultValue(false);
}

void ExN05EventActionMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if( command->GetCommandName() == "Draw" )
  {
    G4int vl;
    const char* t = newValues;
    std::istringstream is(t);
    is >> vl;
    EventAction->SetDrawFlag(vl!=0);
  }
}

