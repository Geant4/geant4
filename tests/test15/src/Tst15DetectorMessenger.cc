// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst15DetectorMessenger.cc,v 1.1 1999-11-18 14:48:14 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst15DetectorMessenger.hh"

#include "Tst15DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

#include "G4ios.hh"

Tst15DetectorMessenger::Tst15DetectorMessenger(Tst15DetectorConstruction * myDC)
:myDetector(myDC)
{
  G4String defParam;

  mydetDir = new G4UIdirectory("/mydet/");
  mydetDir->SetGuidance("Detector setup commands.");

  selMatCmd = new G4UIcmdWithAString("/mydet/SelectMaterial",this);
  selMatCmd->SetGuidance("Select Material of the SimpleBox.");
  selMatCmd->SetGuidance("  Choice : Air, Al, Pb (default)");
  selMatCmd->SetParameterName("choice",true);
  selMatCmd->SetDefaultValue("Pb");
  selMatCmd->SetCandidates("Air Al Pb");
  selMatCmd->AvailableForStates(PreInit,Idle);

  myDetector->SelectMaterial(defParam="Pb");
}

void Tst15DetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command == selMatCmd )
  {
    myDetector->SelectMaterial(newValues);
  }
  return;
}

