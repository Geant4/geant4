//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#include "Sc01DetectorMessenger.hh"

#include "Sc01DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

#include "G4ios.hh"

Sc01DetectorMessenger::Sc01DetectorMessenger(Sc01DetectorConstruction * myDC)
  : myDetector(myDC)
{
  G4String defParam;

  mydetDir = new G4UIdirectory("/mydet/");
  mydetDir->SetGuidance("Detector setup commands.");

  selDetCmd = new G4UIcmdWithAString("/mydet/SelectDetector",this);
  selDetCmd->SetGuidance("Select Detector Setup.");
  selDetCmd->SetGuidance("  Choice : Detector type ");
  selDetCmd->SetParameterName("choice",true);
  selDetCmd->SetDefaultValue("Sphere");
  selDetCmd->SetCandidates("Sphere Orb Box Cone Tube Hype Torus Para Trd");
  selDetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  myDetector->SelectDetector(defParam="Sphere");
}

void Sc01DetectorMessenger::SetNewValue(G4UIcommand * command,
                                         G4String newValues)
{
  if( command == selDetCmd )
  {
    myDetector->SelectDetector(newValues);
    myDetector->SwitchDetector();
  }
  return;
}









