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

#include "SCDetectorMessenger.hh"

#include "SCDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

#include "G4ios.hh"

SCDetectorMessenger::SCDetectorMessenger(SCDetectorConstruction * myDC)
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

 selDetCmd->SetCandidates("Sphere Orb Box Cone manyCons Tube Hype Torus Para Trd b1Ib2 b1Ub2 b1Sb2 b1Ub1 b1Ib1 b1Sb1 TwistedTubs TwistedBox TwistedTrd TwistedTrap TwistedTrap2 TwistedTrap3 Ellipsoid");
  selDetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    myDetector->SelectDetector(defParam="Sphere");
}

void SCDetectorMessenger::SetNewValue(G4UIcommand * command,
                                         G4String newValues)
{
  if( command == selDetCmd )
  {
    myDetector->SelectDetector(newValues);
    myDetector->SwitchDetector();
  }
  return;
}









