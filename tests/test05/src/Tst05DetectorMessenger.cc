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
// $Id: Tst05DetectorMessenger.cc,v 1.4 2001-07-11 10:09:39 gunter Exp $
// ------------------------------------------------------------

#include "Tst05DetectorConstruction.hh"
#include "Tst05DetectorMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

Tst05DetectorMessenger::Tst05DetectorMessenger(Tst05DetectorConstruction * myDC):myDetector(myDC)
{
  G4String defParam;

  mydetDir = new G4UIdirectory("/mydet/");
  mydetDir->SetGuidance("Detector setup commands.");

  selDetCmd = new G4UIcmdWithAString("/mydet/SelectDetector",this);
  selDetCmd->SetGuidance("Select Detector Setup.");
  selDetCmd->SetGuidance("  Choice : rod_place_asm / CMSTracker ");
  selDetCmd->SetParameterName("choice",true);
  selDetCmd->SetDefaultValue("rod_place_asm");
  selDetCmd->SetCandidates("rod_place_asm CMSTracker");
  selDetCmd->AvailableForStates(PreInit,Idle);

  switchCmd = new G4UIcmdWithAString("/mydet/SwitchDetector",this);
  switchCmd->SetGuidance("Assign the selected geometry to G4RunManager.");
  switchCmd->SetGuidance("In cese detector name is associated to this command,");
  switchCmd->SetGuidance("\"/mydet/SelectDetector\" will be invoked and then switched.");
  switchCmd->SetParameterName("choice",true);
  switchCmd->SetDefaultValue(" ");
  switchCmd->SetCandidates("rod_place_asm CMSTracker \" \"");
  switchCmd->AvailableForStates(PreInit,Idle);

  myDetector->SelectDetector(defParam="rod_place_asm");
}

void Tst05DetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command == selDetCmd )
  {
    myDetector->SelectDetector(newValues);
  }
  if( command == switchCmd )
  {
    if(newValues=="rod_place_asm" || newValues=="CMSTracker")
    { myDetector->SelectDetector(newValues); }
  }
  
  return;
}

