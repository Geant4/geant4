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

 selDetCmd->SetCandidates("Sphere Orb Box Cone manyCons Tube Hype Torus Para Trd b1Ib2 b1Ub2 b1Sb2 b1Ub1 b1Ib1 b1Sb1 TwistedBox TwistedTrd TwistedTrap TwistedTrap2");
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









