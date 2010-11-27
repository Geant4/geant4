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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#include "IAEADetectorMessenger.hh"
#include "IAEADetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

IAEADetectorMessenger::IAEADetectorMessenger(IAEADetectorConstruction* detector)
  :IAEADetector(detector)
{

  setIAEAWaterPhantomThicknessCmd = new G4UIcmdWithADoubleAndUnit("/iaea/waterPhantomThickness",this);
  setIAEAWaterPhantomThicknessCmd -> SetGuidance("Set size of water phantom");
  setIAEAWaterPhantomThicknessCmd -> SetParameterName("Size",false);
  setIAEAWaterPhantomThicknessCmd -> SetDefaultUnit("cm");  
  setIAEAWaterPhantomThicknessCmd -> SetUnitCandidates("mm cm m");  
  setIAEAWaterPhantomThicknessCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);
}

IAEADetectorMessenger::~IAEADetectorMessenger()
{ 
  delete setIAEAWaterPhantomThicknessCmd;
}

void IAEADetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == setIAEAWaterPhantomThicknessCmd ){
    IAEADetector->setWaterThickness(setIAEAWaterPhantomThicknessCmd->GetNewDoubleValue(newValue));
  }

}

