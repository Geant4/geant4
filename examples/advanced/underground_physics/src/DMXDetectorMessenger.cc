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
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// SteppingActionMessenger program
// --------------------------------------------------------------

#include "DMXDetectorMessenger.hh"

#include "DMXDetectorConstruction.hh"

#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"


DMXDetectorMessenger::DMXDetectorMessenger
   (DMXDetectorConstruction* DC):detectorConstruction(DC) {

  RoomEKineCutCmd = new G4UIcmdWithADoubleAndUnit("/dmx/RoomMinEnergyCut",this);
  RoomEKineCutCmd->SetGuidance("Minimum Charged particle cut in ROOM");
  RoomEKineCutCmd->SetParameterName("ECut",false,false);
  RoomEKineCutCmd->SetRange("ECut>=250.0*eV");
  RoomEKineCutCmd->SetDefaultUnit("eV");
  RoomEKineCutCmd->SetUnitCategory("Energy");
  RoomEKineCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  EKineCutCmd = new G4UIcmdWithADoubleAndUnit("/dmx/MinEnergyCut",this);
  EKineCutCmd->SetGuidance("Minimum Charged particle cut inside detector");
  EKineCutCmd->SetParameterName("ECut",false,false);
  EKineCutCmd->SetRange("ECut>=250.0*eV");
  EKineCutCmd->SetDefaultUnit("eV");
  EKineCutCmd->SetUnitCategory("Energy");
  EKineCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RoomTimeCutCmd = new G4UIcmdWithADoubleAndUnit("/dmx/RoomTimeCut",this);
  RoomTimeCutCmd->SetGuidance("Set Time Cut (for neutrons) inside ROOM");
  RoomTimeCutCmd->SetParameterName("RTCut",false,false);
  RoomTimeCutCmd->SetRange("RTCut>0.");
  RoomTimeCutCmd->SetDefaultUnit("ns");
  RoomTimeCutCmd->SetUnitCategory("Time");
  RoomTimeCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TimeCutCmd = new G4UIcmdWithADoubleAndUnit("/dmx/TimeCut",this);
  TimeCutCmd->SetGuidance("Set Time Cut (for neutrons) inside detector");
  TimeCutCmd->SetParameterName("TCut",false,false);
  TimeCutCmd->SetRange("TCut>0.");
  TimeCutCmd->SetDefaultUnit("ns");
  TimeCutCmd->SetUnitCategory("Time");
  TimeCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //UpdateCmd = new G4UIcmdWithoutParameter("/dmx/update",this);
  //UpdateCmd->SetGuidance("Update calorimeter geometry.");
  //UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  //UpdateCmd->SetGuidance("if you changed timecut value(s).");
  //UpdateCmd->AvailableForStates(G4State_Idle);

}


//ooooooooooooooooooooooooooooooooooooooooo
DMXDetectorMessenger::~DMXDetectorMessenger() {

  delete RoomTimeCutCmd;
  delete TimeCutCmd;
  //  delete UpdateCmd;

}


//ooooooooooooooooooooooooooooooooooooooooo
void DMXDetectorMessenger::SetNewValue(G4UIcommand* command, 
  G4String newValue) {

  if(command == EKineCutCmd)
   detectorConstruction->
     SetEnergyCut(EKineCutCmd->GetNewDoubleValue(newValue));

  if(command == RoomEKineCutCmd)
   detectorConstruction->
     SetEnergyCut(RoomEKineCutCmd->GetNewDoubleValue(newValue));

  if(command == TimeCutCmd)
   detectorConstruction->
     SetTimeCut(TimeCutCmd->GetNewDoubleValue(newValue));

  if(command == RoomTimeCutCmd)
    detectorConstruction->
      SetRoomTimeCut(RoomTimeCutCmd->GetNewDoubleValue(newValue));

  //  if( command == UpdateCmd )
  //  { detectorConstruction->UpdateGeometry(); }

}



