// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05DetectorMessenger.cc,v 1.1 1999-01-07 16:06:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "ExN05DetectorMessenger.hh"
#include "ExN05DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"
#include <iomanip.h>                

#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

ExN05DetectorMessenger::ExN05DetectorMessenger(ExN05DetectorConstruction * myDet) : myDetector(myDet)
{ 

  mydetDir = new G4UIdirectory("/exN05/detector/");
  mydetDir->SetGuidance("ExN05 detector control.");
  
  SwitchCmd = new G4UIcmdWithAString("/exN05/detector/switchUserLimits",this);
  SwitchCmd->SetGuidance("Switch On/Off UserLimits for crystals");
  SwitchCmd->SetGuidance(" Choice : on / off ");
  SwitchCmd->SetParameterName("choice",true);
  SwitchCmd->SetDefaultValue("on");
  SwitchCmd->SetCandidates("on ON off OFF");
  SwitchCmd->AvailableForStates(PreInit,Idle);
  
  TmaxCmd = new G4UIcmdWithADoubleAndUnit("/exN05/detector/maxTime",this);
  TmaxCmd->SetGuidance("Set maximum time in crystal");
  TmaxCmd->SetParameterName("Tmax",false,false);
  TmaxCmd->SetDefaultUnit("ns");
  TmaxCmd->SetUnitCategory("Time");
  TmaxCmd->AvailableForStates(PreInit,Idle);

  EminCmd = new G4UIcmdWithADoubleAndUnit("/exN05/detector/minEkine",this);
  EminCmd->SetGuidance("Set minimum kinetic energy in crystal");
  EminCmd->SetParameterName("Emin",false,false);
  EminCmd->SetDefaultUnit("MeV");
  EminCmd->SetUnitCategory("Energy");
  EminCmd->AvailableForStates(PreInit,Idle);

  RminCmd = new G4UIcmdWithADoubleAndUnit("/exN05/detector/minRange",this);
  RminCmd->SetGuidance("Set minimum range in crystal");
  RminCmd->SetParameterName("Rmin",false,false);
  RminCmd->SetDefaultUnit("mm");
  RminCmd->SetUnitCategory("Length");
  RminCmd->AvailableForStates(PreInit,Idle);
}

ExN05DetectorMessenger::~ExN05DetectorMessenger()
{
  delete SwitchCmd;
  delete TmaxCmd;
  delete RminCmd;
  delete EminCmd;
}

void ExN05DetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  if( command == SwitchCmd ) { 
    if ( (newValues == "off") || (newValues == "OFF") ){
      myDetector->UseUserLimits(false);
    } else {
      myDetector->UseUserLimits(true);
    } 

  } else if( command == TmaxCmd ) {
    myDetector->SetMaxTimeInCrystal(TmaxCmd->GetNewDoubleValue(newValues));

  } else if( command == EminCmd ) {
    myDetector->SetMinEkineInCrystal(EminCmd->GetNewDoubleValue(newValues));

  } else if( command == RminCmd ) {
    myDetector->SetMinRangeInCrystal(RminCmd->GetNewDoubleValue(newValues));
  }
}

G4String ExN05DetectorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String returnValue('\0');
  char line[255];
  ostrstream os(line,255);
  
  if( command == SwitchCmd ) { 
    if ( myDetector->IsUseUserLimits() )returnValue = "on";
    else returnValue = "off";

  } else if( command == TmaxCmd ) {
    os <<  myDetector->GetMaxTimeInCrystal()/ns << " ns" << '\0';
    returnValue = G4String(line);

  } else if( command == EminCmd ) {
    os <<  myDetector->GetMinEkineInCrystal()/MeV << " MeV" << '\0';
    returnValue = G4String(line);

  } else if( command == RminCmd ) {
    os <<  myDetector->GetMinRangeInCrystal()/mm << " mm" << '\0';
    returnValue = G4String(line);
  }
  return returnValue;
}
