// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelPrimaryGeneratorActionMessenger.cc       *
// * -------                                                            *
// *                                                                    *
// * Version:           0.5                                             *
// * Date:              08/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R. Nartallo
// - First implementation of PrimaryGeneratorMessenger
// - Based on Chandra and XMM models by S Magni and F Lei
// 
// 08.11.2000 R. Nartallo
// - Removed line "delete ErndmCmd" from destructor
//
// **********************************************************************

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "XrayTelPrimaryGeneratorAction.hh"
#include "XrayTelPrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelPrimaryGeneratorMessenger::XrayTelPrimaryGeneratorMessenger(XrayTelPrimaryGeneratorAction* XrayTelGun)
  :XrayTelAction(XrayTelGun)
{ 
  RndmCmd = new G4UIcmdWithAString("/gun/random",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : aperture, point");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("aperture");
  RndmCmd->SetCandidates("aperture point");
  RndmCmd->AvailableForStates(PreInit,Idle);

  SetRmin  = new G4UIcmdWithADoubleAndUnit("/gun/SetRmin",this);
  SetRmin->SetGuidance("Set the minimum radius for random position generation");
  SetRmin->SetParameterName("Rmin",true);
  SetRmin->SetDefaultUnit("cm"); 
  SetRmin->SetDefaultValue(0.0*cm);

  SetRmax = new G4UIcmdWithADoubleAndUnit("/gun/SetRmax",this);
  SetRmax->SetGuidance("Set the maximum radius for random position generation");
  SetRmax->SetParameterName("Rmax",true);
  SetRmax->SetDefaultUnit("cm"); 
  SetRmax->SetDefaultValue(100.0*cm);

  SetTmax = new G4UIcmdWithADoubleAndUnit("/gun/SetTmax",this);
  SetTmax->SetGuidance("Set the maximum aperture angle for direction generation");
  SetTmax->SetParameterName("Tmax",true);
  SetTmax->SetDefaultUnit("degree"); 
  SetTmax->SetDefaultValue(1*degree);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelPrimaryGeneratorMessenger::~XrayTelPrimaryGeneratorMessenger()
{
  delete RndmCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayTelPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == RndmCmd ){ 
    XrayTelAction->SetRndmFlag(newValue);
  } else if( command == ErndmCmd ){ 
    XrayTelAction->SetErndmFlag(newValue);  
  } else if  ( command == SetRmin ){
    XrayTelAction->SetRmin( SetRmin->GetNewDoubleValue ( newValue )  );
  } else if  ( command == SetRmax ){
    XrayTelAction->SetRmax( SetRmax->GetNewDoubleValue ( newValue )  );
  } else if ( command == SetTmax ){
    XrayTelAction->SetTmax( SetTmax->GetNewDoubleValue ( newValue )  );
  }
}








