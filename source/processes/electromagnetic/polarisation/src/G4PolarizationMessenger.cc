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
// $Id: G4PolarizationMessenger.cc 81375 2014-05-27 13:08:40Z gcosmo $
//
//
// GEANT4 Class file
//
// File name:     G4PolarizationManager
//
// Author:        Andreas Schaelicke
//
// Creation date: 01.05.2005
//
// Modifications:
//
// Class Description:
//
// Provides access to general polarization information and to 
// polarization for logical volumes through macro files.

#include "G4PolarizationMessenger.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizationHelper.hh"

#include "G4UIdirectory.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"


G4PolarizationMessenger::G4PolarizationMessenger(G4PolarizationManager * polMgr)
  : polarizationManager(polMgr)
{
  polarizationDirectory = new G4UIdirectory("/polarization/");
  polarizationDirectory->SetGuidance("polarization control commands.");

  managerDirectory = new G4UIdirectory("/polarization/manager/");
  managerDirectory->SetGuidance("general polarization information.");

  verboseCmd = new G4UIcmdWithAnInteger("/polarization/manager/verbose",this);
  verboseCmd->SetGuidance("Set the Verbose level of G4PolarizationManager.");
  verboseCmd->SetGuidance(" 0 : Silent (default)");
  verboseCmd->SetGuidance(" 1 : Verbose");
  verboseCmd->SetParameterName("level",true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->SetRange("level >=0 && level <=1");

  optActivateCmd = new G4UIcmdWithABool("/polarization/manager/activate",this);
  optActivateCmd->SetGuidance("activate/deactivate polarization treatment");
  optActivateCmd->SetParameterName("flag",true);
  optActivateCmd->SetDefaultValue(true);

  volumeDirectory = new G4UIdirectory("/polarization/volume/");
  volumeDirectory->SetGuidance("Status control commands of registered polarized logical volumes.");

  printVolumeListCmd = new G4UIcmdWithoutParameter("/polarization/volume/list",this);
  printVolumeListCmd->SetGuidance("print list of registered polarized logical volumes");
  printVolumeListCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);
  
  setPolarizationCmd = new G4UIcommand("/polarization/volume/set",this);
  setPolarizationCmd->SetGuidance("set or change polarization of a logical volume");
//   setPolarizationCmd->SetParameterName("polarization",true);
//   setPolarizationCmd->SetDefaultValue("worldVolume  0. 0. 0.");
  setPolarizationCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);

  G4UIparameter* param;
  param = new G4UIparameter("logicalVolumeName",'s',false);
  param->SetDefaultValue("worldVolume");
  setPolarizationCmd->SetParameter(param);
  param = new G4UIparameter("px",'d',true);
  param->SetDefaultValue("0.0");
  setPolarizationCmd->SetParameter(param);
  param = new G4UIparameter("py",'d',true);
  param->SetDefaultValue("0.0");
  setPolarizationCmd->SetParameter(param);
  param = new G4UIparameter("pz",'d',true);
  param->SetDefaultValue("0.0");
  setPolarizationCmd->SetParameter(param);

  testDirectory = new G4UIdirectory("/polarization/test/");
  testDirectory->SetGuidance("provides access to some internal test routines.");

  testPolarizationTransformationCmd = new G4UIcmdWithoutParameter("/polarization/test/polarizationTransformation",this);
  testPolarizationTransformationCmd->SetGuidance("checks definition of particle reference frame and corresponding translation routines"); 
  testPolarizationTransformationCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);

  testInteractionFrameCmd = new G4UIcmdWithoutParameter("/polarization/test/interactionFrame",this);
  testInteractionFrameCmd->SetGuidance("checks definition of interaction frame"); 
  testInteractionFrameCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);
}

G4PolarizationMessenger::~G4PolarizationMessenger()
{
  delete verboseCmd;
  delete testInteractionFrameCmd;
  delete testPolarizationTransformationCmd;
  delete testDirectory;
  delete setPolarizationCmd;
  delete printVolumeListCmd;
  delete volumeDirectory;
  delete optActivateCmd;
  delete managerDirectory;
  delete polarizationDirectory;
}

void G4PolarizationMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==verboseCmd ) { 
    polarizationManager->SetVerbose(verboseCmd->GetNewIntValue(newValue)); 
  }
  else if ( command==optActivateCmd ) { 
    polarizationManager->SetActivated(optActivateCmd->GetNewBoolValue(newValue));
  }
  else if ( command==printVolumeListCmd ) { 
    polarizationManager->ListVolumes();
  }
  else if ( command==setPolarizationCmd ) {
    G4Tokenizer next( newValue );
    G4String volumeName=next();
    G4double px=0.,py=0.,pz=0.;
    G4String dvalue=next();
    if (!dvalue.isNull()) {
      px=StoD(dvalue);
      dvalue=next();
      if (!dvalue.isNull()) {
	py=StoD(dvalue);
	dvalue=next();
	if (!dvalue.isNull()) pz=StoD(dvalue);
      }
    }
    G4ThreeVector pol(px,py,pz);
    polarizationManager->SetVolumePolarization(volumeName,pol);
  }
  else if ( command==testPolarizationTransformationCmd ) {
    G4PolarizationHelper::TestPolarizationTransformations();
  }
  else if (command==testInteractionFrameCmd ) {
    G4PolarizationHelper::TestInteractionFrame();
  }
}

G4String G4PolarizationMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==verboseCmd )
  { cv = verboseCmd->ConvertToString(polarizationManager->GetVerbose()); }

  return cv;
}
