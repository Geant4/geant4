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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4RadioactiveDecayBaseMessenger.cc                                //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   29 August 2017                                                    //
//  Description: messenger class for the non-biased version of                //
//               G4RadioactiveDecay.  Based on the code of F. Lei and         //
//               P.R. Truscott.                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4RadioactiveDecayBaseMessenger.hh"
#include "G4NuclearLevelData.hh"
#include <sstream>


G4RadioactiveDecayBaseMessenger::G4RadioactiveDecayBaseMessenger 
(G4RadioactiveDecayBase* theRadioactiveDecayContainer1)
:theRadioactiveDecayContainer(theRadioactiveDecayContainer1)
{
  // main directory for control of the RDM
  grdmDirectory = new G4UIdirectory("/grdm/");
  grdmDirectory->SetGuidance("Controls for the Radioactive Decay Module.");

  // Command to define the limits on nucleus the RDM will treat.
  nucleuslimitsCmd = new
    G4UIcmdWithNucleusLimits("/grdm/nucleusLimits",this);
  nucleuslimitsCmd->SetGuidance 
    ("Set the atomic weight and number limits for the RDM.");
  nucleuslimitsCmd->SetParameterName("aMin","aMax","zMin","zMax",true);

  // Command controlling whether beta decay will be treated faithfully or 
  // in fast mode (obsolete) 
  fbetaCmd = new G4UIcmdWithABool ("/grdm/fBeta",this);
//  fbetaCmd->SetGuidance("false: use 3-body decay, true: use histogram method");
  fbetaCmd->SetGuidance("Fast (approximate) beta decay no longer used.");
  fbetaCmd->SetGuidance("Kept for backward compatibility.");
  fbetaCmd->SetParameterName("fBeta",true);
  fbetaCmd->SetDefaultValue(false);

  // Select a logical volume for RDM
  avolumeCmd = new
    G4UIcmdWithAString("/grdm/selectVolume",this);
  avolumeCmd->SetGuidance 
    ("Suppply a logical volumes name to add it to the RDM apply list");
  avolumeCmd->SetParameterName("aVolume",false);

  // De-select a logical volume for RDM
  deavolumeCmd = new
    G4UIcmdWithAString("/grdm/deselectVolume",this);
  deavolumeCmd->SetGuidance 
    ("Suppply a logical volumes name to remove it from the RDM apply list");
  deavolumeCmd->SetParameterName("aVolume",false);

  // Select all logical volumes for RDM
  allvolumesCmd = new
    G4UIcmdWithoutParameter("/grdm/allVolumes",this);
  allvolumesCmd->SetGuidance 
    (" apply RDM to all logical volumes. No parameter required.");
  //  allvolumeCmd->SetParameterName("AddAVolume",true);

  // De-select all logical volumes for RDM
  deallvolumesCmd = new
    G4UIcmdWithoutParameter("/grdm/noVolumes",this);
  deallvolumesCmd->SetGuidance 
    (" RDM is not applied to any logical volumes");
  //  deallvolumesCmd->SetParameterName("RemoveAVolume",true);z

  // Command to invoke internal conversion or not
  icmCmd = new G4UIcmdWithABool ("/grdm/applyICM",this);
  icmCmd->SetGuidance("Command not active; kept for backward compatibility.");
  icmCmd->SetGuidance("Internal conversion is always turned on.");
  icmCmd->SetParameterName("applyICM",true);
  icmCmd->SetDefaultValue(true);

  // Command to invoke atomic relaxation or not
  armCmd = new G4UIcmdWithABool ("/grdm/applyARM",this);
  armCmd->SetGuidance("True: ARM is applied; false: no");
  armCmd->SetParameterName("applyARM",true);
  armCmd->SetDefaultValue(true);
  //armCmd->AvailableForStates(G4State_PreInit);

  // Command to set the directional bias (collimation) vector
  colldirCmd = new G4UIcmdWith3Vector("/grdm/decayDirection",this);
  colldirCmd->SetGuidance("Supply the direction vector for decay products");
  colldirCmd->SetParameterName("X","Y","Z",false);

  // Command to set the directional bias (collimation) half angle ("cone")
  collangleCmd = new G4UIcmdWithADoubleAndUnit("/grdm/decayHalfAngle",this);
  collangleCmd->SetGuidance
    ("Supply maximum angle from direction vector for decay products");
  collangleCmd->SetParameterName("halfAngle",false);
  collangleCmd->SetUnitCategory("Angle");

  // This command setup the verbose level of radioactive decay
  verboseCmd = new G4UIcmdWithAnInteger("/grdm/verbose",this);
  verboseCmd->SetGuidance("Set verbose level: 0, 1, 2 or 3");
  verboseCmd->SetParameterName("VerboseLevel",true);
  verboseCmd->SetDefaultValue(1);
  verboseCmd->SetRange("VerboseLevel>=0");

  // Use a user-defined decay datafile for a given isotope
  userDecayDataCmd = new G4UIcommand("/grdm/setRadioactiveDecayFile",this);
  userDecayDataCmd->SetGuidance("Supply user-defined radioactive decay data file");

  G4UIparameter*  Z_para= new G4UIparameter("Z_isotope",'i',true);
  Z_para->SetParameterRange("Z_isotope > 0");
  Z_para->SetGuidance("Z: Charge number of isotope");

  G4UIparameter*  A_para= new G4UIparameter("A_isotope",'i',true);
  A_para->SetParameterRange("A_isotope > 1");
  A_para->SetGuidance("A: mass number of isotope");

  G4UIparameter*  FileName_para= new G4UIparameter("file_name",'s',true);
  FileName_para->SetGuidance("Name of the user data file");

  userDecayDataCmd->SetParameter(Z_para);
  userDecayDataCmd->SetParameter(A_para);
  userDecayDataCmd->SetParameter(FileName_para);

  // Use a user-defined evaporation data file for a given isotope
  userEvaporationDataCmd = new G4UIcommand("/grdm/setPhotoEvaporationFile",this);
  userEvaporationDataCmd->SetGuidance("Supply user-defined photon evaporation data file");

  G4UIparameter* Zpara= new G4UIparameter("Z_isotope",'i',true);
  Zpara->SetParameterRange("Z_isotope > 0");
  Zpara->SetGuidance("Z: Charge number of isotope");

  G4UIparameter* Apara= new G4UIparameter("A_isotope",'i',true);
  Apara->SetParameterRange("A_isotope > 1");
  Apara->SetGuidance("A: mass number of isotope");

  G4UIparameter* FileNamepara= new G4UIparameter("file_name",'s',true);
  FileNamepara->SetGuidance("Name of the user data file");

  userEvaporationDataCmd->SetParameter(Zpara);
  userEvaporationDataCmd->SetParameter(Apara);
  userEvaporationDataCmd->SetParameter(FileNamepara);
}


G4RadioactiveDecayBaseMessenger::~G4RadioactiveDecayBaseMessenger ()
{
  delete grdmDirectory;
  delete nucleuslimitsCmd;
  delete fbetaCmd;
  delete verboseCmd;
  delete avolumeCmd;
  delete deavolumeCmd;
  delete allvolumesCmd;
  delete deallvolumesCmd;
  delete icmCmd;
  delete armCmd;
  delete userDecayDataCmd;
  delete userEvaporationDataCmd;
  delete colldirCmd;
  delete collangleCmd;
}


void
G4RadioactiveDecayBaseMessenger::SetNewValue(G4UIcommand *command, G4String newValues)
{
  if (command==nucleuslimitsCmd) {
    theRadioactiveDecayContainer->
    SetNucleusLimits(nucleuslimitsCmd->GetNewNucleusLimitsValue(newValues));

  } else if (command==fbetaCmd) {
    theRadioactiveDecayContainer->
    SetFBeta(fbetaCmd->GetNewBoolValue(newValues));

  } else if (command==avolumeCmd) {
    theRadioactiveDecayContainer->SelectAVolume(newValues);

  } else if (command==deavolumeCmd) {
    theRadioactiveDecayContainer->DeselectAVolume(newValues);

  } else if (command==allvolumesCmd) {
    theRadioactiveDecayContainer->SelectAllVolumes();

  } else if (command==deallvolumesCmd) {
    theRadioactiveDecayContainer->DeselectAllVolumes();

  } else if (command==verboseCmd) {
    theRadioactiveDecayContainer->
    SetVerboseLevel(verboseCmd->GetNewIntValue(newValues));
  } else if (command==icmCmd) {
    theRadioactiveDecayContainer->
    SetICM(icmCmd->GetNewBoolValue(newValues));

  } else if (command==armCmd) {
    theRadioactiveDecayContainer->
    SetARM(armCmd->GetNewBoolValue(newValues));

  } else if (command ==userDecayDataCmd) {
    G4int Z,A;
    G4String file_name;
    const char* nv = (const char*)newValues;
    std::istringstream is(nv);
    is >> Z >> A >> file_name;
    theRadioactiveDecayContainer->AddUserDecayDataFile(Z,A,file_name);

  } else if (command ==userEvaporationDataCmd) {
    G4int Z,A;
    G4String file_name;
    const char* nv = (const char*)newValues;
    std::istringstream is(nv);
    is >> Z >> A >> file_name;
    G4NuclearLevelData::GetInstance()->AddPrivateData(Z,A,file_name);

  } else if (command==colldirCmd) {
    theRadioactiveDecayContainer->
    SetDecayDirection(colldirCmd->GetNew3VectorValue(newValues));

  } else if (command==collangleCmd) {
    theRadioactiveDecayContainer->
    SetDecayHalfAngle(collangleCmd->GetNewDoubleValue(newValues));
  }
}

