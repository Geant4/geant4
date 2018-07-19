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
//  File:   G4RadioactivationMessenger.cc                                     //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   29 August 2017                                                    //
//  Description: messenger class for biased version of G4RadioactiveDecay.    //
//  Based on the code of F. Lei and P.R. Truscott.                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4RadioactivationMessenger.hh"
#include "G4NuclearLevelData.hh"
#include <sstream>


G4RadioactivationMessenger::G4RadioactivationMessenger(G4Radioactivation* theRadioactivationContainer1)
 :theRadioactivationContainer(theRadioactivationContainer1)
{
  grdmDirectory = new G4UIdirectory("/grdm/");
  grdmDirectory->SetGuidance("Controls the biased version of radioactive decay");
/*
  // Command to define the limits on nuclei the RDM will treat
  nucleuslimitsCmd = new G4UIcmdWithNucleusLimits("/grdm/nucleusLimits",this);
  nucleuslimitsCmd->SetGuidance 
    ("Set the atomic weight and number limits for the RDM.");
  nucleuslimitsCmd->SetParameterName("aMin","aMax","zMin","zMax",true);

  // Command to contols whether beta decay will be treated faithfully or 
  // in fast mode (no longer used) 
  fbetaCmd = new G4UIcmdWithABool ("/grdm/fBeta",this);
  fbetaCmd->SetGuidance("false: use 3-body decay, true: use histogram method");
  fbetaCmd->SetParameterName("fBeta",true);
  fbetaCmd->SetDefaultValue(false);

  // Command to select a logical volume for RDM
  avolumeCmd = new G4UIcmdWithAString("/grdm/selectVolume",this);
  avolumeCmd->SetGuidance 
    ("Suppply a logical volumes name to add it to the RDM apply list");
  avolumeCmd->SetParameterName("aVolume",false);

  // Command to de-select a logical volume for RDM
  deavolumeCmd = new G4UIcmdWithAString("/grdm/deselectVolume",this);
  deavolumeCmd->SetGuidance 
    ("Suppply a logical volumes name to remove it from the RDM apply list");
  deavolumeCmd->SetParameterName("aVolume",false);

  // Command to select all logical volumes for RDM
  allvolumesCmd = new G4UIcmdWithoutParameter("/grdm/allVolumes",this);
  allvolumesCmd->SetGuidance 
    (" apply RDM to all logical volumes. No parameter required.");

  // Command to de-selete a logical volume for RDM.
  deallvolumesCmd = new G4UIcmdWithoutParameter("/grdm/noVolumes",this);
  deallvolumesCmd->SetGuidance(" RDM is not applied to any logical volumes");
*/
  // Command to use branching ratio biasing or not
  brbiasCmd = new G4UIcmdWithABool ("/grdm/BRbias",this);
  brbiasCmd->SetGuidance("false: no biasing; true: all branches are treated as equal");
  brbiasCmd->SetParameterName("BRBias",true);
  brbiasCmd->SetDefaultValue(true);
/*
  // Command to apply internal conversion or not
  icmCmd = new G4UIcmdWithABool ("/grdm/applyICM",this);
  icmCmd->SetGuidance("True: ICM is applied; false: no");
  icmCmd->SetParameterName("applyICM",true);
  icmCmd->SetDefaultValue(true);

  // Command to apply atomic relaxation or not
  armCmd = new G4UIcmdWithABool ("/grdm/applyARM",this);
  armCmd->SetGuidance("True: ARM is applied; false: no");
  armCmd->SetParameterName("applyARM",true);
  armCmd->SetDefaultValue(true);
*/
  // Command to set the half-life thresold for isomer production
  hlthCmd = new G4UIcmdWithADoubleAndUnit("/grdm/hlThreshold",this);
  hlthCmd->SetGuidance("Set the h-l threshold for isomer production");
  hlthCmd->SetParameterName("hlThreshold",false);
  hlthCmd->SetUnitCategory("Time");

  // Command to define the incident particle source time profile
  sourcetimeprofileCmd = new G4UIcmdWithAString("/grdm/sourceTimeProfile",this);
  sourcetimeprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the source particle time profile");
  sourcetimeprofileCmd->SetParameterName("STimeProfile",true);
  sourcetimeprofileCmd->SetDefaultValue("source.data");

  // Command to define the incident particle source time profile
  decaybiasprofileCmd = new G4UIcmdWithAString("/grdm/decayBiasProfile",this);
  decaybiasprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the decay bias time profile");
  decaybiasprofileCmd->SetParameterName("DBiasProfile",true);
  decaybiasprofileCmd->SetDefaultValue("bias.data");
/*
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
*/
  // Command to set nuclei spliting parameter
  splitnucleiCmd = new G4UIcmdWithAnInteger("/grdm/splitNuclei",this);
  splitnucleiCmd->SetGuidance("Set number of spliting for the isotopes.");
  splitnucleiCmd->SetParameterName("NSplit",true);
  splitnucleiCmd->SetDefaultValue(1);
  splitnucleiCmd->SetRange("NSplit>=1");
/*
  // Command to set verbosity
  verboseCmd = new G4UIcmdWithAnInteger("/grdm/verbose",this);
  verboseCmd->SetGuidance("Set verbose level: 0, 1, 2 or 3");
  verboseCmd->SetParameterName("VerboseLevel",true);
  verboseCmd->SetDefaultValue(1);
  verboseCmd->SetRange("VerboseLevel>=0");

  //This commansd allows the user to define its own decay datafile for
  // a given isotope
  //
  userDecayDataCmd = new G4UIcommand("/grdm/setRadioactiveDecayFile",this);
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

  // Command allowing user-defined evaporation data file to be used
  // a given isotope

  userEvaporationDataCmd = new G4UIcommand("/grdm/setPhotoEvaporationFile",this);
  userEvaporationDataCmd->SetParameter(Z_para);
  userEvaporationDataCmd->SetParameter(A_para);
  userEvaporationDataCmd->SetParameter(FileName_para);
*/
}


G4RadioactivationMessenger::~G4RadioactivationMessenger()
{
  delete grdmDirectory;
//  delete nucleuslimitsCmd;
  delete sourcetimeprofileCmd;
  delete decaybiasprofileCmd;
//  delete fbetaCmd;
  delete brbiasCmd;
  delete splitnucleiCmd;
//  delete verboseCmd;
//  delete avolumeCmd;
//  delete deavolumeCmd;
//  delete allvolumesCmd;
//  delete deallvolumesCmd;
//  delete icmCmd;
//  delete armCmd;
  delete hlthCmd;
//  delete userDecayDataCmd;
//  delete userEvaporationDataCmd;
//  delete colldirCmd;
//  delete collangleCmd;
}


void G4RadioactivationMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
//  if (command==nucleuslimitsCmd) {theRadioactivationContainer->
//    SetNucleusLimits(nucleuslimitsCmd->GetNewNucleusLimitsValue(newValues));
/*
  } else if (command==fbetaCmd) {theRadioactivationContainer->
    SetFBeta(fbetaCmd->GetNewBoolValue(newValues));

  } else if (command==avolumeCmd) {theRadioactivationContainer->
    SelectAVolume(newValues);

  } else if (command==deavolumeCmd) {theRadioactivationContainer->
    DeselectAVolume(newValues);

  } else if (command==allvolumesCmd) {theRadioactivationContainer->
    SelectAllVolumes();

  } else if (command==deallvolumesCmd) {theRadioactivationContainer->
    DeselectAllVolumes();
*/
  if (command==brbiasCmd) {theRadioactivationContainer->
    SetBRBias(brbiasCmd->GetNewBoolValue(newValues));

  } else if (command==sourcetimeprofileCmd) {theRadioactivationContainer->
    SetSourceTimeProfile(newValues);

  } else if (command==decaybiasprofileCmd) {theRadioactivationContainer->
    SetDecayBias(newValues);

  } else if (command==splitnucleiCmd) {theRadioactivationContainer->
    SetSplitNuclei(splitnucleiCmd->GetNewIntValue(newValues));
/*
  } else if (command==verboseCmd) {theRadioactivationContainer->
    SetVerboseLevel(verboseCmd->GetNewIntValue(newValues));

  } else if (command==icmCmd ) {theRadioactivationContainer->
    SetICM(icmCmd->GetNewBoolValue(newValues));

  } else if (command==armCmd ) {theRadioactivationContainer->
    SetARM(armCmd->GetNewBoolValue(newValues));
*/
  } else if (command==hlthCmd ) {theRadioactivationContainer->
    SetHLThreshold(hlthCmd->GetNewDoubleValue(newValues));
/*
  } else if (command ==userDecayDataCmd) {
    G4int Z,A;
    G4String file_name;
    const char* nv = (const char*)newValues;
    std::istringstream is(nv);
    is >> Z >> A >> file_name;
    theRadioactivationContainer->AddUserDecayDataFile(Z,A,file_name);

  } else if (command ==userEvaporationDataCmd) {
    G4int Z,A;
    G4String file_name;
    const char* nv = (const char*)newValues;
    std::istringstream is(nv);
    is >> Z >> A >> file_name;
    G4NuclearLevelData::GetInstance()->AddPrivateData(Z,A,file_name);

  } else if (command==colldirCmd) {theRadioactivationContainer->
    SetDecayDirection(colldirCmd->GetNew3VectorValue(newValues));

  } else if (command==collangleCmd) {theRadioactivationContainer->
    SetDecayHalfAngle(collangleCmd->GetNewDoubleValue(newValues));
*/
  }
}

