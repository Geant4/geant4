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
#include "G4HadronicException.hh"


G4RadioactiveDecayBaseMessenger::G4RadioactiveDecayBaseMessenger 
(G4RadioactiveDecayBase* theRadioactiveDecayContainer1)
:theRadioactiveDecayContainer(theRadioactiveDecayContainer1)
{
  // main directory for control of the RDM
  old_grdmDirectory = new G4UIdirectory("/grdm/");                                     // To be removed in G4 11.0
  old_grdmDirectory->SetGuidance("Controls for the Radioactive Decay Module.");

  rdmDirectory = new G4UIdirectory("/process/had/rdm/");
  rdmDirectory->SetGuidance("Controls for the Radioactive Decay Module.");
  
  // Command to define the limits on nucleus the RDM will treat.
  old_nucleuslimitsCmd = new G4UIcmdWithNucleusLimits("/grdm/nucleusLimits",this);     // To be removed in G4 11.0
  old_nucleuslimitsCmd->SetGuidance("Set the atomic weight and number limits for the RDM.");
  old_nucleuslimitsCmd->SetParameterName("OldAMin","OldAMax","OldZMin","OldZMax",true);

  nucleuslimitsCmd = new G4UIcmdWithNucleusLimits("/process/had/rdm/nucleusLimits",this);
  nucleuslimitsCmd->SetGuidance("Set the atomic weight and number limits for the RDM.");
  nucleuslimitsCmd->SetParameterName("AMin","AMax","ZMin","ZMax",true);
  
  // Select a logical volume for RDM
  old_avolumeCmd = new G4UIcmdWithAString("/grdm/selectVolume",this);                  // To be removed in G4 11.0
  old_avolumeCmd->SetGuidance("Supply a logical volumes name to add it to the RDM apply list");
  old_avolumeCmd->SetParameterName("OldAVolume",false);

  avolumeCmd = new G4UIcmdWithAString("/process/had/rdm/selectVolume",this);
  avolumeCmd->SetGuidance("Supply a logical volumes name to add it to the RDM apply list");
  avolumeCmd->SetParameterName("AVolume",false);
  
  // De-select a logical volume for RDM
  old_deavolumeCmd = new G4UIcmdWithAString("/grdm/deselectVolume",this);              // To be removed in G4 11.0
  old_deavolumeCmd->SetGuidance("Supply a logical volumes name to remove it from the RDM apply list");
  old_deavolumeCmd->SetParameterName("OldAVolume",false);

  deavolumeCmd = new G4UIcmdWithAString("/process/had/rdm/deselectVolume",this);
  deavolumeCmd->SetGuidance("Supply a logical volumes name to remove it from the RDM apply list");
  deavolumeCmd->SetParameterName("AVolume",false);
  
  // Select all logical volumes for RDM
  old_allvolumesCmd = new G4UIcmdWithoutParameter("/grdm/allVolumes",this);            // To be removed in G4 11.0
  old_allvolumesCmd->SetGuidance(" apply RDM to all logical volumes. No parameter required.");
  //  old_allvolumeCmd->SetParameterName("OldAddAVolume",true);

  allvolumesCmd = new G4UIcmdWithoutParameter("/process/had/rdm/allVolumes",this);
  allvolumesCmd->SetGuidance(" apply RDM to all logical volumes. No parameter required.");
  //  allvolumeCmd->SetParameterName("AddAVolume",true);
  
  // De-select all logical volumes for RDM
  old_deallvolumesCmd = new G4UIcmdWithoutParameter("/grdm/noVolumes",this);           // To be removed in G4 11.0
  old_deallvolumesCmd->SetGuidance(" RDM is not applied to any logical volumes");
  //  old_deallvolumesCmd->SetParameterName("OldRemoveAVolume",true);
  
  deallvolumesCmd = new G4UIcmdWithoutParameter("/process/had/rdm/noVolumes",this);
  deallvolumesCmd->SetGuidance(" RDM is not applied to any logical volumes");
  //  deallvolumesCmd->SetParameterName("RemoveAVolume",true);

  // Command to invoke internal conversion or not
  old_icmCmd = new G4UIcmdWithABool("/grdm/applyICM",this);                           // To be removed in G4 11.0
  old_icmCmd->SetGuidance("Command not active; kept for backward compatibility.");
  old_icmCmd->SetGuidance("Internal conversion is always turned on.");
  old_icmCmd->SetParameterName("OldApplyICM",true);
  old_icmCmd->SetDefaultValue(true);

  icmCmd = new G4UIcmdWithABool("/process/had/rdm/applyICM",this);
  icmCmd->SetGuidance("Command not active; kept for backward compatibility.");
  icmCmd->SetGuidance("Internal conversion is always turned on.");
  icmCmd->SetParameterName("ApplyICM",true);
  icmCmd->SetDefaultValue(true);
  
  // Command to invoke atomic relaxation or not
  old_armCmd = new G4UIcmdWithABool("/grdm/applyARM",this);                           // To be removed in G4 11.0
  old_armCmd->SetGuidance("True: ARM is applied; false: no");
  old_armCmd->SetParameterName("OldApplyARM",true);
  old_armCmd->SetDefaultValue(true);
  //old_armCmd->AvailableForStates(G4State_PreInit);

  armCmd = new G4UIcmdWithABool("/process/had/rdm/applyARM",this);
  armCmd->SetGuidance("True: ARM is applied; false: no");
  armCmd->SetParameterName("ApplyARM",true);
  armCmd->SetDefaultValue(true);
  //armCmd->AvailableForStates(G4State_PreInit);
  
  // Command to set the directional bias (collimation) vector
  old_colldirCmd = new G4UIcmdWith3Vector("/grdm/decayDirection",this);                // To be removed in G4 11.0
  old_colldirCmd->SetGuidance("Supply the direction vector for decay products");
  old_colldirCmd->SetParameterName("OldX","OldY","OldZ",false);

  colldirCmd = new G4UIcmdWith3Vector("/process/had/rdm/decayDirection",this);
  colldirCmd->SetGuidance("Supply the direction vector for decay products");
  colldirCmd->SetParameterName("X","Y","Z",false);
  
  // Command to set the directional bias (collimation) half angle ("cone")
  old_collangleCmd = new G4UIcmdWithADoubleAndUnit("/grdm/decayHalfAngle",this);       // To be removed in G4 11.0
  old_collangleCmd->SetGuidance("Supply maximum angle from direction vector for decay products");
  old_collangleCmd->SetParameterName("OldHalfAngle",false);
  old_collangleCmd->SetUnitCategory("Angle");

  collangleCmd = new G4UIcmdWithADoubleAndUnit("/process/had/rdm/decayHalfAngle",this);
  collangleCmd->SetGuidance("Supply maximum angle from direction vector for decay products");
  collangleCmd->SetParameterName("HalfAngle",false);
  collangleCmd->SetUnitCategory("Angle");
  
  // This command setup the verbose level of radioactive decay
  old_verboseCmd = new G4UIcmdWithAnInteger("/grdm/verbose",this);                     // To be removed in G4 11.0
  old_verboseCmd->SetGuidance("Set verbose level: 0, 1, 2 or 3");
  old_verboseCmd->SetParameterName("OldVerboseLevel",true);
  old_verboseCmd->SetDefaultValue(1);
  old_verboseCmd->SetRange("OldVerboseLevel>=0");

  verboseCmd = new G4UIcmdWithAnInteger("/process/had/rdm/verbose",this);
  verboseCmd->SetGuidance("Set verbose level: 0, 1, 2 or 3");
  verboseCmd->SetParameterName("VerboseLevel",true);
  verboseCmd->SetDefaultValue(1);
  verboseCmd->SetRange("VerboseLevel>=0");
  
  // Use a user-defined decay datafile for a given isotope
  old_userDecayDataCmd = new G4UIcommand("/grdm/setRadioactiveDecayFile",this);        // To be removed in G4 11.0
  old_userDecayDataCmd->SetGuidance("Supply user-defined radioactive decay data file");
  G4UIparameter* old_Z_para= new G4UIparameter("Z_isotope",'i',true);
  old_Z_para->SetParameterRange("Z_isotope > 0");
  old_Z_para->SetGuidance("Z: Charge number of isotope");
  G4UIparameter* old_A_para= new G4UIparameter("A_isotope",'i',true);
  old_A_para->SetParameterRange("A_isotope > 1");
  old_A_para->SetGuidance("A: mass number of isotope");
  G4UIparameter* old_FileName_para= new G4UIparameter("file_name",'s',true);
  old_FileName_para->SetGuidance("Name of the user data file");
  old_userDecayDataCmd->SetParameter(old_Z_para);
  old_userDecayDataCmd->SetParameter(old_A_para);
  old_userDecayDataCmd->SetParameter(old_FileName_para);
  
  userDecayDataCmd = new G4UIcommand("/process/had/rdm/setRadioactiveDecayFile",this);
  userDecayDataCmd->SetGuidance("Supply user-defined radioactive decay data file");
  G4UIparameter* Z_para= new G4UIparameter("Z_isotope",'i',true);
  Z_para->SetParameterRange("Z_isotope > 0");
  Z_para->SetGuidance("Z: Charge number of isotope");
  G4UIparameter* A_para= new G4UIparameter("A_isotope",'i',true);
  A_para->SetParameterRange("A_isotope > 1");
  A_para->SetGuidance("A: mass number of isotope");
  G4UIparameter* FileName_para= new G4UIparameter("file_name",'s',true);
  FileName_para->SetGuidance("Name of the user data file");
  userDecayDataCmd->SetParameter(Z_para);
  userDecayDataCmd->SetParameter(A_para);
  userDecayDataCmd->SetParameter(FileName_para);
  
  // Use a user-defined evaporation data file for a given isotope
  old_userEvaporationDataCmd = new G4UIcommand("/grdm/setPhotoEvaporationFile",this);  // To be removed in G4 11.0
  old_userEvaporationDataCmd->SetGuidance("Supply user-defined photon evaporation data file");
  G4UIparameter* old_Zpara= new G4UIparameter("Z_isotope",'i',true);
  old_Zpara->SetParameterRange("Z_isotope > 0");
  old_Zpara->SetGuidance("Z: Charge number of isotope");
  G4UIparameter* old_Apara= new G4UIparameter("A_isotope",'i',true);
  old_Apara->SetParameterRange("A_isotope > 1");
  old_Apara->SetGuidance("A: mass number of isotope");
  G4UIparameter* old_FileNamepara= new G4UIparameter("file_name",'s',true);
  old_FileNamepara->SetGuidance("Name of the user data file");
  old_userEvaporationDataCmd->SetParameter(old_Zpara);
  old_userEvaporationDataCmd->SetParameter(old_Apara);
  old_userEvaporationDataCmd->SetParameter(old_FileNamepara);
  
  userEvaporationDataCmd = new G4UIcommand("/process/had/rdm/setPhotoEvaporationFile",this);
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
  delete old_grdmDirectory;           // To be removed in G4 11.0
  delete      rdmDirectory;
  delete old_nucleuslimitsCmd;        // To be removed in G4 11.0
  delete     nucleuslimitsCmd;
  delete old_verboseCmd;              // To be removed in G4 11.0
  delete     verboseCmd;
  delete old_avolumeCmd;              // To be removed in G4 11.0
  delete     avolumeCmd;
  delete old_deavolumeCmd;            // To be removed in G4 11.0
  delete     deavolumeCmd;
  delete old_allvolumesCmd;           // To be removed in G4 11.0
  delete     allvolumesCmd;
  delete old_deallvolumesCmd;         // To be removed in G4 11.0
  delete     deallvolumesCmd;
  delete old_icmCmd;                  // To be removed in G4 11.0
  delete     icmCmd;
  delete old_armCmd;                  // To be removed in G4 11.0
  delete     armCmd;
  delete old_userDecayDataCmd;        // To be removed in G4 11.0
  delete     userDecayDataCmd;                                        //<---
  delete old_userEvaporationDataCmd;  // To be removed in G4 11.0
  delete     userEvaporationDataCmd;
  delete old_colldirCmd;              // To be removed in G4 11.0
  delete     colldirCmd;
  delete old_collangleCmd;            // To be removed in G4 11.0 
  delete     collangleCmd;
}


void
G4RadioactiveDecayBaseMessenger::SetNewValue(G4UIcommand *command, G4String newValues)
{
  // Old commands to be removed in G4 11.0
  if ( command == old_nucleuslimitsCmd ) {
    theRadioactiveDecayContainer->
    SetNucleusLimits( old_nucleuslimitsCmd->GetNewNucleusLimitsValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/nucleusLimits in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecayBaseMessenger", "HAD_RDM_671", JustWarning, ed );
  } else if ( command == old_avolumeCmd ) {
    theRadioactiveDecayContainer->SelectAVolume( newValues );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/selectVolume in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecayBaseMessenger", "HAD_RDM_672", JustWarning, ed );
  } else if ( command == old_deavolumeCmd ) {
    theRadioactiveDecayContainer->DeselectAVolume( newValues );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/deselectVolume in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecayBaseMessenger", "HAD_RDM_673", JustWarning, ed );
  } else if ( command == old_allvolumesCmd ) {
    theRadioactiveDecayContainer->SelectAllVolumes();
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/allVolumes in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecayBaseMessenger", "HAD_RDM_674", JustWarning, ed );
  } else if ( command == old_deallvolumesCmd ) {
    theRadioactiveDecayContainer->DeselectAllVolumes();
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/noVolumes in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecayBaseMessenger", "HAD_RDM_675", JustWarning, ed );
  } else if ( command == old_verboseCmd ) {
    theRadioactiveDecayContainer->SetVerboseLevel( old_verboseCmd->GetNewIntValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/verbose in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecayBaseMessenger", "HAD_RDM_676", JustWarning, ed );
  } else if ( command == old_icmCmd ) {
    theRadioactiveDecayContainer->SetICM( old_icmCmd->GetNewBoolValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/applyICM in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecayBaseMessenger", "HAD_RDM_677", JustWarning, ed );
  } else if ( command == old_armCmd ) {
    theRadioactiveDecayContainer->SetARM( old_armCmd->GetNewBoolValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/applyARM in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecayBaseMessenger", "HAD_RDM_678", JustWarning, ed );
  } else if ( command == old_userDecayDataCmd ) {
    G4int Z,A;
    G4String file_name;
    const char* nv = (const char*)newValues;
    std::istringstream is(nv);
    is >> Z >> A >> file_name;
    theRadioactiveDecayContainer->AddUserDecayDataFile(Z,A,file_name);
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/setRadioactiveDecayFile in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecayBaseMessenger", "HAD_RDM_679", JustWarning, ed );
  } else if ( command == old_userEvaporationDataCmd ) {
    G4int Z,A;
    G4String file_name;
    const char* nv = (const char*)newValues;
    std::istringstream is(nv);
    is >> Z >> A >> file_name;
    G4NuclearLevelData::GetInstance()->AddPrivateData(Z,A,file_name);
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/setPhotoEvaporationFile in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecayBaseMessenger", "HAD_RDM_680", JustWarning, ed );
  } else if ( command == old_colldirCmd ) {
    theRadioactiveDecayContainer->SetDecayDirection( old_colldirCmd->GetNew3VectorValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/decayDirection in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecayBaseMessenger", "HAD_RDM_681", JustWarning, ed );
  } else if ( command == old_collangleCmd ) {
    theRadioactiveDecayContainer->SetDecayHalfAngle( old_collangleCmd->GetNewDoubleValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/decayHalfAngle in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecayBaseMessenger", "HAD_RDM_682", JustWarning, ed );
  }

  // New commands
  if ( command == nucleuslimitsCmd ) {
    theRadioactiveDecayContainer->
    SetNucleusLimits( nucleuslimitsCmd->GetNewNucleusLimitsValue( newValues ) );
  } else if ( command == avolumeCmd ) {
    theRadioactiveDecayContainer->SelectAVolume( newValues );
  } else if ( command == deavolumeCmd ) {
    theRadioactiveDecayContainer->DeselectAVolume( newValues );
  } else if ( command == allvolumesCmd ) {
    theRadioactiveDecayContainer->SelectAllVolumes();
  } else if ( command == deallvolumesCmd ) {
    theRadioactiveDecayContainer->DeselectAllVolumes();
  } else if ( command == verboseCmd ) {
    theRadioactiveDecayContainer->SetVerboseLevel( verboseCmd->GetNewIntValue( newValues ) );
  } else if ( command == icmCmd ) {
    theRadioactiveDecayContainer->SetICM( icmCmd->GetNewBoolValue( newValues ) );
  } else if ( command == armCmd ) {
    theRadioactiveDecayContainer->SetARM( armCmd->GetNewBoolValue( newValues ) );
  } else if ( command == userDecayDataCmd ) {
    G4int Z,A;
    G4String file_name;
    const char* nv = (const char*)newValues;
    std::istringstream is(nv);
    is >> Z >> A >> file_name;
    theRadioactiveDecayContainer->AddUserDecayDataFile(Z,A,file_name);
  } else if ( command == userEvaporationDataCmd ) {
    G4int Z,A;
    G4String file_name;
    const char* nv = (const char*)newValues;
    std::istringstream is(nv);
    is >> Z >> A >> file_name;
    G4NuclearLevelData::GetInstance()->AddPrivateData(Z,A,file_name);
  } else if ( command == colldirCmd ) {
    theRadioactiveDecayContainer->SetDecayDirection( colldirCmd->GetNew3VectorValue( newValues ) );
  } else if ( command == collangleCmd ) {
    theRadioactiveDecayContainer->SetDecayHalfAngle( collangleCmd->GetNewDoubleValue( newValues ) );
  }  
}

