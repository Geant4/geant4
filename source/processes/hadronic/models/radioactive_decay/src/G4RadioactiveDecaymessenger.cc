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

#include "G4RadioactiveDecaymessenger.hh"
#include "G4NuclearLevelData.hh"
#include <sstream>
#include "G4HadronicException.hh"


G4RadioactiveDecaymessenger::G4RadioactiveDecaymessenger 
(G4RadioactiveDecay* theRadioactiveDecayContainer1)
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
  
  // The next command contols whether the decay will be treated analoguely or 
  // with variance reduction
  old_analoguemcCmd = new G4UIcmdWithABool("/grdm/analogueMC",this);                  // To be removed in G4 11.0
  old_analoguemcCmd->SetGuidance("false: variance reduction method; true: analogue method");
  old_analoguemcCmd->SetParameterName("OldAnalogueMC",true);
  old_analoguemcCmd->SetDefaultValue(true);

  analoguemcCmd = new G4UIcmdWithABool("/process/had/rdm/analogueMC",this);
  analoguemcCmd->SetGuidance("false: variance reduction method; true: analogue method");
  analoguemcCmd->SetParameterName("AnalogueMC",true);
  analoguemcCmd->SetDefaultValue(true);
  
  // The next command contols whether beta decay will be treated faithfully or in fast mode
  old_fbetaCmd = new G4UIcmdWithABool("/grdm/fBeta",this);                            // To be removed in G4 11.0
  old_fbetaCmd->SetGuidance("false: use 3-body decay, true: use histogram method");
  old_fbetaCmd->SetParameterName("OldFBeta",true);
  old_fbetaCmd->SetDefaultValue(false);

  fbetaCmd = new G4UIcmdWithABool("/process/had/rdm/fBeta",this);
  fbetaCmd->SetGuidance("false: use 3-body decay, true: use histogram method");
  fbetaCmd->SetParameterName("FBeta",true);
  fbetaCmd->SetDefaultValue(false);
 
  // Command to selete a logical volume for RDM.
  old_avolumeCmd = new G4UIcmdWithAString("/grdm/selectVolume",this);                  // To be removed in G4 11.0
  old_avolumeCmd->SetGuidance("Supply a logical volumes name to add it to the RDM apply list");
  old_avolumeCmd->SetParameterName("OldAVolume",false);

  avolumeCmd = new G4UIcmdWithAString("/process/had/rdm/selectVolume",this);
  avolumeCmd->SetGuidance("Supply a logical volumes name to add it to the RDM apply list");
  avolumeCmd->SetParameterName("AVolume",false);
  
  // Command to de-selete a logical volume for RDM.
  old_deavolumeCmd = new G4UIcmdWithAString("/grdm/deselectVolume",this);              // To be removed in G4 11.0
  old_deavolumeCmd->SetGuidance("Supply a logical volumes name to remove it from the RDM apply list");
  old_deavolumeCmd->SetParameterName("OldAVolume",false);

  deavolumeCmd = new G4UIcmdWithAString("/process/had/rdm/deselectVolume",this);
  deavolumeCmd->SetGuidance("Supply a logical volumes name to remove it from the RDM apply list");
  deavolumeCmd->SetParameterName("AVolume",false);
  
  // Command to selete all logical volumes for RDM.
  old_allvolumesCmd = new G4UIcmdWithoutParameter("/grdm/allVolumes",this);            // To be removed in G4 11.0
  old_allvolumesCmd->SetGuidance(" apply RDM to all logical volumes. No parameter required.");
  //  old_allvolumeCmd->SetParameterName("OldAddAVolume",true);
  
  allvolumesCmd = new G4UIcmdWithoutParameter("/process/had/rdm/allVolumes",this);
  allvolumesCmd->SetGuidance(" apply RDM to all logical volumes. No parameter required.");
  //  allvolumeCmd->SetParameterName("AddAVolume",true);

  // Command to de-selete a logical volume for RDM.
  old_deallvolumesCmd = new G4UIcmdWithoutParameter("/grdm/noVolumes",this);           // To be removed in G4 11.0
  old_deallvolumesCmd->SetGuidance(" RDM is not applied to any logical volumes");
  //  old_deallvolumesCmd->SetParameterName("OldRemoveAVolume",true);

  deallvolumesCmd = new G4UIcmdWithoutParameter("/process/had/rdm/noVolumes",this);
  deallvolumesCmd->SetGuidance(" RDM is not applied to any logical volumes");
  //  deallvolumesCmd->SetParameterName("RemoveAVolume",true);
  
  // The next command contols whether the branching ratio biasing will be applied or not
  old_brbiasCmd = new G4UIcmdWithABool("/grdm/BRbias",this);                          // To be removed in G4 11.0
  old_brbiasCmd->SetGuidance("false: no biasing; true: all branches are treated as equal");
  old_brbiasCmd->SetParameterName("OldBRBias",true);
  old_brbiasCmd->SetDefaultValue(true);

  brbiasCmd = new G4UIcmdWithABool("/process/had/rdm/BRbias",this);
  brbiasCmd->SetGuidance("false: no biasing; true: all branches are treated as equal");
  brbiasCmd->SetParameterName("BRBias",true);
  brbiasCmd->SetDefaultValue(true);
  
  // Command controls whether ICM will be applied or not
  old_icmCmd = new G4UIcmdWithABool("/grdm/applyICM",this);                           // To be removed in G4 11.0
  old_icmCmd->SetGuidance("True: ICM is applied; false: no");
  old_icmCmd->SetParameterName("OldApplyICM",true);
  old_icmCmd->SetDefaultValue(true);

  icmCmd = new G4UIcmdWithABool("/process/had/rdm/applyICM",this);
  icmCmd->SetGuidance("True: ICM is applied; false: no");
  icmCmd->SetParameterName("ApplyICM",true);
  icmCmd->SetDefaultValue(true);
  
  // Command contols whether ARM will be applied or not
  old_armCmd = new G4UIcmdWithABool("/grdm/applyARM",this);                           // To be removed in G4 11.0
  old_armCmd->SetGuidance("True: ARM is applied; false: no");
  old_armCmd->SetParameterName("OldApplyARM",true);
  old_armCmd->SetDefaultValue(true);
  //  old_armCmd->AvailableForStates(G4State_PreInit);

  armCmd = new G4UIcmdWithABool("/process/had/rdm/applyARM",this);
  armCmd->SetGuidance("True: ARM is applied; false: no");
  armCmd->SetParameterName("ApplyARM",true);
  armCmd->SetDefaultValue(true);
  //  armCmd->AvailableForStates(G4State_PreInit);

  // Command to set the h-l thresold for isomer production
  old_hlthCmd = new G4UIcmdWithADoubleAndUnit("/grdm/hlThreshold",this);               // To be removed in G4 11.0
  old_hlthCmd->SetGuidance("Set the h-l threshold for isomer production");
  old_hlthCmd->SetParameterName("OldHlThreshold",false);
  //  old_hlthCmd->SetRange("OldHlThreshold>0.");
  old_hlthCmd->SetUnitCategory("Time");
  //  old_hlthCmd->AvailableForStates(G4State_PreInit);

  hlthCmd = new G4UIcmdWithADoubleAndUnit("/process/had/rdm/hlThreshold",this);
  hlthCmd->SetGuidance("Set the h-l threshold for isomer production");
  hlthCmd->SetParameterName("HlThreshold",false);
  //  hlthCmd->SetRange("HlThreshold>0.");
  hlthCmd->SetUnitCategory("Time");
  //  hlthCmd->AvailableForStates(G4State_PreInit);
 
  // Command to define the incident particle source time profile.
  old_sourcetimeprofileCmd = new G4UIcmdWithAString("/grdm/sourceTimeProfile",this);   // To be removed in G4 11.0
  old_sourcetimeprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the source particle time profile");
  old_sourcetimeprofileCmd->SetParameterName("OldSTimeProfile",true);
  old_sourcetimeprofileCmd->SetDefaultValue("source.data");

  sourcetimeprofileCmd = new G4UIcmdWithAString("/process/had/rdm/sourceTimeProfile",this);
  sourcetimeprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the source particle time profile");
  sourcetimeprofileCmd->SetParameterName("STimeProfile",true);
  sourcetimeprofileCmd->SetDefaultValue("source.data");

  // Command to define the incident particle source time profile.
  old_decaybiasprofileCmd = new G4UIcmdWithAString("/grdm/decayBiasProfile",this);     // To be removed in G4 11.0
  old_decaybiasprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the decay bias time profile");
  old_decaybiasprofileCmd->SetParameterName("OldDBiasProfile",true);
  old_decaybiasprofileCmd->SetDefaultValue("bias.data");

  decaybiasprofileCmd = new G4UIcmdWithAString("/process/had/rdm/decayBiasProfile",this);
  decaybiasprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the decay bias time profile");
  decaybiasprofileCmd->SetParameterName("DBiasProfile",true);
  decaybiasprofileCmd->SetDefaultValue("bias.data");
  
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
  
  // This command setup the nuclei spliting parameter
  old_splitnucleiCmd = new G4UIcmdWithAnInteger("/grdm/splitNuclei",this);             // To be removed in G4 11.0
  old_splitnucleiCmd->SetGuidance("Set number of spliting for the isotopes.");
  old_splitnucleiCmd->SetParameterName("OldNSplit",true);
  old_splitnucleiCmd->SetDefaultValue(1);
  old_splitnucleiCmd->SetRange("OldNSplit>=1");

  splitnucleiCmd = new G4UIcmdWithAnInteger("/process/had/rdm/splitNuclei",this);
  splitnucleiCmd->SetGuidance("Set number of spliting for the isotopes.");
  splitnucleiCmd->SetParameterName("NSplit",true);
  splitnucleiCmd->SetDefaultValue(1);
  splitnucleiCmd->SetRange("NSplit>=1");

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
  
  // This command allows the user to define its own decay datafile for a given isotope
  old_userDecayDataCmd = new G4UIcommand("/grdm/setRadioactiveDecayFile",this);        // To be removed in G4 11.0
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
  
  // This commands allows the user to define its own evaporation data file for a given isotope
  old_userEvaporationDataCmd = new G4UIcommand("/grdm/setPhotoEvaporationFile",this);  // To be removed in G4 11.0
  G4UIparameter* old_Z_para_ = new G4UIparameter("Z_isotope",'i',true);
  old_Z_para_->SetParameterRange("Z_isotope > 0");
  old_Z_para_->SetGuidance("Z: Charge number of isotope");
  G4UIparameter* old_A_para_ = new G4UIparameter("A_isotope",'i',true);
  old_A_para_->SetParameterRange("A_isotope > 1");
  old_A_para_->SetGuidance("A: mass number of isotope");
  G4UIparameter* old_FileName_para_ = new G4UIparameter("file_name",'s',true);
  old_FileName_para_->SetGuidance("Name of the user data file");
  old_userEvaporationDataCmd->SetParameter(old_Z_para_);
  old_userEvaporationDataCmd->SetParameter(old_A_para_);
  old_userEvaporationDataCmd->SetParameter(old_FileName_para_);

  userEvaporationDataCmd = new G4UIcommand("/process/had/rdm/setPhotoEvaporationFile",this);
  G4UIparameter* Z_para_ = new G4UIparameter("Z_isotope",'i',true);
  Z_para_->SetParameterRange("Z_isotope > 0");
  Z_para_->SetGuidance("Z: Charge number of isotope");
  G4UIparameter* A_para_ = new G4UIparameter("A_isotope",'i',true);
  A_para_->SetParameterRange("A_isotope > 1");
  A_para_->SetGuidance("A: mass number of isotope");
  G4UIparameter* FileName_para_ = new G4UIparameter("file_name",'s',true);
  FileName_para_->SetGuidance("Name of the user data file");
  userEvaporationDataCmd->SetParameter(Z_para_);
  userEvaporationDataCmd->SetParameter(A_para_);
  userEvaporationDataCmd->SetParameter(FileName_para_);
  
}
////////////////////////////////////////////////////////////////////////////////
//
G4RadioactiveDecaymessenger::~G4RadioactiveDecaymessenger ()
{
  delete old_grdmDirectory;           // To be removed in G4 11.0
  delete     rdmDirectory;
  delete old_nucleuslimitsCmd;        // To be removed in G4 11.0
  delete     nucleuslimitsCmd;
  delete old_sourcetimeprofileCmd;    // To be removed in G4 11.0
  delete     sourcetimeprofileCmd;
  delete old_decaybiasprofileCmd;     // To be removed in G4 11.0
  delete     decaybiasprofileCmd;
  delete old_analoguemcCmd;           // To be removed in G4 11.0
  delete     analoguemcCmd;
  delete old_fbetaCmd;                // To be removed in G4 11.0
  delete     fbetaCmd;
  delete old_brbiasCmd;               // To be removed in G4 11.0
  delete     brbiasCmd;
  delete old_splitnucleiCmd;          // To be removed in G4 11.0
  delete     splitnucleiCmd;
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
  delete old_hlthCmd;                 // To be removed in G4 11.0
  delete     hlthCmd;
  delete old_userDecayDataCmd;        // To be removed in G4 11.0
  delete     userDecayDataCmd;
  delete old_userEvaporationDataCmd;  // To be removed in G4 11.0
  delete     userEvaporationDataCmd;
  delete old_colldirCmd;              // To be removed in G4 11.0
  delete     colldirCmd;
  delete old_collangleCmd;            // To be removed in G4 11.0
  delete     collangleCmd;
}
////////////////////////////////////////////////////////////////////////////////
//
void G4RadioactiveDecaymessenger::SetNewValue (G4UIcommand *command, G4String newValues)
{
  // Old commands to be removed in G4 11.0
  if ( command == old_nucleuslimitsCmd ) {
    theRadioactiveDecayContainer->
      SetNucleusLimits( old_nucleuslimitsCmd->GetNewNucleusLimitsValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/nucleusLimits in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_771", JustWarning, ed );    
  } else if ( command == old_analoguemcCmd ) {
    theRadioactiveDecayContainer->
      SetAnalogueMonteCarlo( old_analoguemcCmd->GetNewBoolValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/analogueMC in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_772", JustWarning, ed );        
  } else if ( command == old_fbetaCmd ) {
    theRadioactiveDecayContainer->SetFBeta( old_fbetaCmd->GetNewBoolValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/fBeta in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_773", JustWarning, ed );        
  } else if ( command == old_avolumeCmd ) {
    theRadioactiveDecayContainer->SelectAVolume( newValues );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/selectVolume in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_774", JustWarning, ed );        
  } else if ( command == old_deavolumeCmd ) {
    theRadioactiveDecayContainer->DeselectAVolume( newValues );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/deselectVolume in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_775", JustWarning, ed );        
  } else if ( command == old_allvolumesCmd ) {
    theRadioactiveDecayContainer->SelectAllVolumes();
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/allVolumes in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_776", JustWarning, ed );        
  } else if ( command == old_deallvolumesCmd ) {
    theRadioactiveDecayContainer->DeselectAllVolumes();
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/noVolumes in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_777", JustWarning, ed );        
  } else if ( command == old_brbiasCmd ) {
    theRadioactiveDecayContainer->SetBRBias( old_brbiasCmd->GetNewBoolValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/BRbias in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_778", JustWarning, ed );        
  } else if ( command == old_sourcetimeprofileCmd ) {
    theRadioactiveDecayContainer->SetSourceTimeProfile( newValues );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/sourceTimeProfile in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_779", JustWarning, ed );        
  } else if ( command == old_decaybiasprofileCmd ) {
    theRadioactiveDecayContainer->SetDecayBias( newValues );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/decayBiasProfile in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_780", JustWarning, ed );        
  } else if ( command == old_splitnucleiCmd ) {
    theRadioactiveDecayContainer->SetSplitNuclei( old_splitnucleiCmd->GetNewIntValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/splitNuclei in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_781", JustWarning, ed );        
  } else if ( command == old_verboseCmd ) {
    theRadioactiveDecayContainer->SetVerboseLevel( old_verboseCmd->GetNewIntValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/verbose in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_782", JustWarning, ed );        
  } else if ( command == old_icmCmd ) {
    theRadioactiveDecayContainer->SetICM( old_icmCmd->GetNewBoolValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/applyICM in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_783", JustWarning, ed );        
  } else if ( command == old_armCmd ) {
    theRadioactiveDecayContainer->SetARM( old_armCmd->GetNewBoolValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/applyARM in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_784", JustWarning, ed );        
  } else if ( command == old_hlthCmd ) {
    theRadioactiveDecayContainer->SetHLThreshold( old_hlthCmd->GetNewDoubleValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/hlThreshold in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_785", JustWarning, ed );        
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
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_786", JustWarning, ed );        
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
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_787", JustWarning, ed );        
  } else if ( command == old_colldirCmd ) {
    theRadioactiveDecayContainer->SetDecayDirection( old_colldirCmd->GetNew3VectorValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/decayDirection in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_788", JustWarning, ed );        
  } else if ( command == old_collangleCmd ) {
    theRadioactiveDecayContainer->SetDecayHalfAngle( old_collangleCmd->GetNewDoubleValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/decayHalfAngle in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactiveDecaymessenger", "HAD_RDM_789", JustWarning, ed );        
  }

  // New commands
  if ( command == nucleuslimitsCmd ) {
    theRadioactiveDecayContainer->
      SetNucleusLimits( nucleuslimitsCmd->GetNewNucleusLimitsValue( newValues ) );
  } else if ( command == analoguemcCmd ) {
    theRadioactiveDecayContainer->
      SetAnalogueMonteCarlo( analoguemcCmd->GetNewBoolValue( newValues ) );
  } else if ( command == fbetaCmd ) {
    theRadioactiveDecayContainer->SetFBeta( fbetaCmd->GetNewBoolValue( newValues ) );
  } else if ( command == avolumeCmd ) {
    theRadioactiveDecayContainer->SelectAVolume( newValues );
  } else if ( command == deavolumeCmd ) {
    theRadioactiveDecayContainer->DeselectAVolume( newValues );
  } else if ( command == allvolumesCmd ) {
    theRadioactiveDecayContainer->SelectAllVolumes();
  } else if ( command == deallvolumesCmd ) {
    theRadioactiveDecayContainer->DeselectAllVolumes();
  } else if ( command == brbiasCmd ) {
    theRadioactiveDecayContainer->SetBRBias( brbiasCmd->GetNewBoolValue( newValues ) );
  } else if ( command == sourcetimeprofileCmd ) {
    theRadioactiveDecayContainer->SetSourceTimeProfile( newValues );
  } else if ( command == decaybiasprofileCmd ) {
    theRadioactiveDecayContainer->SetDecayBias( newValues );
  } else if ( command == splitnucleiCmd ) {
    theRadioactiveDecayContainer->SetSplitNuclei( splitnucleiCmd->GetNewIntValue( newValues ) );
  } else if ( command == verboseCmd ) {
    theRadioactiveDecayContainer->SetVerboseLevel( verboseCmd->GetNewIntValue( newValues ) );
  } else if ( command == icmCmd ) {
    theRadioactiveDecayContainer->SetICM( icmCmd->GetNewBoolValue( newValues ) );
  } else if ( command == armCmd ) {
    theRadioactiveDecayContainer->SetARM( armCmd->GetNewBoolValue( newValues ) );
  } else if ( command == hlthCmd ) {
    theRadioactiveDecayContainer->SetHLThreshold( hlthCmd->GetNewDoubleValue( newValues ) );
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






