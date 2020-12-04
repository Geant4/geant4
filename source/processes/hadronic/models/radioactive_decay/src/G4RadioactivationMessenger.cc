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
#include "G4HadronicException.hh"


G4RadioactivationMessenger::G4RadioactivationMessenger(G4Radioactivation* theRadioactivationContainer1)
 :theRadioactivationContainer(theRadioactivationContainer1)
{
  old_grdmDirectory = new G4UIdirectory("/grdm/");                                    // To be removed in G4 11.0
  old_grdmDirectory->SetGuidance("Controls the biased version of radioactive decay");

  rdmDirectory = new G4UIdirectory("/process/had/rdm/");
  rdmDirectory->SetGuidance("Controls the biased version of radioactive decay");

  // Command to turn on/off variance reduction options
  old_analoguemcCmd = new G4UIcmdWithABool("/grdm/analogueMC",this);                 // To be removed in G4 11.0
  old_analoguemcCmd->SetGuidance("false: variance reduction method; true: analogue method");
  old_analoguemcCmd->SetParameterName("OldAnalogueMC",true);
  old_analoguemcCmd->SetDefaultValue(true);

  analoguemcCmd = new G4UIcmdWithABool("/process/had/rdm/analogueMC",this);
  analoguemcCmd->SetGuidance("false: variance reduction method; true: analogue method");
  analoguemcCmd->SetParameterName("AnalogueMC",true);
  analoguemcCmd->SetDefaultValue(true);
  
  // Command to use branching ratio biasing or not
  old_brbiasCmd = new G4UIcmdWithABool("/grdm/BRbias",this);                         // To be removed in G4 11.0
  old_brbiasCmd->SetGuidance("false: no biasing; true: all branches are treated as equal");
  old_brbiasCmd->SetParameterName("OldBRBias",true);
  old_brbiasCmd->SetDefaultValue(true);

  brbiasCmd = new G4UIcmdWithABool("/process/had/rdm/BRbias",this);
  brbiasCmd->SetGuidance("false: no biasing; true: all branches are treated as equal");
  brbiasCmd->SetParameterName("BRBias",true);
  brbiasCmd->SetDefaultValue(true);
  
  // Command to set the half-life thresold for isomer production
  old_hlthCmd = new G4UIcmdWithADoubleAndUnit("/grdm/hlThreshold",this);              // To be removed in G4 11.0
  old_hlthCmd->SetGuidance("Set the h-l threshold for isomer production");
  old_hlthCmd->SetParameterName("OldhlThreshold",false);
  old_hlthCmd->SetUnitCategory("Time");

  hlthCmd = new G4UIcmdWithADoubleAndUnit("/process/had/rdm/hlThreshold",this);
  hlthCmd->SetGuidance("Set the h-l threshold for isomer production");
  hlthCmd->SetParameterName("hlThreshold",false);
  hlthCmd->SetUnitCategory("Time");
  
  // Command to define the incident particle source time profile
  old_sourcetimeprofileCmd = new G4UIcmdWithAString("/grdm/sourceTimeProfile",this);  // To be removed in G4 11.0
  old_sourcetimeprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the source particle time profile");
  old_sourcetimeprofileCmd->SetParameterName("OldSTimeProfile",true);
  old_sourcetimeprofileCmd->SetDefaultValue("source.data");

  sourcetimeprofileCmd = new G4UIcmdWithAString("/process/had/rdm/sourceTimeProfile",this);
  sourcetimeprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the source particle time profile");
  sourcetimeprofileCmd->SetParameterName("STimeProfile",true);
  sourcetimeprofileCmd->SetDefaultValue("source.data");
  
  // Command to define the incident particle source time profile
  old_decaybiasprofileCmd = new G4UIcmdWithAString("/grdm/decayBiasProfile",this);    // To be removed in G4 11.0
  old_decaybiasprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the decay bias time profile");
  old_decaybiasprofileCmd->SetParameterName("OldDBiasProfile",true);
  old_decaybiasprofileCmd->SetDefaultValue("bias.data");

  decaybiasprofileCmd = new G4UIcmdWithAString("/process/had/rdm/decayBiasProfile",this);
  decaybiasprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the decay bias time profile");
  decaybiasprofileCmd->SetParameterName("DBiasProfile",true);
  decaybiasprofileCmd->SetDefaultValue("bias.data");
  
  // Command to set nuclei splitting parameter
  old_splitnucleiCmd = new G4UIcmdWithAnInteger("/grdm/splitNuclei",this);            // To be removed in G4 11.0
  old_splitnucleiCmd->SetGuidance("Set number of splitting for the isotopes.");
  old_splitnucleiCmd->SetParameterName("OldNSplit",true);
  old_splitnucleiCmd->SetDefaultValue(1);
  old_splitnucleiCmd->SetRange("OldNSplit>=1");

  splitnucleiCmd = new G4UIcmdWithAnInteger("/process/had/rdm/splitNuclei",this);
  splitnucleiCmd->SetGuidance("Set number of splitting for the isotopes.");
  splitnucleiCmd->SetParameterName("NSplit",true);
  splitnucleiCmd->SetDefaultValue(1);
  splitnucleiCmd->SetRange("NSplit>=1");
}


G4RadioactivationMessenger::~G4RadioactivationMessenger()
{
  delete old_grdmDirectory;         // To be removed in G4 11.0
  delete rdmDirectory;
  delete old_analoguemcCmd;         // To be removed in G4 11.0
  delete analoguemcCmd;
  delete old_sourcetimeprofileCmd;  // To be removed in G4 11.0
  delete sourcetimeprofileCmd;
  delete old_decaybiasprofileCmd;   // To be removed in G4 11.0
  delete decaybiasprofileCmd;
  delete old_brbiasCmd;             // To be removed in G4 11.0
  delete brbiasCmd;
  delete old_splitnucleiCmd;        // To be removed in G4 11.0
  delete splitnucleiCmd;
  delete old_hlthCmd;               // To be removed in G4 11.0
  delete hlthCmd;
}


void G4RadioactivationMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  // Old commands to be removed in G4 11.0
  if ( command == old_analoguemcCmd ) { theRadioactivationContainer->
    SetAnalogueMonteCarlo( old_analoguemcCmd->GetNewBoolValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/analogueMC in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactivationMessenger", "HAD_RDM_871", JustWarning, ed );        
  } else if ( command == old_brbiasCmd ) { theRadioactivationContainer->
    SetBRBias( old_brbiasCmd->GetNewBoolValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/BRbias in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactivationMessenger", "HAD_RDM_872", JustWarning, ed );        
  } else if ( command == old_sourcetimeprofileCmd ) { theRadioactivationContainer->
    SetSourceTimeProfile( newValues );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/sourceTimeProfile in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactivationMessenger", "HAD_RDM_873", JustWarning, ed );        
  } else if ( command == old_decaybiasprofileCmd ) { theRadioactivationContainer->
    SetDecayBias( newValues );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/decayBiasProfile in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactivationMessenger", "HAD_RDM_874", JustWarning, ed );        
  } else if ( command == old_splitnucleiCmd ) { theRadioactivationContainer->
    SetSplitNuclei( old_splitnucleiCmd->GetNewIntValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/splitNuclei in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactivationMessenger", "HAD_RDM_875", JustWarning, ed );        
  } else if ( command == old_hlthCmd ) { theRadioactivationContainer->
    SetHLThreshold( old_hlthCmd->GetNewDoubleValue( newValues ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/rdm/hlThreshold in the next major release, Geant4 version 11.0";
    G4Exception( "G4RadioactivationMessenger", "HAD_RDM_876", JustWarning, ed );        
  }

  // New commands
  if ( command == analoguemcCmd ) { theRadioactivationContainer->
    SetAnalogueMonteCarlo( analoguemcCmd->GetNewBoolValue( newValues ) );
  } else if ( command == brbiasCmd ) { theRadioactivationContainer->
    SetBRBias( brbiasCmd->GetNewBoolValue( newValues ) );
  } else if ( command == sourcetimeprofileCmd ) { theRadioactivationContainer->
    SetSourceTimeProfile( newValues );
  } else if ( command == decaybiasprofileCmd ) { theRadioactivationContainer->
    SetDecayBias( newValues );
  } else if ( command == splitnucleiCmd ) { theRadioactivationContainer->
    SetSplitNuclei( splitnucleiCmd->GetNewIntValue( newValues ) );
  } else if ( command == hlthCmd ) { theRadioactivationContainer->
    SetHLThreshold( hlthCmd->GetNewDoubleValue( newValues ) );
  }
}

