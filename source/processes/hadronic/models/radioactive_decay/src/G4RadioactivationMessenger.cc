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
  rdmDirectory = new G4UIdirectory("/process/had/rdm/");
  rdmDirectory->SetGuidance("Controls the biased version of radioactive decay");

  // Command to turn on/off variance reduction options
  analoguemcCmd = new G4UIcmdWithABool("/process/had/rdm/analogueMC",this);
  analoguemcCmd->SetGuidance("false: variance reduction method; true: analogue method");
  analoguemcCmd->SetParameterName("AnalogueMC",true);
  analoguemcCmd->SetDefaultValue(true);
  
  // Command to use branching ratio biasing or not
  brbiasCmd = new G4UIcmdWithABool("/process/had/rdm/BRbias",this);
  brbiasCmd->SetGuidance("false: no biasing; true: all branches are treated as equal");
  brbiasCmd->SetParameterName("BRBias",true);
  brbiasCmd->SetDefaultValue(true);
  
  // Command to set the half-life thresold for isomer production
  hlthCmd = new G4UIcmdWithADoubleAndUnit("/process/had/rdm/hlThreshold",this);
  hlthCmd->SetGuidance("Set the h-l threshold for isomer production");
  hlthCmd->SetParameterName("hlThreshold",false);
  hlthCmd->SetUnitCategory("Time");
  
  // Command to define the incident particle source time profile
  sourcetimeprofileCmd = new G4UIcmdWithAString("/process/had/rdm/sourceTimeProfile",this);
  sourcetimeprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the source particle time profile");
  sourcetimeprofileCmd->SetParameterName("STimeProfile",true);
  sourcetimeprofileCmd->SetDefaultValue("source.data");
  
  // Command to define the incident particle source time profile
  decaybiasprofileCmd = new G4UIcmdWithAString("/process/had/rdm/decayBiasProfile",this);
  decaybiasprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the decay bias time profile");
  decaybiasprofileCmd->SetParameterName("DBiasProfile",true);
  decaybiasprofileCmd->SetDefaultValue("bias.data");
  
  // Command to set nuclei splitting parameter
  splitnucleiCmd = new G4UIcmdWithAnInteger("/process/had/rdm/splitNuclei",this);
  splitnucleiCmd->SetGuidance("Set number of splitting for the isotopes.");
  splitnucleiCmd->SetParameterName("NSplit",true);
  splitnucleiCmd->SetDefaultValue(1);
  splitnucleiCmd->SetRange("NSplit>=1");
}


G4RadioactivationMessenger::~G4RadioactivationMessenger()
{
  delete rdmDirectory;
  delete analoguemcCmd;
  delete sourcetimeprofileCmd;
  delete decaybiasprofileCmd;
  delete brbiasCmd;
  delete splitnucleiCmd;
  delete hlthCmd;
}


void G4RadioactivationMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
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

