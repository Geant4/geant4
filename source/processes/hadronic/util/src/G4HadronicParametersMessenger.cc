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
//---------------------------------------------------------------------------
// ClassName: G4HadronicParametersMessenger
// Author: Alberto Ribon
// Date: May 2020
//---------------------------------------------------------------------------

#include "G4HadronicParametersMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4HadronicParameters.hh"


G4HadronicParametersMessenger::G4HadronicParametersMessenger( G4HadronicParameters* inputHadParam ) :
  theHadronicParameters( inputHadParam )
{
  // Main directory for the whole hadronic
  theDirectory = new G4UIdirectory( "/process/had/" );
  theDirectory->SetGuidance( "Control of general hadronic physics parameters and choices." );

  // This command setup the overall verbose level for hadronics.
  // It is the responsibility of all hadronic models, cross sections and processes to check
  // this verbose level and behave consistently: in particular, if the level is "0" there
  // should be no print out regardless of their eventual internal (local) verbosity level.
  theVerboseCmd = new G4UIcmdWithAnInteger( "/process/had/verbose", this );
  theVerboseCmd->SetGuidance( "Set verbose level: 0 (minimum), 1 (default), 2 (higher), ... (even higher)" );
  theVerboseCmd->SetParameterName( "VerboseLevel", true );
  theVerboseCmd->SetDefaultValue( 1 );
  theVerboseCmd->SetRange( "VerboseLevel>=0" );
  theVerboseCmd->AvailableForStates( G4State_PreInit );

  // This command sets the max energy for hadronics.
  theMaxEnergyCmd = new G4UIcmdWithADoubleAndUnit( "/process/had/maxEnergy", this );
  theMaxEnergyCmd->SetGuidance( "Max energy for hadronics (default: 100 TeV)" );
  theMaxEnergyCmd->SetParameterName( "MaxEnergy", false );
  theMaxEnergyCmd->SetUnitCategory( "Energy" );
  theMaxEnergyCmd->SetRange( "MaxEnergy>0.0" );
  theMaxEnergyCmd->AvailableForStates( G4State_PreInit );

  // This command enable the Cosmic Ray (CR) coalescence
  theCRCoalescenceCmd = new G4UIcmdWithABool( "/process/had/enableCRCoalescence", this );
  theCRCoalescenceCmd->SetGuidance( "Enable Cosmic Ray (CR) coalescence." );
  theCRCoalescenceCmd->SetParameterName( "EnableCRCoalescence", false );
  theCRCoalescenceCmd->SetDefaultValue( false );
}


G4HadronicParametersMessenger::~G4HadronicParametersMessenger() {
  delete theDirectory;
  delete theVerboseCmd;
  delete theMaxEnergyCmd;
  delete theCRCoalescenceCmd;
}


void G4HadronicParametersMessenger::SetNewValue( G4UIcommand *command, G4String newValues ) {
  if ( command == theVerboseCmd ) {
    // Note that the messenger updates only the verbose level of the class G4HadronicParameters:
    // it is the responsibility of each hadronic model, cross section and process to check
    // this verbose level and behave accordingly.
    theHadronicParameters->SetVerboseLevel( theVerboseCmd->GetNewIntValue( newValues ) );
  } else if ( command == theMaxEnergyCmd ) { 
    theHadronicParameters->SetMaxEnergy( theMaxEnergyCmd->GetNewDoubleValue( newValues ) );
  } else if ( command == theCRCoalescenceCmd ) {
    theHadronicParameters->SetEnableCRCoalescence( theCRCoalescenceCmd->GetNewBoolValue( newValues ) );
  }
}
