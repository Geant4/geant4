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
// ClassName: G4FTFTuningsMessenger
// Author: Alberto Ribon
// Date: August 2022
//---------------------------------------------------------------------------

#include "G4FTFTuningsMessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4FTFTunings.hh"


G4FTFTuningsMessenger::G4FTFTuningsMessenger() {
  // These two commands select an alternative set of FTF parameters (called "tune"):
  // either via its integer index
  theFTFTuneIndexCmd = new G4UIcmdWithAnInteger( "/process/had/models/ftf/selectTuneByIndex", this );
  theFTFTuneIndexCmd->SetGuidance( "Select one FTF set of parameters (tune) via its index: 0 (default), 1, 2, ..." );
  theFTFTuneIndexCmd->SetParameterName( "indexFTFTune", true );
  theFTFTuneIndexCmd->SetDefaultValue( 0 );
  theFTFTuneIndexCmd->SetRange( "indexFTFTune>=0" );
  theFTFTuneIndexCmd->AvailableForStates( G4State_PreInit );
  // or via its string name
  theFTFTuneNameCmd = new G4UIcmdWithAString( "/process/had/models/ftf/selectTuneByName", this );
  theFTFTuneNameCmd->SetGuidance( "Select one FTF set of parametes (tune) via its name (string)." );
  theFTFTuneNameCmd->SetGuidance( " (default) is the default." );
  theFTFTuneNameCmd->SetParameterName( "nameFTFTune", true );
  theFTFTuneNameCmd->SetDefaultValue( "default" );
  theFTFTuneNameCmd->AvailableForStates( G4State_PreInit );
}


G4FTFTuningsMessenger::~G4FTFTuningsMessenger() {
  delete theFTFTuneIndexCmd;
  delete theFTFTuneNameCmd;
}


void G4FTFTuningsMessenger::SetNewValue( G4UIcommand *command, G4String newValues ) {
  if ( command == theFTFTuneIndexCmd  ||  command == theFTFTuneNameCmd ) {
    G4int index = -999;
    if ( command == theFTFTuneIndexCmd ) {
      G4int value = theFTFTuneIndexCmd->GetNewIntValue( newValues );
      if ( value >= 0  &&  value < G4FTFTunings::sNumberOfTunes ) {
        index = value;
      } else {
        G4ExceptionDescription ed;
        ed << "The FTF tune index=" << value << " value is wrong!";
        command->CommandFailed( ed );
      }
    } else {
      for ( G4int i = 0; i < G4FTFTunings::sNumberOfTunes; ++i ) {
        if ( newValues == G4FTFTunings::Instance()->GetTuneName(i) ) {
          index = i;
          break;
        }
      }
      if ( index < 0 ) {
        G4ExceptionDescription ed;
        ed << "The FTF tune name=" << newValues << " is not found!";
        command->CommandFailed( ed );
      }
    }
    if ( index >= 0 ) {
      // Tune applicability state: 0 means switched off; 1 means switched on.
      G4FTFTunings::Instance()->SetTuneApplicabilityState(index, 1);
    }
  }
}
