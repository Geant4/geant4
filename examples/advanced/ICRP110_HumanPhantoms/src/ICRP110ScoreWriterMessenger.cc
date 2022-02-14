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
// Code developed by:
// S.Guatelli, M. Large and A. Malaroda, University of Wollongong
//
//
#include "ICRP110ScoreWriterMessenger.hh"
#include "ICRP110UserScoreWriter.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"
#include "G4RunManager.hh"

ICRP110ScoreWriterMessenger::ICRP110ScoreWriterMessenger(ICRP110UserScoreWriter* myUsrScWriter)
  :fUserScoreWriter(myUsrScWriter)
{ 
  fPhantomDir = new G4UIdirectory("/phantom/");
  fPhantomDir -> SetGuidance("Set Your Phantom.");
  
  fSexCmd = new G4UIcmdWithAString("/phantom/setScoreWriterSex",this);
  fSexCmd -> SetGuidance("Set sex of Phantom: Male or Female.");
  fSexCmd -> SetParameterName("ScoreWriterSex",true);
  fSexCmd -> SetDefaultValue("female");
  fSexCmd -> SetCandidates("male female");
  fSexCmd -> AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  fSectionCmd = new G4UIcmdWithAString("/phantom/setScoreWriterSection",this);
  fSectionCmd -> SetGuidance("Set section of Phantom: Head, Trunk or Full.");
  fSectionCmd -> SetParameterName("ScoreWriterSection",true);
  fSectionCmd -> SetDefaultValue("head");
  fSectionCmd -> SetCandidates("head trunk full");
  fSectionCmd -> AvailableForStates(G4State_PreInit,G4State_Idle); 
}

ICRP110ScoreWriterMessenger::~ICRP110ScoreWriterMessenger()
{
  delete  fSexCmd;
  delete  fSectionCmd;
  delete  fPhantomDir;
}

void ICRP110ScoreWriterMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{     
  if( command == fSexCmd )
    { 
     fUserScoreWriter -> SetPhantomSex(newValue); 
    }
  if( command == fSectionCmd )
    {
     fUserScoreWriter -> SetPhantomSection(newValue);
    }
}



