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
//
// $Id: DetectorMessenger.cc,v 1.3 2009-04-15 13:26:26 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class DetectorMessenger implementation
//
// ----------------------------------------------------------------------------

#include "globals.hh"

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

// ----------------------------------------------------------------------------

DetectorMessenger::DetectorMessenger( DetectorConstruction* myDet )
  : theDetector( myDet )
{ 
  theDetectorDir = new G4UIdirectory( "/mydet/" );
  theDetectorDir->SetGuidance("Detector control.");

  theReadCommand = new G4UIcmdWithAString("/mydet/readFile", this);
  theReadCommand ->SetGuidance("READ GDML file with given name");
  theReadCommand ->SetParameterName("FileRead", false);
  theReadCommand ->SetDefaultValue("color_extension.gdml");
  theReadCommand ->AvailableForStates(G4State_PreInit);

  theWriteCommand = new G4UIcmdWithAString("/mydet/writeFile", this);
  theWriteCommand ->SetGuidance("WRITE GDML file with given name");
  theWriteCommand ->SetParameterName("FileWrite", false);
  theWriteCommand ->SetDefaultValue("color_extension_test.gdml");
  theWriteCommand ->AvailableForStates(G4State_PreInit);
}

// ----------------------------------------------------------------------------

DetectorMessenger::~DetectorMessenger()
{
  delete theReadCommand;
  delete theWriteCommand;
  delete theDetectorDir;
}

// ----------------------------------------------------------------------------

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if ( command == theReadCommand )
  { 
    theDetector->SetReadFile(newValue);
  }
  if ( command == theWriteCommand )
  { 
    theDetector->SetWriteFile(newValue);
  }
}
