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
// $Id: DetectorMessenger.cc,v 1.2 2008/12/18 12:57:14 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
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
}

// ----------------------------------------------------------------------------

DetectorMessenger::~DetectorMessenger()
{
  delete theReadCommand;
  delete theDetectorDir;
}

// ----------------------------------------------------------------------------

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if ( command == theReadCommand )
  { 
    theDetector->SetReadFile(newValue );
  }
}
