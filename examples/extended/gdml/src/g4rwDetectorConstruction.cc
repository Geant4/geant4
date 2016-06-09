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
// $Id: g4rwDetectorConstruction.cc,v 1.3 2006/06/29 17:21:47 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------

#include "globals.hh"
#include "G4VisAttributes.hh"

#include "g4rwDetectorConstruction.hh"

#include "G4Processor/GDMLProcessor.h"

// Added here just to help resolve properly dependencies
#include "G4BooleanSolid.hh"
#include "G4CSGSolid.hh"

gogdmlDetectorConstruction::gogdmlDetectorConstruction()
{
  sxp.Initialize();
  config.SetURI( "test.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );
}

gogdmlDetectorConstruction::~gogdmlDetectorConstruction()
{
  sxp.Finalize();
}

G4VPhysicalVolume* gogdmlDetectorConstruction::Construct()
{ 
  sxp.Run();
  
  fWorld =  (G4VPhysicalVolume *)GDMLProcessor::GetInstance()->GetWorldVolume();

  fWorld->GetLogicalVolume()->SetVisAttributes (G4VisAttributes::Invisible);

  
  if( fWorld == 0 ) {
    G4Exception(
      "World volume not set properly check your setup selection criteria or GDML input!"
      );
  }

  return fWorld;
}

