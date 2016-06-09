//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: g4rwDetectorConstruction.cc,v 1.1 2004/12/06 11:01:14 radoone Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------

#include "globals.hh"
#include "G4VisAttributes.hh"

#include "g4rwDetectorConstruction.hh"

#include "Processor/GDMLProcessor.h"

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

