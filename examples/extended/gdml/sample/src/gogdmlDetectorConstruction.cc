// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: gogdmlDetectorConstruction.cc,v 1.1.1.1 2002-05-31 00:34:43 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "gogdmlDetectorConstruction.hh"

#include "GDMLProcessor.hh"

// Added here just to help resolve properly dependencies
#include "G4BooleanSolid.hh"
#include "G4CSGSolid.hh"

gogdmlDetectorConstruction::gogdmlDetectorConstruction()
{
  sxp.Initialize();
  config.SetURI( "test.gdml" );
  config.SetSetupName( "Test1" );
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
  
  if( fWorld == 0 ) {
    G4Exception(
        "World volume not set properly check your setup selection criteria or GDML input!"
               );
  }

  return fWorld;
}

