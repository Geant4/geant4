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
// $Id: ExN03DetectorConstruction.cc,v 1.2 2002-06-03 12:09:31 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//

#include "ExN03DetectorConstruction.hh"
#include "ExN03DetectorMessenger.hh"
#include "ExN03CalorimeterSD.hh"

//#include "G4UniformMagField.hh"
//#include "G4FieldManager.hh"
//#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
//#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "GDMLProcessor.hh"
#include "GDMLExpressionEvaluator.hh"

// Added here just to help resolve properly dependencies
#include "G4BooleanSolid.hh"

ExN03DetectorConstruction::ExN03DetectorConstruction()
: calorimeterSD( 0 ) {
  // create commands for interactive definition of the calorimeter  
  detectorMessenger = new ExN03DetectorMessenger(this);
  
  sxp.Initialize();
  config.SetURI( "NO3.gdml" );
  config.SetSetupName( "N03" );
  sxp.Configure( &config );
}

ExN03DetectorConstruction::~ExN03DetectorConstruction() {
  sxp.Finalize();
}

G4VPhysicalVolume* ExN03DetectorConstruction::Construct() { 
  sxp.Run();
  
  fWorld =  (G4VPhysicalVolume *)GDMLProcessor::GetInstance()->GetWorldVolume();
  
  if( fWorld == 0 ) {
    G4Exception(
        "World volume not set properly check your setup selection criteria or GDML input!"
               );
  }
  
  //                               
  // Sensitive Detectors: Absorber and Gap
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!calorimeterSD)
  {
    calorimeterSD = new ExN03CalorimeterSD("CalorSD",this);
    SDman->AddNewDetector( calorimeterSD );
  }
  G4LogicalVolume* lv = 0;
  lv = FindLogicalVolume( "Absorber" );
  if ( lv )
      lv->SetSensitiveDetector(calorimeterSD);
  lv = FindLogicalVolume( "Gap" );
  if ( lv )
      lv->SetSensitiveDetector(calorimeterSD);
  
  //                                        
  // Visualization attributes
  //
  lv = FindLogicalVolume( "World" );
  lv->SetVisAttributes (G4VisAttributes::Invisible);
  lv = FindLogicalVolume( "Calorimeter" );
  lv->SetVisAttributes (G4VisAttributes::Invisible);
  
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  lv = FindLogicalVolume( "Layer" );
  lv->SetVisAttributes(simpleBoxVisAtt);

  return fWorld;
}

G4LogicalVolume* ExN03DetectorConstruction::FindLogicalVolume( const G4String& vn ) {
  return const_cast<G4LogicalVolume*>( GDMLProcessor::GetInstance()->GetLogicalVolume( vn ) );
}

const G4VPhysicalVolume* ExN03DetectorConstruction::GetAbsorber() {
  return GDMLProcessor::GetInstance()->GetPhysicalVolume( "pv_Absorber_0" );
}

const G4VPhysicalVolume* ExN03DetectorConstruction::GetGap() {
  return GDMLProcessor::GetInstance()->GetPhysicalVolume( "pv_Gap_0" );
}

G4double ExN03DetectorConstruction::GetWorldSizeX() {
  GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
  return( calc->Eval( "WorldSizeX" ) );
}

G4double ExN03DetectorConstruction::GetWorldSizeYZ() {
  GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
  return( calc->Eval( "WorldSizeYZ" ) );
}

G4double ExN03DetectorConstruction::GetCalorThickness() {
  GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
  return( calc->Eval( "CalorThickness" ) );
}

G4double ExN03DetectorConstruction::GetCalorSizeYZ() {
  GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
  return( calc->Eval( "CalorSizeYZ" ) );
}

G4int ExN03DetectorConstruction::GetNbOfLayers() {
  GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
  return( (G4int)calc->Eval( "NbOfLayers" ) );
}

