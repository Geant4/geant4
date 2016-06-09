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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#include "G4ORNLUpperLargeIntestine.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"

#include "G4VisAttributes.hh"

G4ORNLUpperLargeIntestine::G4ORNLUpperLargeIntestine()
{
}

G4ORNLUpperLargeIntestine::~G4ORNLUpperLargeIntestine()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4ORNLUpperLargeIntestine::ConstructUpperLargeIntestine(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/ORNLUpperLargeIntestine.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 
  G4LogicalVolume* logicUpperLargeIntestine = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("UpperLargeIntestineVolume");

  G4PhysicalVolumeStore::DeRegister((G4VPhysicalVolume*)GDMLProcessor::GetInstance()->GetWorldVolume());

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("UpperLargeIntestinePos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("UpperLargeIntestineRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physUpperLargeIntestine = new G4PVPlacement(rm,position,
      			       "physicalUpperLargeIntestine",
  			       logicUpperLargeIntestine,
			       mother,
			       false,
			       0);


  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicUpperLargeIntestine->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* UpperLargeIntestineVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  UpperLargeIntestineVisAtt->SetForceSolid(true);
  logicUpperLargeIntestine->SetVisAttributes(UpperLargeIntestineVisAtt);

  G4cout << "UpperLargeIntestine created !!!!!!" << G4endl;
  
  return physUpperLargeIntestine;
}
