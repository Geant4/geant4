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
//
#include "G4ORNLMaleBodyFactory.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomColour.hh"
#include "Saxana/SAXProcessor.h"
#include "Saxana/ProcessingConfigurator.h"
G4ORNLMaleBodyFactory::G4ORNLMaleBodyFactory()
{
}

G4ORNLMaleBodyFactory::~G4ORNLMaleBodyFactory()
{
}

G4VPhysicalVolume* G4ORNLMaleBodyFactory::CreateOrgan(const G4String& gdmlFile,
						      G4VPhysicalVolume* motherVolume,
						      const G4String& colourName, G4bool visAttribute,G4bool sensitivity)
{

 G4cout<< "ORNLBodyFactory: "<< "gdmlData/Male/ORNL"<< gdmlFile <<".gdml" << G4endl;
  SAXProcessor sxp;
  ProcessingConfigurator config;
  // Initialize GDML Processor
  sxp.Initialize();
  G4String name = "gdmlData/Male/ORNL"+ gdmlFile + ".gdml";
  config.SetURI(name);
  G4cout << "config.SetURI(name);" << G4endl;
  config.SetSetupName( "Default" );
  G4cout << "config.SetSetUpName;" << G4endl;
  sxp.Configure( &config );
  G4cout << "sxp.Configure( &config );" << G4endl;
  // Run GDML Processor
  sxp.Run(); 
  G4cout << "sxp.Run();" << G4endl;
 
  G4String logicalVolumeName = gdmlFile + "Volume"; 
  G4LogicalVolume* logicOrgan = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume(logicalVolumeName);

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("OrganPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("OrganRot");
 
  G4PhysicalVolumeStore::DeRegister((G4VPhysicalVolume*)GDMLProcessor::GetInstance()->GetWorldVolume());
 
  //   sxp.Finalize();
  // Define rotation and position here!
  G4VPhysicalVolume* physOrgan = new G4PVPlacement(rm,position,
      			       "physicalOrgan",
  			       logicOrgan,
			       motherVolume,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicOrgan->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }
// Visualization Attributes
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* organVisAtt = new G4VisAttributes(colour);
  // Visualization Attributes
  organVisAtt->SetForceSolid(visAttribute);
  logicOrgan->SetVisAttributes(organVisAtt);

  G4cout << "Organ created !!!!!!  from " << name <<G4endl;
  sxp.Finalize();

  return physOrgan;
}


