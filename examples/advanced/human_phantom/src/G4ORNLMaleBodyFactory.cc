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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
// 
// Contribution: F. Ambroglini INFN Perugia, Italy
//
//
#ifdef G4LIB_USE_GDML

#include "G4ORNLMaleBodyFactory.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GDMLParser.hh"
#include "globals.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomColour.hh"

G4ORNLMaleBodyFactory::G4ORNLMaleBodyFactory()
{
}

G4ORNLMaleBodyFactory::~G4ORNLMaleBodyFactory()
{
}

G4VPhysicalVolume* G4ORNLMaleBodyFactory::CreateOrgan(const G4String& gdmlFile,
						      G4VPhysicalVolume* motherVolume,
						      const G4String& colourName, G4bool visAttribute, G4bool)
{

  G4cout<< "ORNLBodyFactory: "<< "gdmlData/Male/ORNL"<< gdmlFile <<".gdml" << G4endl;
  G4GDMLParser parser;
  G4String filename = "gdmlData/Male/ORNL"+ gdmlFile + ".gdml";
  parser.Read(filename);

  G4String logicalVolumeName = gdmlFile + "Volume"; 
  G4LogicalVolume* logicOrgan = parser.GetVolume(logicalVolumeName);
  G4ThreeVector position = parser.GetPosition("OrganPos");
  G4ThreeVector rot = parser.GetRotation("OrganRot");
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateX(rot.x()); rm->rotateY(rot.y()); rm->rotateZ(rot.z()); 

  G4PhysicalVolumeStore::DeRegister(parser.GetWorldVolume());
 
  // Define rotation and position here!
  G4VPhysicalVolume* physOrgan = new G4PVPlacement(rm,position,
      			       "physicalOrgan",
  			       logicOrgan,
			       motherVolume,
			       false,
			       0);

  
// Visualization Attributes
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* organVisAtt = new G4VisAttributes(colour);
  // Visualization Attributes
  organVisAtt->SetForceSolid(visAttribute);
  logicOrgan->SetVisAttributes(organVisAtt);

  G4cout << "Organ created !!!!!!  from " << filename <<G4endl;

  return physOrgan;
}
#endif
