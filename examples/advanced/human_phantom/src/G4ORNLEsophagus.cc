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
#include "G4ORNLEsophagus.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"

#include "G4VisAttributes.hh"

G4ORNLEsophagus::G4ORNLEsophagus()
{
}

G4ORNLEsophagus::~G4ORNLEsophagus()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4ORNLEsophagus::ConstructEsophagus(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/ORNLEsophagus.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicEsophagus = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("EsophagusVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("EsophagusPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("EsophagusRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physEsophagus = new G4PVPlacement(rm,position,
      			       "physicalEsophagus",
  			       logicEsophagus,
			       mother,
			       false,
			       0);

  // Visualization Attributes
  G4VisAttributes* EsophagusVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  EsophagusVisAtt->SetForceSolid(true);
  logicEsophagus->SetVisAttributes(EsophagusVisAtt);

  G4cout << "Esophagus created !!!!!!" << G4endl;
  
  return physEsophagus;
}
