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
#include "G4ORNLBrain.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"

#include "G4VisAttributes.hh"

G4ORNLBrain::G4ORNLBrain()
{
}

G4ORNLBrain::~G4ORNLBrain()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4ORNLBrain::ConstructBrain(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/ORNLBrain.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicBrain = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("BrainVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("BrainPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("BrainRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physBrain = new G4PVPlacement(rm,position,
      			       "physicalBrain",
  			       logicBrain,
			       mother,
			       false,
			       0);

  // Visualization Attributes
  G4VisAttributes* BrainVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  BrainVisAtt->SetForceSolid(true);
  logicBrain->SetVisAttributes(BrainVisAtt);

  G4cout << "Brain created !!!!!!" << G4endl;
  
  return physBrain;
}
