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
#include "G4ORNLTrunk.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"

#include "G4VisAttributes.hh"

G4ORNLTrunk::G4ORNLTrunk()
{
}

G4ORNLTrunk::~G4ORNLTrunk()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4ORNLTrunk::ConstructTrunk(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/ORNLTrunk.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicTrunk = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("TrunkVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("TrunkPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("TrunkRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physTrunk = new G4PVPlacement(rm,position,
      			       "physicalTrunk",
  			       logicTrunk,
			       mother,
			       false,
			       0);

  // Visualization Attributes
  G4VisAttributes* TrunkVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  TrunkVisAtt->SetForceSolid(true);
  logicTrunk->SetVisAttributes(TrunkVisAtt);

  G4cout << "Trunk created !!!!!!" << G4endl;
  
  return physTrunk;
}
