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
#include "G4MIRDPelvis.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

G4MIRDPelvis::G4MIRDPelvis()
{
}

G4MIRDPelvis::~G4MIRDPelvis()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4MIRDPelvis::ConstructPelvis(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/MIRDPelvis.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicPelvis = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("PelvisVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("PelvisPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("PelvisRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physPelvis = new G4PVPlacement(rm,position,
      			       "physicalPelvis",
  			       logicPelvis,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicPelvis->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* PelvisVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  PelvisVisAtt->SetForceSolid(true);
  logicPelvis->SetVisAttributes(PelvisVisAtt);

  G4cout << "Pelvis created !!!!!!" << G4endl;
  
  return physPelvis;
}
