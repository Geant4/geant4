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
#include "G4MIRDHead.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"

#include "G4VisAttributes.hh"

G4MIRDHead::G4MIRDHead():physHead(0)
{
 // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/MIRDHead.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

}

G4MIRDHead::~G4MIRDHead()
{
  sxp.Finalize();
}

void G4MIRDHead::ConstructHead(G4VPhysicalVolume* mother, G4String sex)
{

  //  if (sex == "female") 
  //  else 

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicHead = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("HeadVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("HeadPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("HeadRot");
  
  // Define rotation and position here!
  physHead = new G4PVPlacement(rm,position,
      			       "physicalHead",
  			       logicHead,
			       mother,
			       false,
			       0);

  // Visualization Attributes
  G4VisAttributes* HeadVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  HeadVisAtt->SetForceSolid(true);
  logicHead->SetVisAttributes(HeadVisAtt);

  G4cout << "Head created !!!!!!" << G4endl;
}
