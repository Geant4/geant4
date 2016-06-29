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
// Authors: S. Guatelli , M. G. Pia, INFN Genova and F. Ambroglini INFN Perugia, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//

#include"G4CustomFemaleBuilder.hh"

#include "G4VBodyFactory.hh"
#include "G4VoxelBreastFactory.hh"
//#include "G4MIRDBodyFactory.hh"
//#include "G4ORNLCustomFemaleBodyFactory.hh"


G4CustomFemaleBuilder::G4CustomFemaleBuilder()
{  
}

G4CustomFemaleBuilder::~G4CustomFemaleBuilder()

{ 
  delete body;
} 

void G4CustomFemaleBuilder::BuildUterus(const G4String& colourName, G4bool solidVis, G4bool sensitivity )

{
if (trunkVolume == 0)
  G4Exception("G4CustomFemaleBuilder::BuildUterus()", "human_phantom0001", FatalException, "The trunk volume is missing !!!!!");

 body -> CreateOrgan("Uterus",trunkVolume, colourName, solidVis, sensitivity);  }

void G4CustomFemaleBuilder::BuildLeftOvary(const G4String& colourName, G4bool solidVis, G4bool sensitivity )

{
  if (trunkVolume == 0)
    G4Exception("G4CustomFemaleBuilder::BuildLeftOvary()", "human_phantom0002", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("LeftOvary",trunkVolume, colourName,
 		     solidVis, sensitivity);  
}

void G4CustomFemaleBuilder::BuildRightOvary(const G4String& colourName, G4bool solidVis, G4bool sensitivity )

{
  if (trunkVolume == 0)
   G4Exception("G4CustomFemaleBuilder::BuildRightOvary()", "human_phantom0003", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("RightOvary",trunkVolume, colourName,
 		     solidVis, sensitivity);  
}

void G4CustomFemaleBuilder::BuildVoxelLeftBreast(const G4String& colourName,
				      G4bool solidVis, G4bool sensitivity)
{ 
  G4cout << "BuildVoxelLeftBreast" << G4endl;
  if (motherVolume == 0)
    G4Exception("G4CustomFemaleBuilder::BuildVoxelLeftBreast()", "human_phantom0004", FatalException, "The world volume is missing !!!!!");
 
  G4VBodyFactory* customBody = new G4VoxelBreastFactory();
  customBody -> CreateOrgan("LeftBreast",motherVolume, colourName,
		      solidVis, sensitivity);  
  delete customBody;
}
  

void G4CustomFemaleBuilder::BuildVoxelRightBreast(const G4String& colourName, 
				       G4bool solidVis, 
				       G4bool sensitivity)
{
  if (motherVolume == 0)
    G4Exception("G4CustomFemaleBuilder::BuildVoxelRightBreast()", "human_phantom0005", FatalException, "The world volume is missing !!!!!");

  G4VBodyFactory* customBody = new G4VoxelBreastFactory();
  customBody -> CreateOrgan("RightBreast",motherVolume, colourName,
                       solidVis, sensitivity);  
  delete customBody;
}
