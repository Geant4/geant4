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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
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

   G4Exception("The trunk volume is missing !!!!!");

 body -> CreateOrgan("Uterus",trunkVolume, colourName, solidVis, sensitivity);  }

void G4CustomFemaleBuilder::BuildLeftOvary(const G4String& colourName, G4bool solidVis, G4bool sensitivity )

{
if (trunkVolume == 0)
    G4Exception("The trunk volume is missing !!!!!");

  body -> CreateOrgan("LeftOvary",trunkVolume, colourName,
 		     solidVis, sensitivity);  
}

void G4CustomFemaleBuilder::BuildRightOvary(const G4String& colourName, G4bool solidVis, G4bool sensitivity )

{
 if (trunkVolume == 0)

    G4Exception("The trunk volume is missing !!!!!");

  body -> CreateOrgan("RightOvary",trunkVolume, colourName,
 		     solidVis, sensitivity);  
}

void G4CustomFemaleBuilder::BuildVoxelLeftBreast(const G4String& colourName,
				      G4bool solidVis, G4bool sensitivity)
{ 
  G4cout << "BuildVoxelLeftBreast" << G4endl;
  if (motherVolume == 0)
    G4Exception("The world volume is missing !!!!!");
 
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
  G4Exception("The world volume is missing !!!!!");

  G4VBodyFactory* customBody = new G4VoxelBreastFactory();
  customBody -> CreateOrgan("RightBreast",motherVolume, colourName,
                       solidVis, sensitivity);  
  delete customBody;
}
