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
#include "G4VoxelBreastFactory.hh"
#include "G4VoxelLeftBreast.hh"
//#include "G4ParameterisedRightBreast.hh"

G4VoxelBreastFactory::G4VoxelBreastFactory()
{
  // Map with name of the organ and pointer to the MIRDOrgan class
  //  organ["ParameterisedRightBreast"] = new G4ParameterisedRightBreast();
  organ["VoxelLeftBreast"] = new G4VoxelLeftBreast();
}

G4VoxelBreastFactory::~G4VoxelBreastFactory()
{ 
  //delete organ["ParameterisedRightBreast"]; organ["ParameterisedRightBreast"]=0;
  delete organ["VoxelLeftBreast"]; 
  organ["VoxelLeftBreast"]=0;
}



G4VPhysicalVolume* G4VoxelBreastFactory::CreateOrgan(const G4String& organ_name,
						     G4VPhysicalVolume* motherVolume,
						     const G4String& colourName, 
						     G4bool visAttribute,
						     G4bool sensitivity)
{
 return organ[organ_name] -> Construct(organ_name,motherVolume,colourName, visAttribute, sensitivity);
}


