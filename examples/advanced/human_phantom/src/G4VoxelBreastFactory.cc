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
#include "G4VoxelBreastFactory.hh"
#include "G4VoxelLeftBreast.hh"
#include "G4VoxelRightBreast.hh"

G4VoxelBreastFactory::G4VoxelBreastFactory()
{
  // Map with name of the organ and pointer to the MIRDOrgan class
  //  organ["ParameterisedRightBreast"] = new G4ParameterisedRightBreast();
  organ["LeftBreast"] = new G4VoxelLeftBreast();
  organ["RightBreast"] = new G4VoxelRightBreast();
}

G4VoxelBreastFactory::~G4VoxelBreastFactory()
{ 
  delete organ["RightBreast"]; 
  organ["RightBreast"]=0;

  delete organ["LeftBreast"]; 
  organ["LeftBreast"]=0;
}



G4VPhysicalVolume* G4VoxelBreastFactory::CreateOrgan(const G4String& organ_name,
						     G4VPhysicalVolume* motherVolume,
						     const G4String& colourName, 
						     G4bool visAttribute,
						     G4bool sensitivity)
{
 return organ[organ_name] -> Construct(organ_name,motherVolume,colourName, visAttribute, sensitivity);
}


