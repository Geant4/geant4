 
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
// $Id: HadrontherapyDetectorConstruction.hh,v 3.0, September 2004
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G. Candiano, G.A.P. Cirrone, F. Di Rosa, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "HadrontherapyMaterial.hh"
#include "HadrontherapyModulator.hh"
#include "HadrontherapyMaterial.hh"
#include "G4Transform3D.hh"
#include "G4ios.hh"
#include <fstream>
#include <strstream>

HadrontherapyModulator::HadrontherapyModulator()
  : numberOfModules(72), physiMotherMod(0)
{
  solidMod = new G4Tubs* [numberOfModules];
  for (G4int i=0; i< numberOfModules; i++)
    solidMod[i]=0;

  logicMod = new G4LogicalVolume* [numberOfModules];
  for (G4int i =0; i<numberOfModules; i++)
    logicMod[i]=0;
 
  physiMod = new G4VPhysicalVolume* [numberOfModules];
  for (G4int i=0; i<numberOfModules; i++)
    physiMod[i]=0; 
   
  startAngle = new G4DataVector;
  spanningAngle = new G4DataVector;
  Ztranslation = new G4DataVector;

  // rotation matrix of the mother volume
  rotationMatrix = new G4RotationMatrix(); 
  G4double angle = 90. *deg;             
  rotationMatrix -> rotateY(angle);
}

HadrontherapyModulator::~HadrontherapyModulator() 
{
  for (G4int i = 0; i < numberOfModules; i++)
    if (physiMod[i])
      delete physiMod[i];
 delete[] physiMod;

  for (G4int j = 0; j < numberOfModules; j++)
    if (logicMod[j])
      delete logicMod[j];
 delete[] logicMod;

  for (G4int k = 0; k<numberOfModules; k++)
    if (solidMod[k])
      delete solidMod[k];
 delete[] solidMod;

 delete Ztranslation;
 delete spanningAngle;
 delete startAngle;
 delete rotationMatrix;
}

void HadrontherapyModulator::BuildModulator(G4VPhysicalVolume* motherVolume)
{

  G4VPhysicalVolume* mother = motherVolume;
  ReadFile("modulator.dat");
 
  //Materials
  HadrontherapyMaterial* material = new HadrontherapyMaterial();
  G4Material* MotherModMater = material -> GetMat("Air");  
  G4Material* Mod0Mater = material -> GetMat("Air");
  G4Material* ModMater = material -> GetMat("PMMA");
  delete material;

 G4double innerRadiusOfTheTube = 2.5 *cm;
 G4double outerRadiusOfTheTube = 9.5 *cm;
 G4double hightOfTheTube = 0.03*cm;

  //----------------------------------------------------------
  //Mother volume  
  //----------------------------------------------------------
  
  G4double hightOfTheTube0 = 5.0 *cm;
  G4double startAngleOfTheTube0 = 45 *deg;
  G4double spanningAngleOfTheTube0 = 90 *deg;

 
 G4Box* solidMotherMod = new G4Box("MotherMod", 12. *cm, 12. *cm, 12. *cm);
 

 G4LogicalVolume* logicMotherMod = new G4LogicalVolume(solidMotherMod, MotherModMater,
                                      "MotherMod",0,0,0);
  
 G4ThreeVector positionMotherMod = G4ThreeVector(-2260.50 *mm, 30 *mm, 50 *mm);

 physiMotherMod = new G4PVPlacement(rotationMatrix,
                                    positionMotherMod,
                                     "MotherMod",        // its name
				   logicMotherMod,     // its logical volume				  
				   mother, // its mother  volume
				   false,           // no boolean operations
				   0);              // no particular field 
 

  G4ThreeVector positionMod0 = G4ThreeVector(0*cm,0*cm,0*cm);
 
 
  G4Tubs* solidMod0 = new G4Tubs("Mod0",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube0,
			 startAngleOfTheTube0, 
			 spanningAngleOfTheTube0);
  
  G4LogicalVolume* logicMod0 = new G4LogicalVolume(solidMod0, Mod0Mater, "Mod0",0,0,0);
  
  
  G4VPhysicalVolume* physiMod0 = new G4PVPlacement(0,positionMod0, 
                                                   "Mod0",        // its name 
				logicMod0,    // its logical volume			  
				physiMotherMod,  // its mother  volume
				false,         // no boolean operations
				0);            // no particular field
  

  G4RotationMatrix rm20;
  rm20.rotateZ(90 *deg);
  
  G4ThreeVector positionMod20 = G4ThreeVector(0*cm,0*cm,0*cm);
  
  G4Tubs* solidMod20 = new G4Tubs("Mod20",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube0,
			  startAngleOfTheTube0, 
			  spanningAngleOfTheTube0);
  
  G4LogicalVolume* logicMod20 = new G4LogicalVolume(solidMod20, Mod0Mater, "Mod0",0,0,0);
  
  
  G4VPhysicalVolume* physiMod20 = new G4PVPlacement(G4Transform3D(rm20, positionMod20), 
				 logicMod20,    // its logical volume			  
				 "Mod20",        // its name
				 logicMotherMod,  // its mother  volume
				 false,         // no boolean operations
				 0);            // no particular field
    
 
 G4RotationMatrix rm40;
 rm40.rotateZ(180 *deg);
  
 G4ThreeVector positionMod40 = G4ThreeVector(0*cm,0*cm,0*cm);
  
 G4Tubs*  solidMod40 = new G4Tubs("Mod40",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube0,
			  startAngleOfTheTube0, 
			  spanningAngleOfTheTube0);
  
  G4LogicalVolume* logicMod40 = new G4LogicalVolume(solidMod40, Mod0Mater, "Mod40",0,0,0);
  
  G4VPhysicalVolume* physiMod40 = new G4PVPlacement(G4Transform3D(rm40, positionMod40), 
				 logicMod40,    // its logical volume			  
				 "Mod40",        // its name
				 logicMotherMod,  // its mother  volume
				 false,         // no boolean operations
				 0);            // no particular field
  
  G4RotationMatrix rm60;
  rm60.rotateZ(270 *deg);
  
  G4ThreeVector positionMod60 = G4ThreeVector(0*cm,0*cm,0*cm);
  
  G4Tubs* solidMod60 = new G4Tubs("Mod60",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube0,
			  startAngleOfTheTube0, 
			  spanningAngleOfTheTube0);
  
  G4LogicalVolume* logicMod60 = new G4LogicalVolume(solidMod60, Mod0Mater, "Mod60",0,0,0);
  
  
  G4VPhysicalVolume* physiMod60 = new G4PVPlacement(G4Transform3D(rm60, positionMod60), 
				 logicMod60,    // its logical volume			  
				 "Mod60",        // its name
				 logicMotherMod,  // its mother  volume
				 false,         // no boolean operations
				 0);            // no particular field
  char name[6];

  G4VisAttributes * red = new G4VisAttributes( G4Colour(1. ,0. ,0.));
  red -> SetVisibility(true);
  red -> SetForceSolid(true);

  G4int numberSlices = 18;
 
  G4VPhysicalVolume* physiMother = 0;

 for (G4int i = 0; i < numberOfModules; i++)
 {
   if (i % numberSlices ==0 ) physiMother = physiMod0;
   if (i % numberSlices == 1) physiMother = physiMod20;
   if (i % numberSlices == 2) physiMother = physiMod40;
   if (i % numberSlices == 3) physiMother = physiMod60;

   sprintf(name, "Mod%2d", i);
   solidMod[i] = new G4Tubs(G4String(name), innerRadiusOfTheTube, outerRadiusOfTheTube, 
                            hightOfTheTube, (*startAngle)[i], (*spanningAngle)[i]);
   
   logicMod[i] = new G4LogicalVolume(solidMod[i], ModMater, G4String(name), 0,0,0);

   physiMod[i] = new G4PVPlacement(0,G4ThreeVector(0,0,(*Ztranslation)[i]),G4String(name),
                                    logicMod[i],
                                   physiMother, false,0);
				   
   logicMod[i] -> SetVisAttributes(red);
 }
}
void HadrontherapyModulator::SetModulatorAngle(G4double angle)
{
  G4double rotationAngle = angle;
  rotationMatrix -> rotateX(rotationAngle);
  physiMotherMod -> SetRotation(rotationMatrix);  
  G4cout << "MODULATOR HAS BEEN ROTATED OF   " << rotationAngle/deg << " deg" << G4endl;
}

void HadrontherapyModulator::ReadFile(G4String name)
{ 
  file = name;   
  ReadData(file);
  G4cout << file << " is the input file to model the modulator!" << G4endl;
}

void HadrontherapyModulator::ReadData(G4String fileName)
{
  char nameChar[100] = {""};
  std::ostrstream ost(nameChar, 100, std::ios::out);
 
  ost << fileName;
  
  G4String name(nameChar);
  
  std::ifstream file(fileName);
  std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
    {
	  G4String excep = "HadrontherapyModulator - data file: not found";
	  G4Exception(excep);
    }

  G4double a = 0;
  G4int k = 1;
  
  do
    {
      file >> a;
      // The file is organized into three columns:
      // 1st column is the start angle
      // 2nd column is the spanning angle
      // 3rd column is the translation along z
      // The file terminates with the pattern: -1   -1 -1
      
      if (a == -1)
	{
	  G4cout<<" The file is read!"<<G4endl;
	}
      else
	{
	  if (k == 1)
	    {	
	      G4double angle = a * deg;
	      startAngle -> push_back(angle);  
	      //G4cout<< angle/deg <<"deg";
	      k++;
	      
	    }
	  else if (k == 2)
	    {
	      G4double spanning = a * deg;
	      spanningAngle -> push_back(spanning);
	      //G4cout<<" "<< spanning/deg <<"spanning"<<G4endl;
	      k ++;
	    } 
	  else if (k == 3)
	    {
	      G4double translation = a * cm;
	      Ztranslation -> push_back(translation);
	      //G4cout<<" "<< translation/cm <<": x translation"<<G4endl;
	      k = 1;
	    }
	} 
    } while (a != -1); // end of file
  
  file.close();
}








