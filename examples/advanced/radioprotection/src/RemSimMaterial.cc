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
// Code developed by:
//  S.Guatelli
//
//    *******************************
//    *                             *
//    *    RemSimMaterial.cc        *
//    *                             *
//    *******************************
//
// $Id: RemSimMaterial.cc,v 1.3 2004-03-12 10:55:55 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "globals.hh"
#include "G4MaterialTable.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "RemSimMaterial.hh"

RemSimMaterial::RemSimMaterial():
  matPb(0),matAir(0),matH2O(0),vacuum(0) 
{;}

RemSimMaterial::~RemSimMaterial()
{
  delete vacuum;
  delete matH2O;
  delete matAir;
  delete matPb;
}

void RemSimMaterial::DefineMaterials()
{
  // Define required materials

  G4double A;  // atomic mass
  G4double Z;  // atomic number
  G4double d;  // density
 
  // General elements
 
  A = 1.01*g/mole;
  G4Element* elH = new G4Element ("Hydrogen","H",Z = 1.,A);
  
  A = 14.01*g/mole;
  G4Element* elN = new G4Element("Nitrogen","N",Z = 7.,A);

  A = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A);
  
  // Lead material
  A = 207.19*g/mole;
  Z = 82;
  d = 11.35*g/cm3;
  matPb = new G4Material("Lead",Z,A,d);
 
  // Air material
  G4double airDensity = 1.290*mg/cm3;
  G4Material* matAir = new G4Material("Air",airDensity,2);
  matAir->AddElement(elN,0.7);
  matAir->AddElement(elO,0.3);

  // Water
  d = 1.000*g/cm3;
  matH2O = new G4Material("Water",d,2);
  matH2O->AddElement(elH,2);
  matH2O->AddElement(elO,1);
  matH2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);
   
  //vacuum
  G4double vacuumDensity = 1.e-25 *g/cm3;
  G4double pressure = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  vacuum = new G4Material("Galactic", Z=1., A=1.01*g/mole,
			  vacuumDensity,kStateGas,temperature,pressure);
}

G4Material* RemSimMaterial::GetMaterial(G4String material)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(material); 
  return pttoMaterial; 
}
