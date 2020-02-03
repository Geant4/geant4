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
// Previous authors: G. Guerrieri, S. Guatelli, and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
//

#include "G4HumanPhantomMaterial.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialTable.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

G4HumanPhantomMaterial::G4HumanPhantomMaterial(): 
  soft(0),  skeleton(0),lung(0), adipose(0), glandular(0),
  adipose_glandular(0)
{;}

G4HumanPhantomMaterial::~G4HumanPhantomMaterial()
{;}

void G4HumanPhantomMaterial::DefineMaterials()
{
  // Define required materials

  G4double A;  // atomic mass
  G4double Z;  // atomic number
  G4double d;  // density
 
  // General elements
 
  A = 1.01*g/mole;
  G4Element* elH = new G4Element ("Hydrogen","H",Z = 1.,A);

  A = 12.011*g/mole;
  G4Element* elC = new G4Element("Carbon","C",Z = 6.,A);  

  A = 14.01*g/mole;
  G4Element* elN = new G4Element("Nitrogen","N",Z = 7.,A);

  A = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A);

  A = 22.99*g/mole;
  G4Element* elNa = new G4Element("Sodium","Na",Z = 11.,A);

  A = 24.305*g/mole;
  G4Element* elMg = new G4Element("Magnesium","Mg",Z = 12.,A);

  A = 30.974*g/mole;
  G4Element* elP = new G4Element("Phosphorus","P",Z = 15.,A);
 
  A = 32.064*g/mole;
  G4Element* elS = new G4Element("Sulfur","S",Z = 16.,A);
 
  A = 35.453*g/mole;
  G4Element* elCl = new G4Element("Chlorine","Cl",Z = 17.,A);
 
  A = 39.098*g/mole;
  G4Element* elK = new G4Element("Potassium","K",Z = 19.,A);

  A = 40.08*g/mole;
  G4Element* elCa = new G4Element("Calcium","Ca",Z = 20.,A);

  A = 55.85*g/mole;
  G4Element* elFe  = new G4Element("Iron","Fe",Z = 26.,A);
 
  A = 65.38*g/mole;
  G4Element* elZn = new G4Element("Zinc","Zn",Z = 30.,A);

  A = 85.47 *g/mole;
  G4Element* elRb = new G4Element("Rb","Rb",Z = 37.,A);

  A = 87.62 *g/mole;
  G4Element* elSr = new G4Element("Sr","Sr",Z = 38.,A);

  A = 91.22 *g/mole;
  G4Element* elZr = new G4Element("Zr","Zr",Z = 40.,A);

  A = 207.19 *g/mole;
  G4Element* elPb = new G4Element("Lead","Pb", Z = 82.,A);

  // Water
  d = 1.000*g/cm3;
  matH2O = new G4Material("Water",d,2);
  matH2O->AddElement(elH,2);
  matH2O->AddElement(elO,1);
  matH2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

  // MIRD soft tissue
  d = 0.9869 *g/cm3;
  soft = new G4Material("soft_tissue",d,16);
  soft->AddElement(elH,0.1047);
  soft->AddElement(elC,0.2302);
  soft->AddElement(elN,0.0234);
  soft->AddElement(elO,0.6321);
  soft->AddElement(elNa,0.0013);
  soft->AddElement(elMg,0.00015);
  soft->AddElement(elP,0.0024);
  soft->AddElement(elS,0.0022);
  soft->AddElement(elCl,0.0014);
  soft->AddElement(elK,0.0021);
  soft->AddElement(elFe,0.000063);
  soft->AddElement(elZn,0.000032);
  soft->AddElement(elRb,0.0000057);
  soft->AddElement(elSr,0.00000034);
  soft->AddElement(elZr,0.000008);
  soft->AddElement(elPb,0.00000016);
 
  // MIRD Skeleton

  d = 1.4862*g/cm3;
  skeleton = new G4Material("skeleton",d,15);
  skeleton -> AddElement(elH,0.0704);
  skeleton -> AddElement(elC,0.2279);
  skeleton -> AddElement(elN,0.0387);
  skeleton -> AddElement(elO,0.4856);
  skeleton -> AddElement(elNa,0.0032); 
  skeleton -> AddElement(elMg,0.0011); 
  skeleton -> AddElement(elP,0.0694);
  skeleton -> AddElement(elS,0.0017);
  skeleton -> AddElement(elCl,0.0014);
  skeleton -> AddElement(elK,0.0015);
  skeleton -> AddElement(elCa,0.0991);
  skeleton -> AddElement(elFe,0.00008);
  skeleton -> AddElement(elZn,0.000048);
  skeleton -> AddElement(elSr,0.000032);
  skeleton -> AddElement(elPb,0.000011);
 
  // MIRD lung material
  d = 0.2958 *g/cm3;
  lung = new G4Material("lung_material", d,16);
  lung -> AddElement(elH, 0.1021);
  lung -> AddElement(elC, 0.1001);
  lung -> AddElement(elN,0.028);
  lung -> AddElement(elO,0.7596);
  lung -> AddElement(elNa,0.0019);
  lung -> AddElement(elMg,0.000074);
  lung -> AddElement(elP,0.00081);
  lung -> AddElement(elS,0.0023);
  lung -> AddElement(elCl,0.0027);
  lung -> AddElement(elK,0.0020);
  lung -> AddElement(elCa,0.00007);
  lung -> AddElement(elFe,0.00037);
  lung -> AddElement(elZn,0.000011);
  lung -> AddElement(elRb,0.0000037);
  lung -> AddElement(elSr,0.000000059);
  lung -> AddElement(elPb,0.00000041);

  G4double density_adipose = 0.93 *g/cm3;
  adipose = new G4Material("adipose", density_adipose,8);  
  adipose -> AddElement(elH, 0.112);
  adipose -> AddElement(elC, 0.619);
  adipose -> AddElement(elN, 0.017);
  adipose -> AddElement(elO, 0.251);
  adipose -> AddElement(elS, 0.00025);
  adipose -> AddElement(elP, 0.00025);
  adipose -> AddElement(elK, 0.00025);
  adipose -> AddElement(elCa,0.00025);

  G4double density_glandular = 1.04 * g/cm3;
  glandular = new G4Material("glandular", density_glandular,8);
  glandular -> AddElement(elH, 0.1);
  glandular -> AddElement(elC,0.184);
  glandular -> AddElement(elN, 0.032);
  glandular -> AddElement(elO, 0.679);
  glandular -> AddElement(elS, 0.00125);
  glandular -> AddElement(elP, 0.00125);
  glandular -> AddElement(elK, 0.00125);
  glandular -> AddElement(elCa,0.00125);


  d = (density_adipose * 0.5) + (density_glandular * 0.5);
  adipose_glandular = new G4Material("adipose_glandular", d, 2);
  adipose_glandular -> AddMaterial(adipose, 0.5);
  adipose_glandular -> AddMaterial(glandular, 0.5);

  // Air 
  d = 1.290*mg/cm3;
  G4Material* matAir = new G4Material("Air",d,2);
  matAir->AddElement(elN,0.7);
  matAir->AddElement(elO,0.3);
}

G4Material* G4HumanPhantomMaterial::GetMaterial(G4String material)
{
  // Returns a material 
  G4Material* pttoMaterial = G4Material::GetMaterial(material); 
  return pttoMaterial; 
}
