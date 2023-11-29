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
 fWater(nullptr), fSoft(nullptr),  fSkeleton(nullptr), fLung(nullptr), fAdipose(nullptr), fGlandular(nullptr),
  fAdipose_glandular(nullptr)
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
  auto* elH = new G4Element ("Hydrogen","H",Z = 1.,A);

  A = 12.011*g/mole;
  auto* elC = new G4Element("Carbon","C",Z = 6.,A);  

  A = 14.01*g/mole;
  auto* elN = new G4Element("Nitrogen","N",Z = 7.,A);

  A = 16.00*g/mole;
  auto* elO = new G4Element("Oxygen","O",Z = 8.,A);

  A = 22.99*g/mole;
  auto* elNa = new G4Element("Sodium","Na",Z = 11.,A);

  A = 24.305*g/mole;
  auto* elMg = new G4Element("Magnesium","Mg",Z = 12.,A);

  A = 30.974*g/mole;
  auto* elP = new G4Element("Phosphorus","P",Z = 15.,A);
 
  A = 32.064*g/mole;
  auto* elS = new G4Element("Sulfur","S",Z = 16.,A);
 
  A = 35.453*g/mole;
  auto* elCl = new G4Element("Chlorine","Cl",Z = 17.,A);
 
  A = 39.098*g/mole;
  auto* elK = new G4Element("Potassium","K",Z = 19.,A);

  A = 40.08*g/mole;
  auto* elCa = new G4Element("Calcium","Ca",Z = 20.,A);

  A = 55.85*g/mole;
  auto* elFe  = new G4Element("Iron","Fe",Z = 26.,A);
 
  A = 65.38*g/mole;
  auto* elZn = new G4Element("Zinc","Zn",Z = 30.,A);

  A = 85.47 *g/mole;
  auto* elRb = new G4Element("Rb","Rb",Z = 37.,A);

  A = 87.62 *g/mole;
  auto* elSr = new G4Element("Sr","Sr",Z = 38.,A);

  A = 91.22 *g/mole;
  auto* elZr = new G4Element("Zr","Zr",Z = 40.,A);

  A = 207.19 *g/mole;
  auto* elPb = new G4Element("Lead","Pb", Z = 82.,A);

  // Water
  d = 1.000*g/cm3;
  fWater = new G4Material("Water",d,2);
  fWater->AddElement(elH,2);
  fWater->AddElement(elO,1);
  fWater->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

  // MIRD soft tissue
  d = 0.9869 *g/cm3;
  fSoft = new G4Material("soft_tissue",d,16);
  fSoft->AddElement(elH,0.1047);
  fSoft->AddElement(elC,0.2302);
  fSoft->AddElement(elN,0.0234);
  fSoft->AddElement(elO,0.6321);
  fSoft->AddElement(elNa,0.0013);
  fSoft->AddElement(elMg,0.00015);
  fSoft->AddElement(elP,0.0024);
  fSoft->AddElement(elS,0.0022);
  fSoft->AddElement(elCl,0.0014);
  fSoft->AddElement(elK,0.0021);
  fSoft->AddElement(elFe,0.000063);
  fSoft->AddElement(elZn,0.000032);
  fSoft->AddElement(elRb,0.0000057);
  fSoft->AddElement(elSr,0.00000034);
  fSoft->AddElement(elZr,0.000008);
  fSoft->AddElement(elPb,0.00000016);
 
  // MIRD Skeleton

  d = 1.4862*g/cm3;
  fSkeleton = new G4Material("skeleton",d,15);
  fSkeleton -> AddElement(elH,0.0704);
  fSkeleton -> AddElement(elC,0.2279);
  fSkeleton -> AddElement(elN,0.0387);
  fSkeleton -> AddElement(elO,0.4856);
  fSkeleton -> AddElement(elNa,0.0032); 
  fSkeleton -> AddElement(elMg,0.0011); 
  fSkeleton -> AddElement(elP,0.0694);
  fSkeleton -> AddElement(elS,0.0017);
  fSkeleton -> AddElement(elCl,0.0014);
  fSkeleton -> AddElement(elK,0.0015);
  fSkeleton -> AddElement(elCa,0.0991);
  fSkeleton -> AddElement(elFe,0.00008);
  fSkeleton -> AddElement(elZn,0.000048);
  fSkeleton -> AddElement(elSr,0.000032);
  fSkeleton -> AddElement(elPb,0.000011);
 
  // MIRD lung material
  d = 0.2958 *g/cm3;
  fLung = new G4Material("lung_material", d,16);
  fLung -> AddElement(elH, 0.1021);
  fLung -> AddElement(elC, 0.1001);
  fLung -> AddElement(elN,0.028);
  fLung -> AddElement(elO,0.7596);
  fLung -> AddElement(elNa,0.0019);
  fLung -> AddElement(elMg,0.000074);
  fLung -> AddElement(elP,0.00081);
  fLung -> AddElement(elS,0.0023);
  fLung -> AddElement(elCl,0.0027);
  fLung -> AddElement(elK,0.0020);
  fLung -> AddElement(elCa,0.00007);
  fLung -> AddElement(elFe,0.00037);
  fLung -> AddElement(elZn,0.000011);
  fLung -> AddElement(elRb,0.0000037);
  fLung -> AddElement(elSr,0.000000059);
  fLung -> AddElement(elPb,0.00000041);

  G4double density_adipose = 0.93 *g/cm3;
  fAdipose = new G4Material("adipose", density_adipose,8);  
  fAdipose -> AddElement(elH, 0.112);
  fAdipose -> AddElement(elC, 0.619);
  fAdipose -> AddElement(elN, 0.017);
  fAdipose -> AddElement(elO, 0.251);
  fAdipose -> AddElement(elS, 0.00025);
  fAdipose -> AddElement(elP, 0.00025);
  fAdipose -> AddElement(elK, 0.00025);
  fAdipose -> AddElement(elCa,0.00025);

  G4double density_glandular = 1.04 * g/cm3;
  fGlandular = new G4Material("glandular", density_glandular,8);
  fGlandular -> AddElement(elH, 0.1);
  fGlandular -> AddElement(elC,0.184);
  fGlandular -> AddElement(elN, 0.032);
  fGlandular -> AddElement(elO, 0.679);
  fGlandular -> AddElement(elS, 0.00125);
  fGlandular -> AddElement(elP, 0.00125);
  fGlandular -> AddElement(elK, 0.00125);
  fGlandular -> AddElement(elCa,0.00125);


  d = (density_adipose * 0.5) + (density_glandular * 0.5);
  fAdipose_glandular = new G4Material("adipose_glandular", d, 2);
  fAdipose_glandular -> AddMaterial(fAdipose, 0.5);
  fAdipose_glandular -> AddMaterial(fGlandular, 0.5);

  // Air 
  d = 1.290*mg/cm3;
  auto* matAir = new G4Material("Air",d,2);
  matAir->AddElement(elN,0.7);
  matAir->AddElement(elO,0.3);
}

G4Material* G4HumanPhantomMaterial::GetMaterial(G4String material)
{
  // Returns a material 
  G4Material* pttoMaterial = G4Material::GetMaterial(material); 
  return pttoMaterial; 
}
