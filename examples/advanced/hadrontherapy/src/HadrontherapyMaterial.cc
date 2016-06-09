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
// $Id: HadrontherapyMaterial.cc; Last modified: G.A.P.Cirrone, March 2008;
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------
#include "globals.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "HadrontherapyMaterial.hh"
#include "G4NistManager.hh"

HadrontherapyMaterial::HadrontherapyMaterial():
  matW(0), matAl(0), matSi(0), matTa(0), matCu(0),
  vacuum(0), matplexiglass(0), brass(0),
  kapton(0), matPb(0), titanium(0), matAir(0),
  matH2O(0), soft(0), gold(0), bone(0), muscle(0)
{;}

HadrontherapyMaterial::~HadrontherapyMaterial()
{
  delete muscle;
  delete bone;
  delete gold;
  delete soft;
  delete matH2O;
  delete matAir;
  delete titanium;
  delete matPb;
  delete kapton;
  delete brass;
  delete matplexiglass;
  delete vacuum;
  delete matCu;
  delete matTa;
  delete matSi;
  delete matAl;
  delete matW;
}

void HadrontherapyMaterial::DefineMaterials()
{
  // Define required materials

  G4double z; // Atomic numebr
  G4double a; // Atomic mass
  G4double d; // Density
  G4int nComponents;// Number of components 
  G4double fractionmass; // Fraction in mass of an element in a material
  G4int nAtoms; // Number of atoms in a molecuule

  /////////////////////////////////////////////////////////////////////////////
  // MATERIA DEFINITION FOLLOWING THE NIST DATABASE

  // Pointer to the G4Nist manager
  // for the material definition following
  // the NIST database.
  // It recommended to use this when 
  // the Standard models for electromagnetic physic
  // are called

  G4NistManager* nistMaterialManager = G4NistManager::Instance();
  G4bool isotopes = false;

  // Material NIST definition
  nistMaterialManager -> FindOrBuildMaterial("G4_AIR"  , isotopes);
  nistMaterialManager -> FindOrBuildMaterial("G4_WATER", isotopes);
  nistMaterialManager -> FindOrBuildMaterial("G4_PMMA", isotopes);
  nistMaterialManager -> FindOrBuildMaterial("G4_MYLAR", isotopes);
  /////////////////////////////////////////////////////////////////////////////

  // Elements 
  a = 1.01*g/mole;
  G4Element* elH = new G4Element ("Hydrogen","H",z = 1.,a);
  
  a = 14.01*g/mole;
  G4Element* elN = new G4Element("Nitrogen","N",z = 7.,a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O",z = 8.,a);

  a = 12.011*g/mole;
  G4Element* elC = new G4Element("Carbon","C",z = 6.,a);
 
  a = 22.99*g/mole;
  G4Element* elNa = new G4Element("Sodium","Na",z = 11.,a);
 
  a = 24.305*g/mole;
  G4Element* elMg = new G4Element("Magnesium","Mg",z = 12.,a);

  a = 30.974*g/mole;
  G4Element* elP = new G4Element("Phosphorus","P",z = 15.,a);
 
  a = 32.06*g/mole;
  G4Element* elS = new G4Element("Sulfur","S",z = 16.,a);
 
  a = 35.453*g/mole;
  G4Element* elCl = new G4Element("Chlorine","Cl",z = 17.,a);
 
  a = 39.098*g/mole;
  G4Element* elK = new G4Element("Potassium","K",z = 19.,a);

  a = 40.08*g/mole;
  G4Element* elCa = new G4Element("Calcium","Ca",z = 20.,a);
 
  a = 65.38*g/mole;
  G4Element* elZn = new G4Element("Zinc","Zn",z = 30.,a);
 
  a = 55.85*g/mole;
  G4Element* elFe = new G4Element("Iron","Fe",z = 26.,a);
  
  a = 63.546*g/mole;
  d = 8.90*g/cm3;
  G4Element* elCu = new G4Element("Copper","Cu", z = 29., a);  

  // Materials 
  
  // Tungsten
  z = 74.;
  a = 183.84* g/mole;
  d = 19.3*g/cm3;
  matW = new G4Material("Tungsten", z, a, d);
 
  // Aluminum
  z = 13.;
  a = 26.98*g/mole;
  d = 2.700*g/cm3;
  matAl = new G4Material("MatAluminum",z, a, d);
 
  // Silicon
  z = 14.;
  a = 28.09*g/mole;
  d = 2.330*g/cm3;
  matSi = new G4Material("MatSilicon",z, a, d);
 
  // Tantalum
  z = 73.;
  a = 180.948*g/mole;
  d = 16.6*g/cm3;
  matTa = new G4Material("MatTantalum", z, a, d);
 
  // Copper
  z = 29.;
  a = 63.546*g/mole;
  d = 8.90*g/cm3;
  matCu = new G4Material("MatCopper", z, a, d);
  
  // Vacuum
  G4double density = universe_mean_density;
  G4double pressure = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  a = 1.01*g/mole;
  z = 1;
  vacuum = new G4Material("Galactic", z, a,
			  density,kStateGas,temperature,pressure);
  
  // Plexiglass
  d = 1.18*g/cm3;
  matplexiglass = new G4Material("PMMA",d,3);
  matplexiglass -> AddElement(elH,0.08);
  matplexiglass -> AddElement(elC,0.60);
  matplexiglass -> AddElement(elO,0.32);
 
  // Brass
  d = 8.40*g/cm3;
  nComponents = 2;
  G4Material* brass = new G4Material("Brass", d, nComponents);        
  brass -> AddElement(elZn, fractionmass = 30*perCent);
  brass -> AddElement(elCu, fractionmass = 70*perCent);

  // Kapton
  d = 1.43*g/cm3;
  nComponents = 4;
  G4Material* kapton = new G4Material("Kapton", d, nComponents);
  kapton -> AddElement(elH, nAtoms = 10);
  kapton -> AddElement(elO, nAtoms = 5);
  kapton -> AddElement(elC, nAtoms = 22);
  kapton -> AddElement(elN, nAtoms = 2);

  // Lead 
  a = 207.19*g/mole;
  z = 82.;
  d = 11.35*g/cm3;
  matPb = new G4Material("Lead", z, a, d);

  // Titanium
  z = 22.;
  a = 47.88*g/mole;
  d = 4.50*g/cm3;
  titanium = new G4Material("titanium", z, a, d);
 
  // Air material
  d = 1.290*mg/cm3;
  nComponents = 2;
  G4Material* matAir = new G4Material("Air", d, nComponents);
  matAir -> AddElement(elN,0.7);
  matAir -> AddElement(elO,0.3);

  // Water by "hand"
  d = 1.000*g/cm3;
  nComponents = 2;
  matH2O = new G4Material("Water", d, nComponents);
  matH2O -> AddElement(elH,2);
  matH2O -> AddElement(elO,1);
  matH2O -> GetIonisation()->SetMeanExcitationEnergy(75.0*eV);
  matH2O -> SetChemicalFormula("H_2O");
  G4cout << "-----------> CHEMICAL FORMULA FOR WATER FIXED <----------"<< G4endl;
 
 
  //soft tissue(http://www.nist.gov)
  d = 1.0*g/cm3;
  nComponents = 13;
  soft = new G4Material("tissue",d, nComponents);
  soft -> AddElement(elH,0.104472);
  soft -> AddElement(elC,0.23219);
  soft -> AddElement(elN,0.02488);
  soft -> AddElement(elO,0.630238);
  soft -> AddElement(elNa,0.00113);
  soft -> AddElement(elMg,0.00013);
  soft -> AddElement(elP,0.00133);
  soft -> AddElement(elS,0.00199);
  soft -> AddElement(elCl,0.00134);
  soft -> AddElement(elK,0.00199);
  soft -> AddElement(elCa,0.00023);
  soft -> AddElement(elFe,0.00005);
  soft -> AddElement(elZn,0.00003); 
  
  // Gold
  z = 79;
  a = 196.97*g/mole;
  d = 19.32*g/cm3;
  gold = new G4Material("gold", z, a, d);

  // Compact bone 
  d = 1.85*g/cm3;
  nComponents = 8;
  bone = new G4Material("bone", d, nComponents);
  bone -> AddElement(elH,0.063984);
  bone -> AddElement(elC,0.278);
  bone -> AddElement(elN,0.027);
  bone -> AddElement(elO,0.410016);
  bone -> AddElement(elMg,0.002);
  bone -> AddElement(elP,0.07);
  bone -> AddElement(elS,0.002);
  bone -> AddElement(elCa,0.147);

  //muscle
  nComponents = 9;
  muscle = new G4Material("muscle", d, nComponents);
  muscle -> AddElement(elH,0.101997);
  muscle -> AddElement(elC,0.123);
  muscle -> AddElement(elN,0.035);
  muscle -> AddElement(elNa,0.0008);
  muscle -> AddElement(elO,0.729);
  muscle -> AddElement(elMg,0.0002);
  muscle -> AddElement(elP,0.002);
  muscle -> AddElement(elS,0.005);
  muscle -> AddElement(elK,0.003);
}

G4Material* HadrontherapyMaterial::GetMat(G4String material)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(material); 
  return pttoMaterial; 
}

