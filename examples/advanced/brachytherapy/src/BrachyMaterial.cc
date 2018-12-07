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
// Code developed by:
//  S.Guatelli, A. Le
//
//    *******************************
//    *                             *
//    *    BrachyMaterial.cc        *
//    *                             *
//    *******************************
//
//
#include "globals.hh"
#include "Randomize.hh"  
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "BrachyMaterial.hh"

BrachyMaterial::BrachyMaterial()
{;}

BrachyMaterial::~BrachyMaterial()
{;}

void BrachyMaterial::DefineMaterials()
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

  A=26.98*g/mole;
  G4Element* elAl = new G4Element("Aluminum","Al", Z = 13.,A);

  A = 12.011*g/mole;
  G4Element* elC = new G4Element("Carbon","C",Z = 6.,A);
 
  A = 22.99*g/mole;
  G4Element* elNa = new G4Element("Sodium","Na",Z = 11.,A);
 
  A = 24.305*g/mole;
  G4Element* elMg = new G4Element("Magnesium","Mg",Z = 12.,A);

  A = 30.974*g/mole;
  G4Element* elP = new G4Element("Phosphorus","P",Z = 15.,A);
 
  A = 32.06*g/mole;
  G4Element* elS = new G4Element("Sulfur","S",Z = 16.,A);
 
  A = 35.453*g/mole;
  G4Element* elCl = new G4Element("Chlorine","Cl",Z = 17.,A);
 
  A = 39.098*g/mole;
  G4Element* elK = new G4Element("Potassium","K",Z = 19.,A);

  A = 40.08*g/mole;
  G4Element* elCa = new G4Element("Calcium","Ca",Z = 20.,A);
  
  A = 65.38*g/mole;
  G4Element* elZn = new G4Element("Zinc","Zn",Z = 30.,A);

  A  =  54.94*g/mole;
  G4Element* elMn   =  new G4Element("Manganese","Mn",Z = 25.,A);
 
  A = 28.09*g/mole;
  G4Element* elSi  = new G4Element("Silicon","Si",Z = 14.,A);

  A = 52.00*g/mole;
  G4Element* elCr  = new G4Element("Chromium","Cr",Z = 24.,A);

  A = 58.70*g/mole;
  G4Element* elNi  = new G4Element("Nickel","Ni",Z = 28.,A);

  A = 55.85*g/mole;
  G4Element* elFe  = new G4Element("Iron","Fe",Z = 26.,A);
 
  A = 183.84* g/mole;
  d = 19.3*g/cm3;
  matW = new G4Material("Tungsten",Z = 74.,A,d);

   // Perspex, plexiglass, lucite 
  d = 1.19*g/cm3;
  matplexiglass = new G4Material("Plexiglass",d,3);
  matplexiglass->AddElement(elH,0.08);
  matplexiglass->AddElement(elC,0.60);
  matplexiglass->AddElement(elO,0.32);
 
  // Lead material
  A = 207.19*g/mole;
  Z = 82;
  d = 11.35*g/cm3;
  matPb = new G4Material("Lead",Z,A,d);

  // Iridium (Medical Physics, Vol 25, No 10, Oct 1998)
  d = 22.42*g/cm3;
  A = 191.96260*g/mole ;
  Z = 77;
  matir192 = new G4Material("Iridium",Z,A,d);

  //titanium
  A = 47.88*g/mole;
  d = 4.50*g/cm3;
  Titanium = new G4Material("titanium" ,Z = 22.,A,d);
 
  //silver
  A = 107.87*g/mole;
  d = 10.49*g/cm3;
  Z = 22.0;
  matAg = new G4Material("Silver", Z, A, d);

  // Air material
  d = 1.290*mg/cm3;
  matAir = new G4Material("Air",d,2);
  matAir->AddElement(elN,0.7);
  matAir->AddElement(elO,0.3);

  // Water
  d = 1.000*g/cm3;
  matH2O = new G4Material("Water",d,2);
  matH2O->AddElement(elH,2);
  matH2O->AddElement(elO,1);
  matH2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

  //soft tissue(http://www.nist.gov)
  d = 1.0*g/cm3;
  soft = new G4Material("tissue",d,13);
  soft->AddElement(elH,0.104472);
  soft->AddElement(elC,0.23219);
  soft->AddElement(elN,0.02488);
  soft->AddElement(elO,0.630238);
  soft->AddElement(elNa,0.00113);
  soft->AddElement(elMg,0.00013);
  soft->AddElement(elP,0.00133);
  soft->AddElement(elS,0.00199);
  soft->AddElement(elCl,0.00134);
  soft->AddElement(elK,0.00199);
  soft->AddElement(elCa,0.00023);
  soft->AddElement(elFe,0.00005);
  soft->AddElement(elZn,0.00003); 
 
  // Stainless steel (Medical Physics, Vol 25, No 10, Oct 1998)
  d = 8.02*g/cm3 ;
  matsteel = new G4Material("Stainless steel",d,5);
  matsteel->AddElement(elMn, 0.02);
  matsteel->AddElement(elSi, 0.01);
  matsteel->AddElement(elCr, 0.19);
  matsteel->AddElement(elNi, 0.10);
  matsteel->AddElement(elFe, 0.68);

 //Define Stainless-steel-304 - Flexi source
  d = 7.999*g/cm3 ;
  mat304steel = new G4Material("Stainless steel 304",d,6);
  mat304steel->AddElement(elMn, 0.02);
  mat304steel->AddElement(elSi, 0.01);
  mat304steel->AddElement(elCr, 0.19);
  mat304steel->AddElement(elNi, 0.10);
  mat304steel->AddElement(elFe, 0.6792);
  mat304steel->AddElement(elC, 0.0008);
 
  //gold 
  A = 196.97*g/mole;
  d = 19.32*g/cm3;
  gold = new G4Material("gold",Z = 79.,A,d);

  //Iodine Core
  A = 124.9*g/mole;
  d = 4.862*g/cm3;
  matI = new G4Material("Iodine",Z = 53.,A,d);

  //ceramic(Medical Physics, May 2000)
  d = 2.88*g/cm3;
  ceramic = new G4Material("allumina",d,2);
  ceramic->AddElement(elAl,2);
  ceramic->AddElement(elO,3);

  G4double density = universe_mean_density;
  G4double pressure = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  A=1.01*g/mole;
  Vacuum = new G4Material("Galactic", Z = 1., A,density,kStateGas,temperature,pressure);
  
  //compact bone (http://www.NIST.gov)
  d = 1.85*g/cm3;
  bone = new G4Material("bone",d,8);
  bone->AddElement(elH,0.063984);
  bone->AddElement(elC,0.278);
  bone->AddElement(elN,0.027);
  bone->AddElement(elO,0.410016);
  bone->AddElement(elMg,0.002);
  bone->AddElement(elP,0.07);
  bone->AddElement(elS,0.002);
  bone->AddElement(elCa,0.147);

  //muscle(http://www.NIST.gov)
  muscle = new G4Material("muscle",d,9);
  muscle->AddElement(elH,0.101997);
  muscle->AddElement(elC,0.123);
  muscle->AddElement(elN,0.035);
  muscle->AddElement(elNa,0.0008);
  muscle->AddElement(elO,0.729);
  muscle->AddElement(elMg,0.0002);
  muscle->AddElement(elP,0.002);
  muscle->AddElement(elS,0.005);
  muscle->AddElement(elK,0.003);
}

G4Material* BrachyMaterial::GetMat(G4String material)
{
  // Returns a material 
  G4Material* pttoMaterial = G4Material::GetMaterial(material); 
  return pttoMaterial; 
}
