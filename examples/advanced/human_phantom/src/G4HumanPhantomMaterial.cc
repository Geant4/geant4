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
//  S.Guatelli
//
//    *******************************
//    *                             *
//    *    G4HumanPhantomMaterial.cc        *
//    *                             *
//    *******************************
//
// $Id: G4HumanPhantomMaterial.cc,v 1.1 2006-09-26 17:29:00 guatelli Exp $
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
#include "G4HumanPhantomMaterial.hh"

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

   A = 18.998*g/mole;
  G4Element* elF = new G4Element("Fluorine","F",Z = 9.,A);

  A = 22.99*g/mole;
  G4Element* elNa = new G4Element("Sodium","Na",Z = 11.,A);

  A = 24.305*g/mole;
  G4Element* elMg = new G4Element("Magnesium","Mg",Z = 12.,A);

  //  A=26.98*g/mole;
  //G4Element* elAl = new G4Element("Aluminum","Al", Z = 13.,A);

  A = 28.09*g/mole;
  G4Element* elSi  = new G4Element("Silicon","Si",Z = 14.,A);

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


  //  A  =  54.94*g/mole;
  //G4Element* elMn   =  new G4Element("Manganese","Mn",Z = 25.,A);
 

  //A = 52.00*g/mole;
  // G4Element* elCr  = new G4Element("Chromium","Cr",Z = 24.,A);

  //A = 58.70*g/mole;
  // G4Element* elNi  = new G4Element("Nickel","Ni",Z = 28.,A);

 
  //  A = 183.84* g/mole;
  //d = 19.3*g/cm3;
  //matW = new G4Material("Tungsten",Z = 74.,A,d);

   // Perspex, plexiglass, lucite 
  //d = 1.19*g/cm3;
  //matplexiglass = new G4Material("Plexiglass",d,3);
  //matplexiglass->AddElement(elH,0.08);
  //matplexiglass->AddElement(elC,0.60);
  //matplexiglass->AddElement(elO,0.32);
 
  // Lead material
  //A = 207.19*g/mole;
  // Z = 82;
  //d = 11.35*g/cm3;
  //matPb = new G4Material("Lead",Z,A,d);

  // Air material
  // d = 1.290*mg/cm3;
  // G4Material* matAir = new G4Material("Air",d,2);
  //matAir->AddElement(elN,0.7);
  //matAir->AddElement(elO,0.3);

  // Water
  d = 1.000*g/cm3;
  matH2O = new G4Material("Water",d,2);
  matH2O->AddElement(elH,2);
  matH2O->AddElement(elO,1);
  matH2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);


  //soft tissue(http://www.nist.gov)
  d = 1.04*g/cm3;
  soft = new G4Material("soft_tissue",d,16);
  soft->AddElement(elH,0.10454);
  soft->AddElement(elC,0.22663);
  soft->AddElement(elN,0.02490);
  soft->AddElement(elO,0.63525);
  soft->AddElement(elNa,0.00112);
  soft->AddElement(elMg,0.00013);
  soft->AddElement(elSi,0.00030);
  soft->AddElement(elP,0.00134);
  soft->AddElement(elS,0.00204);
  soft->AddElement(elCl,0.00133);
  soft->AddElement(elK,0.00208);
  soft->AddElement(elCa,0.00024);
  soft->AddElement(elFe,0.00005);
  soft->AddElement(elZn,0.00003);
  soft->AddElement(elRb,0.00001);
  soft->AddElement(elZr,0.00001);

  /*
  // Stainless steel (Medical Physics, Vol 25, No 10, Oct 1998)
  d = 8.02*g/cm3 ;
  matsteel = new G4Material("Stainless steel",d,5);
  matsteel->AddElement(elMn, 0.02);
  matsteel->AddElement(elSi, 0.01);
  matsteel->AddElement(elCr, 0.19);
  matsteel->AddElement(elNi, 0.10);
  matsteel->AddElement(elFe, 0.68);
 
  //gold(chimica degli elementi N.N Greenwood,A.Earnshaw)
  A = 196.97*g/mole;
  d = 19.32*g/cm3;
  gold = new G4Material("gold",Z = 79.,A,d);

  //IodiumCore(chimica degli elementi N.N Greenwood,A.Earnshaw)
  A = 124.9*g/mole;
  d = 4.862*g/cm3;
  matI = new G4Material("Iodium",Z = 53.,A,d);

  //ceramic(Medical Physics, May 2000)
  d = 2.88*g/cm3;
  ceramic = new G4Material("allumina",d,2);
  ceramic->AddElement(elAl,2);
  ceramic->AddElement(elO,3);

  G4double density = universe_mean_density;
  G4double pressure = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  A=1.01*g/mole;
  Vacuum = new G4Material("Galactic", Z = 1., A,
				     density,kStateGas,temperature,pressure);

  */ 
 
  // Skeleton

  d = 1.4*g/cm3;
  skeleton = new G4Material("skeleton",d,18);
  skeleton -> AddElement(elH,0.07337);
  skeleton -> AddElement(elC,0.25475);
  skeleton -> AddElement(elN,0.03057);
  skeleton -> AddElement(elO,0.47893);
  skeleton -> AddElement(elF,0.00025);
  skeleton -> AddElement(elNa,0.00326); 
  skeleton -> AddElement(elMg,0.00112);
  skeleton -> AddElement(elSi,0.00002); 
  skeleton -> AddElement(elP,0.05095);
  skeleton -> AddElement(elS,0.00173);
  skeleton -> AddElement(elCl,0.00143);
  skeleton -> AddElement(elK,0.00153);
  skeleton -> AddElement(elCa,0.10190);
  skeleton -> AddElement(elFe,0.00008);
  skeleton -> AddElement(elZn,0.00005);
  skeleton -> AddElement(elRb,0.00002);
  skeleton -> AddElement(elSr,0.00003);
  skeleton -> AddElement(elPb,0.00001);
 
  d = 0.296 *g/cm3;
  lung = new G4Material("lung_material", d,15);
  lung -> AddElement(elH, 0.10134);
  lung -> AddElement(elC, 0.10238);
  lung -> AddElement(elN,0.02866);
  lung -> AddElement(elO,0.75752);
  lung -> AddElement(elNa,0.00184);
  lung -> AddElement(elMg,0.00007);
  lung -> AddElement(elSi,0.00006);
  lung -> AddElement(elP,0.00080);
  lung -> AddElement(elS,0.00225);
  lung -> AddElement(elCl,0.00266);
  lung -> AddElement(elK,0.00194);
  lung -> AddElement(elCa,0.00009);
  lung -> AddElement(elFe,0.00037);
  lung -> AddElement(elZn,0.00001);
  lung -> AddElement(elRb,0.00001);

  
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
  
 /*
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
  */
}

G4Material* G4HumanPhantomMaterial::GetMaterial(G4String material)
{
  // Returns a material 
  G4Material* pttoMaterial = G4Material::GetMaterial(material); 
  return pttoMaterial; 
}
