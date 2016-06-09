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
//    *    RemSimMaterial.cc        *
//    *                             *
//    *******************************
//
// $Id: RemSimMaterial.cc,v 1.7 2006-06-29 16:23:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "globals.hh"
#include "G4MaterialTable.hh"
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "RemSimMaterial.hh"

RemSimMaterial::RemSimMaterial():
  matPb(0), matAir(0), matH2O(0), Al(0), nylon(0), mylar(0), 
  beta(0), nextel(0), kevlar(0),
  vacuum(0), betaCloth(0), eterogeneousNextel(0), kevlarVacuum(0),
  polyethylene(0), polyacrylate(0), evoh(0), nomex(0), nomexAir(0), 
  kevlarAir(0), moon(0)
{;}

RemSimMaterial::~RemSimMaterial()
{
  delete moon;
  delete kevlarAir;
  delete nomexAir;
  delete nomex;
  delete evoh;
  delete polyacrylate;
  delete polyethylene;
  delete kevlarVacuum;
  delete eterogeneousNextel;
  delete betaCloth;
  delete vacuum;
  delete kevlar;
  delete nextel;
  delete beta;
  delete mylar;
  delete nylon;
  delete Al;
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
  
  A = 10.811*g/mole;
  G4Element* elB = new G4Element ("Boro","B",Z = 5.,A);
 
  A = 14.01*g/mole;
  G4Element* elN = new G4Element("Nitrogen","N",Z = 7.,A);

  A = 16.00*g/mole;

  G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A);

  A=26.98*g/mole;
  G4Element* elAl = new G4Element("Aluminum","Al", Z = 13.,A);

  A = 12.011*g/mole;
  G4Element* elC = new G4Element("Carbon","C",Z = 6.,A);
  
  A = 24.305*g/mole;
  G4Element* elMg = new G4Element("Magnesium","Mg",Z = 12.,A);
 
  A = 35.453*g/mole;
  G4Element* elCl = new G4Element("Chlorine","Cl",Z = 17.,A);

  A = 40.08*g/mole;
  G4Element* elCa = new G4Element("Calcium","Ca",Z = 20.,A);
 
  A = 28.09*g/mole;
  G4Element* elSi  = new G4Element("Silicon","Si",Z = 14.,A);

  A = 55.85*g/mole;
  G4Element* elFe  = new G4Element("Iron","Fe",Z = 26.,A);

  //Lead material
  A = 207.19*g/mole;
  Z = 82;
  d = 11.35*g/cm3;
  matPb = new G4Material("Lead",Z,A,d);
 
  // Air material
  G4double airDensity = 1.290*mg/cm3;
  G4Material* matAir = new G4Material("Air",airDensity,2);
  matAir -> AddElement(elN,0.7);
  matAir -> AddElement(elO,0.3);

  // Water
  d = 1.000*g/cm3;
  matH2O = new G4Material("Water",d,2);
  matH2O -> AddElement(elH,2);
  matH2O -> AddElement(elO,1);
  matH2O -> GetIonisation() -> SetMeanExcitationEnergy(75.0*eV);
  
  Al = new G4Material("Aluminium", Z= 13., A= 26.98*g/mole, d = 2.7*g/cm3);

  //nylon (alenia spazio)
  d = 1.14 *g/cm3;
  nylon = new G4Material("nylon",d,4);
  nylon -> AddElement(elH,0.108);
  nylon -> AddElement(elC,0.68);
  nylon -> AddElement(elN,0.099);
  nylon -> AddElement(elO,0.113);

  //mylar (alenia spazio)
  d= 1.4 *g/cm3;
  mylar = new G4Material("mylar",d,3);
  mylar -> AddElement(elH,0.042);
  mylar -> AddElement(elC,0.625);
  mylar -> AddElement(elO,0.333);

  //beta cloth
  G4double betaDensity = 2.3 *g/cm3;
  beta = new G4Material("beta",betaDensity,2);
  beta -> AddElement(elO,0.53);
  beta -> AddElement(elSi,0.47);
 
  G4double nextelDensity = 2.7 * g/cm3;
  nextel = new G4Material("nextel",nextelDensity,4);
  nextel -> AddElement(elB,0.04);
  nextel -> AddElement(elO,0.52);
  nextel -> AddElement(elAl,0.33);
  nextel -> AddElement(elSi,0.11);

  //kevlar
  G4double kevlarDensity = 1.44 *g/cm3;
  kevlar = new G4Material("kevlar",d,4);
  kevlar -> AddElement(elH,0.04);
  kevlar -> AddElement(elC,0.71);
  kevlar -> AddElement(elO,0.12);
  kevlar -> AddElement(elN,0.13);
  
  //vacuum
  G4double vacuumDensity = 1.e-25 *g/cm3;
  G4double pressure = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  vacuum = new G4Material("Galactic", Z=1., A=1.01*g/mole,
			   vacuumDensity,kStateGas,temperature,pressure);

  d = (vacuumDensity*0.44)+ (betaDensity*0.56);
  betaCloth = new G4Material("betacloth",d,2);
  betaCloth -> AddMaterial (beta, 0.56);
  betaCloth -> AddMaterial (vacuum,0.44);

  d = (vacuumDensity*0.73)+ (nextelDensity*0.27); 
  eterogeneousNextel = new G4Material("Nextel312AF62",d,2);
  eterogeneousNextel -> AddMaterial (vacuum, 0.73);
  eterogeneousNextel -> AddMaterial (nextel, 0.27);

  d = (vacuumDensity*0.44)+ (kevlarDensity*0.56);
  kevlarVacuum = new G4Material("kevlarVacuum",d,2);
  kevlarVacuum -> AddMaterial (vacuum, 0.44);
  kevlarVacuum -> AddMaterial (kevlar, 0.56);

  d = 0.94 * g/cm3;
  polyethylene = new G4Material("polyethylene",d,2);
  polyethylene -> AddElement(elH,0.14);
  polyethylene -> AddElement(elC,0.86);
  
  d = 1.19 * g/cm3;
  polyacrylate = new G4Material("polyacrylate",d,3);
  polyacrylate -> AddElement(elH,0.08);
  polyacrylate -> AddElement(elC,0.60);
  polyacrylate -> AddElement(elO,0.32);
 
  d = 1.17 * g/cm3;
  evoh = new G4Material("evoh",d,3);
  evoh -> AddElement(elH,0.11);
  evoh -> AddElement(elC,0.67);
  evoh -> AddElement(elO,0.22);
  
  G4double nomexDensity = 0.98 * g/cm3;
  nomex = new G4Material("nomex",nomexDensity,5);
  nomex -> AddElement(elH,0.04);
  nomex -> AddElement(elC,0.54);
  nomex -> AddElement(elN,0.09);
  nomex -> AddElement(elO,0.10);
  nomex -> AddElement(elCl,0.23);

  d = 0.45*(nomexDensity)+ 0.55*(airDensity);
  nomexAir = new G4Material("nomexAir",d,2);
  nomexAir -> AddMaterial(nomex,0.45);
  nomexAir -> AddMaterial(matAir,0.55);

  d =0.56*(kevlarDensity)+ 0.44*(airDensity);
  kevlarAir = new G4Material("kevlarAir",d,2);
  kevlarAir -> AddMaterial(kevlar,0.56);
  kevlarAir -> AddMaterial(matAir,0.44);

  d = 1.5*g/cm3;
  moon =  new G4Material("moon",d,6);
  moon -> AddElement(elO,0.45);
  moon -> AddElement(elMg,0.05);
  moon -> AddElement(elAl,0.13);
  moon -> AddElement(elSi,0.21);
  moon -> AddElement(elCa,0.10);
  moon -> AddElement(elFe,0.06); 
}

G4Material* RemSimMaterial::GetMaterial(G4String material)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(material); 
  return pttoMaterial; 
}
