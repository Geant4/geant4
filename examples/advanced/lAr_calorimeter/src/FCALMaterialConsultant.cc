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
//   Author:           Mathieu Fontaine, Rachid Mazini
//                     fontaine@lps.umontreal.ca  Rachid.Mazini@cern.ch
//
//   Language:         C++
//   Tested on:        g++
//   Prerequisites:    None
//   Purpose:          This is the place, where all the materials get defined.
//                     Instead of coding those materials locally, where they
//                     are needed, it is much easier to maintain, if we keep
//                     all materials for a detector component in one place.
//                     Everybody who needs some of these parameters, can
//                     query the FCALMaterialConsultant.
//                 --> This class is made a singleton by making the
//                     constructor private and hiding it behind the
//                     construct() method, which creates a first instance
//                     if it does not exist. This is to prevent multiple
//                     copies of this consultant with potentially different
//                     contents (once the data is loaded from files and/or
//                     can be changed by user interaction).
//                 --> The method Material is provided to access to the data
//                     stored, a routine ShowMeAllYouKnow can be queried to
//                     dump the entire knowledge of this consultant.
//
//                   * Ideas on how the theFCALMaterialConsultant pointer
//                     is made static are borrowed from G4VisManager.
//
//----------------------------------------------------------------------------------

#include "FCALMaterialConsultant.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

FCALMaterialConsultant *FCALMaterialConsultant::theFCALMaterialConsultant = NULL;

FCALMaterialConsultant::FCALMaterialConsultant()
{
  G4double a,z,density,fractionmass;
  G4String name,symbol;
  G4int nel,natoms;

  //------------
  // elements
  //------------

  a=1.01*g/mole;
  elH=new G4Element(name="Hydrogen",symbol="H2",z=1.,a);

  a=2.01*g/mole;
  elD=new G4Element(name="Deuterium",symbol="D",z=1.,a);

  a=4.*g/mole;
  elHe=new G4Element(name="Helium",symbol="He",z=2.,a);

  a=6.94*g/mole;
  elLi=new G4Element(name="Lithium",symbol="Li",z=3.,a);

  a=9.01*g/mole;
  elBe=new G4Element(name="Berillium",symbol="Be",z=4.,a);

  a=12.01*g/mole;
  elC=new G4Element(name="Carbon",symbol="C",z=6.,a);

  a=14.01*g/mole;
  elN=new G4Element(name="Nitrogen",symbol="N2",z=7.,a);

  a=16.*g/mole;
  elO=new G4Element(name="Oxygen",symbol="O2",z=8.,a);

  a=20.18*g/mole;
  elNe=new G4Element(name="Neon",symbol="Ne",z=10.,a);

  a=22.99*g/mole;
  elNa=new G4Element(name="Sodium",symbol="Na",z=11.,a);

  a=26.98*g/mole;
  elAl=new G4Element(name="Aluminium",symbol="Al",z=13.,a);

  a=28.085*g/mole;
  elSi=new G4Element(name="Silicon",symbol="Si",z=14.,a);

  a=40.08*g/mole;
  elCa=new G4Element(name="Calcium",symbol="Ca",z=20.,a);

  a=55.850*g/mole;
  elFe=new G4Element(name="Iron",symbol="Fe",z=26.,a);

  a=63.54*g/mole;
  elCu=new G4Element(name="Copper",symbol="Cu",z=29.,a);

  a=183.85*g/mole;
  elW=new G4Element(name="Tungstenm",symbol="W",z=74.,a);

  a=207.19*g/mole;
  elPb=new G4Element(name="Lead",symbol="Pb",z=82.,a);

  a=238.03*g/mole;
  elU=new G4Element(name="Uranium",symbol="U",z=92.,a);


  //-------------------
  // simple materials
  //-------------------

  density = 2.7*g/cm3;
  a = 26.98*g/mole;
  Aluminium = new G4Material(name="Aluminium",z=13.,a,density);
  
  density = 7.87*g/cm3;
  a = 55.85*g/mole;
  Iron = new G4Material(name="Iron",z=26.,a,density);

  density = 8.96*g/cm3;
  a = 63.54*g/mole;
  Copper = new G4Material(name="Copper",z=29.,a,density);

  density = 19.3*g/cm3;
  a = 183.85*g/mole;
  Tungsten = new G4Material(name="Tungsten",z=74.,a,density);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  Lead = new G4Material(name="Lead",z=82.,a,density);

  density = 1.4*g/cm3;
  a = 39.95*g/mole;
  LiquidArgon = new G4Material(name="LiquidArgon",z=18.,a,density);

  density = 0.002*g/cm3;
  a = 39.95*g/mole;
  ArgonGas = new G4Material(name="ArgonGas",z=18.,a,density);

  density = 8.96*g/cm3;
  a = 58.69*g/mole;
  Nickel = new G4Material(name="Nickel",z=28.,a,density);

  
  //------------------
  // mixtures
  //------------------

  density = 1.290*mg/cm3;
  Air = new G4Material(name="Air",density, nel=2);
  Air->AddElement(elN, 0.7);
  Air->AddElement(elO, 0.3);

  RhoaCell = Air;


  density              = 1.e-5*g/cm3;
  G4double pressure    = 2.e-2*bar;
  G4double temperature = STP_Temperature;         //from PhysicalConstants.h
  Vacuum = new G4Material(name="Vacuum", density, nel=1,
			  kStateGas,temperature,pressure);
  Vacuum->AddMaterial(Air, fractionmass=1.);


  density = 0.002*g/cm3;
  CO2 = new G4Material(name="CO2",density,nel=2);
  CO2->AddElement(elC, natoms=1);
  CO2->AddElement(elO, natoms=2);

  density = 1.42*g/cm3;
  Kapton = new G4Material(name="Kapton",density, nel=4);
  Kapton->AddElement(elH, fractionmass = 0.0273);
  Kapton->AddElement(elC, fractionmass = 0.7213);
  Kapton->AddElement(elN, fractionmass = 0.0765);
  Kapton->AddElement(elO, fractionmass = 0.1749);

  density = 1.032*g/cm3;
  Polystyrene = new G4Material(name="Polystyrene",density,nel=2);
  Polystyrene->AddElement(elC, natoms=8);
  Polystyrene->AddElement(elH, natoms=8);

  density = 5.185*g/cm3;
  FCAL1CuArKap = new G4Material(name="FCAL1CuArKap",density,nel=3);
  FCAL1CuArKap->AddMaterial(Copper, fractionmass = 0.864);
  FCAL1CuArKap->AddMaterial(Kapton, fractionmass = 0.068);
  FCAL1CuArKap->AddMaterial(LiquidArgon, fractionmass = 0.068);

  density = 8.701*g/cm3;
  FCAL1CuAr = new G4Material(name="FCAL1CuAr",density,nel=2);
  FCAL1CuAr->AddMaterial(Copper, fractionmass = 0.994); 
  FCAL1CuAr->AddMaterial(LiquidArgon, fractionmass = 0.006);

  density = 5.185*g/cm3;
  FCAL2CuArKap = new G4Material(name="FCAL2CuArKap",density,nel=3);
  FCAL2CuArKap->AddMaterial(Copper, fractionmass = 0.864);
  FCAL2CuArKap->AddMaterial(Kapton, fractionmass = 0.068);
  FCAL2CuArKap->AddMaterial(LiquidArgon, fractionmass = 0.068);

  density = 18.6*g/cm3;
  FCAL2WFeNi = new G4Material(name="FCAL2WFeNi",density,nel=3);
  FCAL2WFeNi->AddMaterial(Tungsten, fractionmass = 0.97);
  FCAL2WFeNi->AddMaterial(Iron, fractionmass = 0.01);
  FCAL2WFeNi->AddMaterial(Nickel, fractionmass = 0.02);
  
  density = 15.366*g/cm3;
  FCAL2WFeNiCuAr = new G4Material(name="FCAL2WFeNiCuAr",density,nel=3);
  FCAL2WFeNiCuAr->AddMaterial(FCAL2WFeNi, fractionmass = 0.913);
  FCAL2WFeNiCuAr->AddMaterial(Copper, fractionmass = 0.077);
  FCAL2WFeNiCuAr->AddMaterial(LiquidArgon, fractionmass = 0.01);

  density = 0.002*g/cm3;
  MWPCArCO2 = new G4Material(name="MWPCArCO2",density,nel=2);
  MWPCArCO2->AddMaterial(CO2, fractionmass = 0.2);
  MWPCArCO2->AddMaterial(ArgonGas, fractionmass = 0.8);


  // must  check recipe for concrete

  density = 2.5*g/cm3;
  ShieldingConcrete = new G4Material(name="ShieldingConcrete",density,nel=6);
  ShieldingConcrete->AddElement(elO, fractionmass = 0.52);
  ShieldingConcrete->AddElement(elSi, fractionmass = 0.325);
  ShieldingConcrete->AddElement(elCa, fractionmass = 0.06);
  ShieldingConcrete->AddElement(elNa, fractionmass = 0.015);
  ShieldingConcrete->AddElement(elFe, fractionmass = 0.04);
  ShieldingConcrete->AddElement(elAl, fractionmass = 0.04);

  // must have the right composition for stainless steel

  density = 8.96*g/cm3;
  StainlessSteel = new G4Material(name="StainlessSteel",density,nel=1);
  StainlessSteel->AddElement(elO, fractionmass = 1.);

}

FCALMaterialConsultant * FCALMaterialConsultant::GetInstance()
{
  if (theFCALMaterialConsultant == NULL) {
    theFCALMaterialConsultant = new FCALMaterialConsultant();
  }
  return theFCALMaterialConsultant;
}

G4Material * FCALMaterialConsultant::Material(G4String what)
{
  G4Material* material = 0;
  if(what == "Air")               material = Air;
  if(what == "Vacuum")            material = Vacuum;
  if(what == "LiquidArgon")       material = LiquidArgon;
  if(what == "Aluminium")         material = Aluminium;
  if(what == "Iron")              material = Iron;
  if(what == "Copper")            material = Copper;
  if(what == "Tungsten")          material = Tungsten;
  if(what == "Lead")              material = Lead;
  if(what == "CO2")               material = CO2;
  if(what == "ArgonGas")          material = ArgonGas;
  if(what == "ShieldingConcrete") material = ShieldingConcrete;
  if(what == "Polystyrene")       material = Polystyrene;
  if(what == "StainlessSteel")    material = StainlessSteel;
  if(what == "Nickel")            material = Nickel;
  if(what == "FCAL1CuArKap")      material = FCAL1CuArKap;
  if(what == "FCAL1CuAr")         material = FCAL1CuAr;
  if(what == "FCAL2CuArKap")      material = FCAL2CuArKap;
  if(what == "FCAL2WFeNi")        material = FCAL2WFeNi;
  if(what == "FCAL2WFeNiCuAr")    material = FCAL2WFeNiCuAr;
  if(what == "MWPCArCO2")         material = MWPCArCO2;
  if(what == "RhoaCell")          material = RhoaCell;

  return material;
}
				  
