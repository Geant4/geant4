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
//
// $Id: nistComparison.cc,v 1.1 2005-03-01 12:14:41 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//
//
//  Test the construction of materials from the NIST data base
//

#include "G4NistManager.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"

int main() {

G4UnitDefinition::BuildUnitsTable();

// initialise NIST data base
//
G4NistManager*  nistMat = G4NistManager::Instance();
//
// define Elements
//
//G4Element* H  = new G4Element("Hydrogen" ,  "H",  1.,   1.01*g/mole);
//G4Element* C  = new G4Element("Carbon"   ,  "C",  6.,   2.01*g/mole);
//G4Element* N  = new G4Element("Nitrogen" ,  "N",  7.,  14.01*g/mole);
//G4Element* O  = new G4Element("Oxygen"   ,  "O",  8.,  16.00*g/mole);
//G4Element* Si = new G4Element("Silicon"  , "Si", 14.,  28.09*g/mole);
//G4Element* Ge = new G4Element("Germanium", "Ge", 32.,  72.59*g/mole);
//G4Element* Bi = new G4Element("Bismuth"  , "Bi", 83., 208.98*g/mole);

G4bool buildIsotopes;

G4Element* H  = nistMat->FindOrBuildElement ( "H", buildIsotopes=false);
G4Element* C  = nistMat->FindOrBuildElement ( "C", buildIsotopes=false);
G4Element* N  = nistMat->FindOrBuildElement ( "N", buildIsotopes=false);
G4Element* O  = nistMat->FindOrBuildElement ( "O", buildIsotopes=false);
G4Element* Si = nistMat->FindOrBuildElement ("Si", buildIsotopes=false);
G4Element* Ge = nistMat->FindOrBuildElement ("Ge", buildIsotopes=false);
G4Element* Bi = nistMat->FindOrBuildElement ("Bi", buildIsotopes=false);

//----------------------------------------------------------------------
G4int ncomponents, natoms;
G4double density, temperature, pressure;
 
// CH4
//
density     = 0.717*mg/cm3;
pressure    = 1.*atmosphere;
temperature = 273.15*kelvin;
G4Material* CH4 = new G4Material("Methane", density, ncomponents=2,
                                     kStateGas,temperature,pressure);
CH4->AddElement(H, natoms=4);
CH4->AddElement(C, natoms=1);

nistMat->FindOrBuildMaterial ("G4_METHANE");

// Air
//
density     = 1.205*mg/cm3;
pressure    = 1.*atmosphere;
temperature = 293.15*kelvin;
G4Material* Air = new G4Material("dry_Air", density, ncomponents=2,
                                     kStateGas,temperature,pressure);
Air->AddElement(N, 75*perCent);
Air->AddElement(O, 25*perCent);

nistMat->FindOrBuildMaterial ("G4_AIR");

// CO2
//
density     = 1.977*mg/cm3;
pressure    = 1.*atmosphere;
temperature = 273.15*kelvin;
G4Material* CO2 = new G4Material("Carbonic gas", density, ncomponents=2,
                                     kStateGas,temperature,pressure);
CO2->AddElement(C, natoms=1);
CO2->AddElement(O, natoms=2);

nistMat->FindOrBuildMaterial ("G4_CARBON_DIOXIDE");

// C8H18
//
density = 0.703*g/cm3;
G4Material* C8H18 = new G4Material("liquid Octane", density, ncomponents=2);
C8H18->AddElement(H, natoms=18);
C8H18->AddElement(C, natoms=8);

nistMat->FindOrBuildMaterial ("G4_OCTANE");

// H2O
//
density = 1.000*g/cm3;
G4Material* H2O = new G4Material("Water", density, ncomponents=2);
H2O->AddElement(H, natoms=2);
H2O->AddElement(O, natoms=1);

nistMat->FindOrBuildMaterial ("G4_WATER");

// Mylar
//
density = 1.390*g/cm3;
G4Material* mylar = new G4Material("Mylar", density, ncomponents=3);
mylar->AddElement(H, natoms=4);
mylar->AddElement(C, natoms=5);
mylar->AddElement(O, natoms=2);

nistMat->FindOrBuildMaterial ("G4_MYLAR");

// SiO2
// 
density = 2.320*g/cm3;
G4Material* SiO2 = new G4Material("Quartz", density, ncomponents=2);
SiO2->AddElement(O , natoms=2);
SiO2->AddElement(Si, natoms=1);

nistMat->FindOrBuildMaterial ("G4_SILICON_DIOXIDE");

// BGO
//
density = 7.130*g/cm3;
G4Material* BGO = new G4Material("BGO", density, ncomponents=3);
BGO->AddElement(O , natoms=12);
BGO->AddElement(Ge, natoms=3);
BGO->AddElement(Bi, natoms=4);

nistMat->FindOrBuildMaterial ("G4_BGO");

//----------------------------------------------------------------------

// Print the table of materials
//
G4cout << *(G4Material::GetMaterialTable()) << G4endl;
             
return EXIT_SUCCESS;
}
