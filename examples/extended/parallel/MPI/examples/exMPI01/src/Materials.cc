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
//
/// @file Materials.cc
/// @brief Define materials

#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Materials.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Materials::Materials()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Materials::~Materials()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Materials::Construct()
{
  G4double A, Z;

  // ------------------------------------------------------------------------
  // Elements
  // ------------------------------------------------------------------------
  G4Element* elH  = new G4Element("Hydrogen","H",  Z=1.,  A=1.00794*g/mole);
  G4Element* elC  = new G4Element("Carbon",  "C",  Z=6.,  A= 12.011 *g/mole);
  G4Element* elN  = new G4Element("Nitrogen","N",  Z=7.,  A= 14.00674*g/mole);
  G4Element* elO  = new G4Element("Oxygen",  "O",  Z=8.,  A= 15.9994*g/mole);
  G4Element* elNa = new G4Element("Sodium",  "Na", Z=11., A= 22.989768*g/mole);
  G4Element* elSi = new G4Element("Silicon", "Si", Z=14., A= 28.0855*g/mole);
  G4Element* elAr = new G4Element("Argon",   "Ar", Z=18., A= 39.948*g/mole);
  G4Element* elI  = new G4Element("Iodine",  "I",  Z=53., A= 126.90447*g/mole);
  G4Element* elCs = new G4Element("Cesium",  "Cs", Z=55., A= 132.90543*g/mole);

  // ------------------------------------------------------------------------
  // Materials
  // ------------------------------------------------------------------------
  G4double density, massfraction;
  G4int natoms, nel;

  // temperature of experimental hall is controlled at 20 degree.
  const G4double expTemp = STP_Temperature+20.*kelvin;

  // vacuum
  density = universe_mean_density;
  G4Material* Vacuum = new G4Material("Vacuum", density, nel=2);
  Vacuum-> AddElement(elN, .7);
  Vacuum-> AddElement(elO, .3);

  // air
  density = 1.2929e-03 *g/cm3;  // at 20 degree
  G4Material* Air = new G4Material("Air", density, nel=3,
                                   kStateGas, expTemp);
  G4double ttt = 75.47+23.20+1.28;
  Air-> AddElement(elN,  massfraction= 75.47/ttt);
  Air-> AddElement(elO,  massfraction= 23.20/ttt);
  Air-> AddElement(elAr, massfraction=  1.28/ttt);

  // Ar gas
  A = 39.948 *g/mole;
  const G4double denAr = 1.782e-03 *g/cm3 * STP_Temperature/expTemp;
  G4Material* Ar= new G4Material("ArgonGas", Z=18., A, denAr,
                                 kStateGas, expTemp);

  // ethane (C2H6)
  const G4double denEthane = 1.356e-3 *g/cm3 * STP_Temperature/expTemp;
  G4Material* Ethane= new G4Material("Ethane", denEthane, nel=2,
                                     kStateGas, expTemp);
  Ethane-> AddElement(elC, natoms=2);
  Ethane-> AddElement(elH, natoms=6);

  // Ar(50%) + ethane(50%) mixture
  density =  (denAr+denEthane)/2.;
  G4Material* ArEthane = new G4Material("ArEthane", density, nel=2,
                                        kStateGas, expTemp);
  ArEthane-> AddMaterial(Ar, massfraction= denAr/density/2.);
  ArEthane-> AddMaterial(Ethane, massfraction= denEthane/density/2.);

  // silicon
  A = 28.0855 *g/mole;
  density = 2.33 *g/cm3;
  new G4Material("SiliconWafer", Z=14., A, density);

  // alminium
  A = 26.98 *g/mole;
  density = 2.70 *g/cm3;
  new G4Material("Al", Z=13., A, density);

  // iron
  A = 55.847 *g/mole;
  density = 7.87 *g/cm3;
  new G4Material("Iron", Z=26., A, density);

  // lead
  A = 207.2 *g/mole;
  density = 11.35 *g/cm3;
  new G4Material("Lead", Z=82., A, density);

  // scintillator (Polystyene(C6H5CH=CH2))
  density = 1.032 *g/cm3;
  G4Material* Scinti = new G4Material("Scinti", density, nel=2);
  Scinti-> AddElement(elC, natoms=8);
  Scinti-> AddElement(elH, natoms=8);

  // quartz (SiO2, crystalline)
  density = 2.64 *g/cm3;
  G4Material* Quartz = new G4Material("Quartz", density, nel= 2);
  Quartz-> AddElement(elSi, natoms=1);
  Quartz-> AddElement(elO,  natoms=2);

  // NaI crystal
  density = 3.67 *g/cm3;
  G4Material* NaI = new G4Material("NaI", density, nel= 2);
  NaI-> AddElement(elNa, natoms=1);
  NaI-> AddElement(elI,  natoms=1);

  // CsI crystal
  density = 4.51 *g/cm3;
  G4Material* CsI = new G4Material("CsI", density, nel= 2);
  CsI-> AddElement(elCs, natoms=1);
  CsI-> AddElement(elI,  natoms=1);

}
