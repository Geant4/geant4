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
#include "MaterialTest47.hh"

#include "G4UnitsTable.hh"
#include "G4Material.hh"

MaterialTest47::MaterialTest47() {}

MaterialTest47::~MaterialTest47() {}

G4Material* MaterialTest47::GetMaterial(const G4String& name) { 

  G4Material* ma = 0;
  if (name == "Al") {
    ma = new G4Material("Aluminium", 13, 26.98*g/mole, 2.7*g/cm3);
  } else if (name == "Au") {
    ma = new G4Material("Gold", 79, 196.97*g/mole, 18.85*g/cm3);
  } else if (name == "Be") {
    ma = new G4Material("Beryllium", 4, 9.0122*g/mole, 1.848*g/cm3);
  } else if (name == "Brass") {
    G4Element *elCu = new G4Element("Copper", "Cu", 29, 63.546*g/mole);
    G4Element *elZn = new G4Element("Zinc",   "Zn", 30, 65.39*g/mole);
    ma = new G4Material("Brass", 8.53*g/cm3, 2);
    ma->AddElement(elCu, 0.70);
    ma->AddElement(elZn, 0.30);
  } else if (name == "Cd") {
    ma = new G4Material("Cadmium", 48, 112.41*g/mole, 8.63*g/cm3);
  } else if (name == "C") {
    ma = new G4Material("Carbon", 6, 12.011*g/mole, 2.265*g/cm3);
  } else if (name == "Cu") {
    ma = new G4Material("Copper", 29, 63.546*g/mole, 8.96*g/cm3);
  } else if (name == "D") {
    ma = new G4Material("Deuterium", 1, 2.01*g/mole, 162*mg/cm3);
  } else if (name == "H") {
    ma = new G4Material("Hydrogen", 1, 1.00794*g/mole, 70.8*mg/cm3);
  } else if (name == "Fe") {
    ma = new G4Material("Iron", 26, 55.85*g/mole, 7.87*g/cm3);
  } else if (name == "Nb") {
    ma = new G4Material("Niobium", 41, 92.906*g/mole, 8.55*g/cm3);
  } else if (name == "Pb") {
    ma = new G4Material("Lead", 82, 207.19*g/mole, 11.35*g/cm3);
  } else if (name == "PbWO4") {
    G4Element *elPb = new G4Element("Lead", "Pb", 82, 207.19*g/mole);
    G4Element *elW  = new G4Element("Tungsten", "W", 74, 183.85*g/mole);
    G4Element *elO  = new G4Element("Oxygen", "O", 8, 15.999*g/mole);
    ma = new G4Material("PbWO4", 8.28*g/cm3, 3);
    ma->AddElement(elPb, 0.45532661);
    ma->AddElement(elW,  0.40403397);
    ma->AddElement(elO,  0.14063942);
  } else if (name == "Sn") {
    ma = new G4Material("Tin", 50, 118.69*g/mole, 7.31*g/cm3);
  } else if (name == "Ta") {
    ma = new G4Material("Tantalum", 73, 180.9479*g/mole, 16.65*g/cm3);
  } else if (name == "Ti") {
    ma = new G4Material("Titanium", 22, 47.88*g/mole, 4.53*g/cm3);
  } else if (name == "U") {
    ma = new G4Material("Uranium", 92, 238.03*g/mole, 18.95*g/cm3);
  }

  if (ma != 0) G4cout << "Material is selected: " << ma->GetName() << G4endl;
  else         G4cout << "Cannot initialize material " << name << G4endl;

  return ma;
}	

  






