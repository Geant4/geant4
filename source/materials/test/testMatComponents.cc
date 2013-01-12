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
// $Id: G4MaterialsTest.cc,v 1.9 2006-06-29 19:13:14 gunter Exp $
//
// 
// ------------------------------------------------------------
//
//
//  This program illustrates the different ways to define materials
//

#include "G4ios.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <map>
#include <iomanip>

int main() {

G4String name, symbol;             // a=mass of a mole;
G4double a, z, density;            // z=mean number of protons;  


G4int ncomponents, natoms;
G4double fractionmass;

G4UnitDefinition::BuildUnitsTable();

//
// define Elements
//

a = 1.01*g/mole;
G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

a = 12.01*g/mole;
G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

a = 14.01*g/mole;
G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

a = 16.00*g/mole;
G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

a = 28.09*g/mole;
G4Element* elSi = new G4Element(name="Silicon",symbol="Si" , z= 14., a);

//
// define a material from elements.   case 1: chemical molecule
//
 
density = 1.000*g/cm3;
G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
H2O->AddElement(elH, natoms=2);
H2O->AddElement(elO, natoms=1);

density = 2.200*g/cm3;
G4Material* SiO2 = new G4Material(name="quartz", density, ncomponents=2);
SiO2->AddElement(elSi, natoms=1);
SiO2->AddElement(elO , natoms=2);

//
// define a material from elements.   case 2: mixture by fractional mass
//

density = 1.290*mg/cm3;
G4Material* Air = new G4Material(name="Air  "  , density, ncomponents=2);
Air->AddElement(elN, fractionmass=0.7);
Air->AddElement(elO, fractionmass=0.3);

//
// define a material from elements and/or others materials (mixture of mixtures)
//

density = 0.200*g/cm3;
G4Material* Aerog = new G4Material(name="Aerogel", density, ncomponents=3);
Aerog->AddMaterial(SiO2, fractionmass=0.625);
Aerog->AddMaterial(H2O , fractionmass=0.374);
Aerog->AddElement (elC , fractionmass=0.1*perCent);



//
// Print the table of materials
//
G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//
// mass of molecule
//
G4Material* material = H2O;
G4String matName = material->GetName();
G4double mass = material->GetMassOfMolecule();

 G4cout << "\n---> " << matName << ": mass of molecule "
  << mass << " = " << mass/g << " g = " << mass/CLHEP::amu << " amu"
        << G4endl;

//
// print MatComponents
//
std::map<G4Material*,G4double> matComponents;
std::map<G4Material*,G4double>::iterator it;

material = Aerog;
matName = material->GetName();

matComponents = material->GetMatComponents();

 G4cout << "\n---> " << matName << ": material Components" << G4endl;

 for (it = matComponents.begin(); it != matComponents.end(); it++) {
    G4Material* mat  = it->first;
    G4String name = mat->GetName();
    G4double fraction   = it->second;         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << fraction
           << G4endl;	   
 }
             
return EXIT_SUCCESS;
}
