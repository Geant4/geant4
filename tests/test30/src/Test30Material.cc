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
// -------------------------------------------------------------
//      GEANT 4 class 
//
//      History: based on object model of
//      ---------- Test30Material -------
//                by Vladimir Ivanchenko, 12 March 2002 
// 
//    Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test30Material.hh"

#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test30Material::Test30Material()
{
	Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test30Material::~Test30Material()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test30Material::Initialise()
{ 

  G4std::vector<G4Material*> list;
	G4Material* ma;

  ma  = new G4Material("Be",    4.,  9.01*g/mole, 1.848*g/cm3);
	list.push_back(ma);
  ma  = new G4Material("C",     6.,  12.00*g/mole, 2.0*g/cm3);
	list.push_back(ma);
	ma  = new G4Material("Graphite",6., 12.00*g/mole, 2.265*g/cm3 );
  ma->SetChemicalFormula("Graphite");
	list.push_back(ma);
  ma  = new G4Material("Al",    13.,  26.98*g/mole,  2.7 *g/cm3);
	list.push_back(ma);
	ma  = new G4Material("Si",    14.,  28.055*g/mole, 2.33*g/cm3);  
	list.push_back(ma);
	ma  = new G4Material("LAr",   18.,  39.95*g/mole,  1.393*g/cm3);
	list.push_back(ma);
	ma  = new G4Material("LXe",   54., 131.29*g/mole,  3.02*g/cm3);
	list.push_back(ma);
	ma  = new G4Material("Fe",    26.,  55.85*g/mole,  7.87*g/cm3);
	list.push_back(ma);
	ma  = new G4Material("Cu",    29.,  63.55*g/mole,  8.96*g/cm3);
	list.push_back(ma);  
	ma  = new G4Material("Au",    79., 196.97*g/mole, 19.32*g/cm3);
  list.push_back(ma);
	ma  = new G4Material("W",     74., 183.85*g/mole, 19.30*g/cm3);
  list.push_back(ma);
	ma  = new G4Material("Pb",    82., 207.19*g/mole, 11.35*g/cm3);
	list.push_back(ma);  
	ma  = new G4Material("U",     92., 238.03*g/mole, 18.95*g/cm3);

  G4Element*   H  = new G4Element ("Hydrogen", "H",   1. ,  1.01*g/mole);
  G4Element*   N  = new G4Element ("Nitrigen", "N",   7. , 14.00*g/mole);
	G4Element*   O  = new G4Element ("Oxygen"  , "O",   8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C",   6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I  = new G4Element ("Iodide"  , "I",  53. , 126.9044*g/mole);

  ma = new G4Material("O2", 8., 16.00*g/mole, 1.1*g/cm3);
  ma->SetChemicalFormula("O_2");
	list.push_back(ma);	
  ma = new G4Material ("Water" , 1.*g/cm3, 2);
  ma->AddElement(H,2);
  ma->AddElement(O,1);
  ma->SetChemicalFormula("H_2O");
	list.push_back(ma);
  ma = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  ma->AddElement(H,6);
  ma->AddElement(C,2);
  ma->SetChemicalFormula("C_2H_6");
	list.push_back(ma);  
  ma = new G4Material ("CsI" , 4.53*g/cm3, 2);
  ma->AddElement(Cs,1);
  ma->AddElement(I,1);
  ma->SetChemicalFormula("CsI");
	list.push_back(ma);
  ma = new G4Material("Air"  , 1.290*mg/cm3, 2);
	// use fraction in mass
  ma->AddElement(N, 0.7);
  ma->AddElement(O, 0.3);
	list.push_back(ma);
	
//  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Material* Test30Material::GetMaterial(const G4String& name)
{ 

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

	G4Material* ma = G4Material::GetMaterial(name);
	
  G4cout << "Material is selected: " << ma->GetName() << G4endl;
  return ma;
}	

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  






