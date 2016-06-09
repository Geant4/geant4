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
// $Id: TiaraMaterials.cc,v 1.4 2003/06/25 09:13:07 gunter Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
//

#include "TiaraMaterials.hh"
#include "G4Element.hh"
#include "G4Material.hh"


TiaraMaterials::TiaraMaterials()
{

  G4String name, symbol;

  FillElementMap("Hydrogen", "H", 1, 1.01*g/mole);
  FillElementMap("Carbon", "C", 6, 12.01*g/mole);
  FillElementMap("Oxygen", "O", 8,  16.00*g/mole);
  FillElementMap("Natrium", "Na", 11, 22.99*g/mole);
  FillElementMap("Magnesium", "Mg", 12, 24.305*g/mole);  
  //  FillElementMap("Hg", "Hg", 80, 200.59*g/mole);
  FillElementMap("Aluminium", "Al", 14, 26.98*g/mole);
  FillElementMap("Silicon", "Si", 14, 28.09*g/mole);
  FillElementMap("K", "K", 19, 39.1*g/mole);
  FillElementMap("Calzium", "Ca", 31, 69.72*g/mole);
  FillElementMap("Iron", "Fe", 26, 55.85*g/mole);
  FillElementMap("Nitro", "N", 7 , 14.006*g/mole);

  
  G4Material *mat = 0;
  mat = CreateConcrete();
  fMapNameMaterial[mat->GetName()] = mat;
  mat = CreateIron();
  fMapNameMaterial[mat->GetName()] = mat;
  mat = CreateAir();
  fMapNameMaterial[mat->GetName()] = mat;
  mat = CreateVakuum();
  fMapNameMaterial[mat->GetName()] = mat;
}

TiaraMaterials::~TiaraMaterials(){
}


G4Material *TiaraMaterials::GetMaterial(const G4String &matName) const {
  G4Material *mat =0;

  TiaraMapNameMaterial::const_iterator matElm = 
    fMapNameMaterial.find(matName);

  if (matElm == fMapNameMaterial.end()) {
    G4cout << "TiaraMaterials::GetMaterial: material named \"" <<
      matName << "\", not known!" << G4endl;
  }
  else {
    mat = matElm->second;
  }
  return mat;
}


void TiaraMaterials::FillElementMap(const G4String &name, 
				   const G4String &symbol,
				   G4int Z,
				   G4double A) {
  TiaraMapSymbolElement::iterator it = fMapSymbolElement.find(symbol);
  if (it!=fMapSymbolElement.end()) {
    G4cout << "TiaraMaterials::FillElementMap: symbol: " 
	   << symbol << ", already defined" << G4endl;
  }
  else {
    fMapSymbolElement[symbol] = new G4Element(name, symbol, Z, A);
    if (!fMapSymbolElement[symbol]) {
      G4Exception("TiaraMaterials::FillElementMap: new failed to create G4Element!");
    }
  }
  return;
}

G4Material *TiaraMaterials::CreateMCNPConcrete(){
  typedef std::map< G4Element* , G4double > TiaraMapElementFraction;

  TiaraMapElementFraction concreteFractions; 
  concreteFractions[fMapSymbolElement["H"]] = 0.01;
  concreteFractions[fMapSymbolElement["O"]] = 0.529;
  concreteFractions[fMapSymbolElement["Na"]] = 0.016;
  //  concreteFractions[fMapSymbolElement["Hg"]] = 0.002 ;
  concreteFractions[fMapSymbolElement["Al"]] = 0.034;
  concreteFractions[fMapSymbolElement["Si"]] = 0.337 + 0.002;
  concreteFractions[fMapSymbolElement["K"]] = 0.013;
  concreteFractions[fMapSymbolElement["Ca"]] = 0.044;
  concreteFractions[fMapSymbolElement["Fe"]] = 0.014;
  concreteFractions[fMapSymbolElement["C"]] = 0.001;

  G4Material *mat =0;
  G4double density = 2.03*g/cm3;
  mat = new G4Material("concrete", density, 
		       concreteFractions.size());
  for (TiaraMapElementFraction::iterator it = concreteFractions.begin();
       it != concreteFractions.end(); ++it) {
    mat->AddElement(it->first , it->second);
  }
  
  return mat;

}

G4Material *TiaraMaterials::CreateConcrete(){
  typedef std::map< G4Element* , G4double > TiaraMapElementFraction;

  TiaraMapElementFraction concreteFractions; 
  concreteFractions[fMapSymbolElement["H"]] = 0.010853422;
  concreteFractions[fMapSymbolElement["O"]] = 0.48165594;
  concreteFractions[fMapSymbolElement["Na"]] = 0.020327181;
  concreteFractions[fMapSymbolElement["Mg"]] = 0.010832401;
  concreteFractions[fMapSymbolElement["Al"]] = 0.060514396;
  concreteFractions[fMapSymbolElement["Si"]] = 0.22409956;
  concreteFractions[fMapSymbolElement["K"]] = 0.010680188;
  concreteFractions[fMapSymbolElement["Ca"]] = 0.12388306;
  concreteFractions[fMapSymbolElement["Fe"]] = 0.056603178;

  G4Material *mat =0;
  G4double density = 2.31*g/cm3;
  mat = new G4Material("concrete", density, 
		       concreteFractions.size());
  for (TiaraMapElementFraction::iterator it = concreteFractions.begin();
       it != concreteFractions.end(); ++it) {
    mat->AddElement(it->first , it->second);
  }
  
  return mat;

}

G4Material *TiaraMaterials::CreateAir(){
  G4Material *mat =0;
  G4double density = 1.29*mg/cm3;
  mat = new G4Material("air", density, 2);

  mat->AddElement(fMapSymbolElement["N"], 0.75);
  mat->AddElement(fMapSymbolElement["O"], 0.25);

  return mat;
}

G4Material *TiaraMaterials::CreateIron(){
  G4Material *mat =0;
  G4double density = 7.87*g/cm3;
  mat = new G4Material("iron", density, 1);
  mat->AddElement(fMapSymbolElement["Fe"],1);
  return mat;
}

G4Material *TiaraMaterials::CreateVakuum(){
  G4Material *mat =0;
  G4double density(universe_mean_density);
  G4double temperature(2.73*kelvin);
  G4double pressure(3.e-18*pascal);
  mat = new G4Material("vacuum", 1., 1.01*g/mole, 
		       density,
		       kStateGas,temperature,pressure);
  return mat;
}
