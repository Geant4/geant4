#include "B01MaterialFactory.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "PhysicalConstants.h"


B01MaterialFactory::B01MaterialFactory(){

  G4String name, symbol;

  FillElementMap("Hydrogen", "H", 1, 1.01*g/mole);
  FillElementMap("Carbon", "C", 6, 12.01*g/mole);
  FillElementMap("Oxygen", "O", 8,  16.00*g/mole);
  FillElementMap("Natrium", "Na", 11, 22.99*g/mole);
  FillElementMap("Hg", "Hg", 80, 200.59*g/mole);
  FillElementMap("Aluminium", "Al", 14, 26.98*g/mole);
  FillElementMap("Silicon", "Si", 14, 28.09*g/mole);
  FillElementMap("K", "K", 19, 39.1*g/mole);
  FillElementMap("Calzium", "Ca", 31, 69.72*g/mole);
  FillElementMap("Iron", "Fe", 26, 55.85*g/mole);

  fConcreteFractions[fMapSymbolElement["H"]] = 0.01;
  fConcreteFractions[fMapSymbolElement["O"]] = 0.529;
  fConcreteFractions[fMapSymbolElement["Na"]] = 0.016;
  fConcreteFractions[fMapSymbolElement["Hg"]] = 0.002;
  fConcreteFractions[fMapSymbolElement["Al"]] = 0.034;
  fConcreteFractions[fMapSymbolElement["Si"]] = 0.337;
  fConcreteFractions[fMapSymbolElement["K"]] = 0.013;
  fConcreteFractions[fMapSymbolElement["Ca"]] = 0.044;
  fConcreteFractions[fMapSymbolElement["Fe"]] = 0.014;
  fConcreteFractions[fMapSymbolElement["C"]] = 0.001;

}

B01MaterialFactory::~B01MaterialFactory(){
}

void B01MaterialFactory::FillElementMap(const G4String &name, 
				   const G4String &symbol,
				   G4int Z,
				   G4double A) {
  B01MapSymbolElement::iterator it = fMapSymbolElement.find(symbol);
  if (it!=fMapSymbolElement.end()) {
    G4cout << "B01MaterialFactory::FillElementMap: symbol: " 
	   << symbol << ", already defined" << G4endl;
    return;
  }
  fMapSymbolElement[symbol] = new G4Element(name, symbol, Z, A);
  return;
}

G4Material *B01MaterialFactory::CreateConcrete(){
  G4double density = 2.03*g/cm3;
  G4Material* Concrete = new G4Material("Concrete", density, 10);
  for (B01MapElementFraction::iterator it = fConcreteFractions.begin();
       it != fConcreteFractions.end(); it++) {
    Concrete->AddElement(it->first , it->second);
  }

  return Concrete;

}

G4Material *B01MaterialFactory::CreateLightConcrete(){
  G4double density = 0.0203*g/cm3;
  G4Material* LightConcrete = new G4Material("LightConcrete", density, 10);
  for (B01MapElementFraction::iterator it = fConcreteFractions.begin();
       it != fConcreteFractions.end(); it++) {
    LightConcrete->AddElement(it->first , it->second);
  }
  return LightConcrete;
}

G4Material *B01MaterialFactory::CreateGalactic(){
  G4double density = universe_mean_density;  //from PhysicalConstants.h
  G4double pressure    = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  return  new G4Material("Galactic", 1., 1.01*g/mole, density,
			 kStateGas,temperature,pressure);
}

