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
// $Id: Tst33MaterialFactory.cc,v 1.3 2002-10-31 08:32:44 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33MaterialFactory.cc
//
// ----------------------------------------------------------------------

#include "Tst33MaterialFactory.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "PhysicalConstants.h"


Tst33MaterialFactory::Tst33MaterialFactory(){

  G4String name, symbol;

  FillElementMap("Hydrogen", "H", 1, 1.01*G4std::g/G4std::mole);
  FillElementMap("Carbon", "C", 6, 12.01*G4std::g/G4std::mole);
  FillElementMap("Oxygen", "O", 8,  16.00*G4std::g/G4std::mole);
  FillElementMap("Natrium", "Na", 11, 22.99*G4std::g/G4std::mole);
  //  FillElementMap("Hg", "Hg", 80, 200.59*G4std::g/G4std::mole);
  FillElementMap("Aluminium", "Al", 14, 26.98*G4std::g/G4std::mole);
  FillElementMap("Silicon", "Si", 14, 28.09*G4std::g/G4std::mole);
  FillElementMap("K", "K", 19, 39.1*G4std::g/G4std::mole);
  FillElementMap("Calzium", "Ca", 31, 69.72*G4std::g/G4std::mole);
  FillElementMap("Iron", "Fe", 26, 55.85*G4std::g/G4std::mole);

  fConcreteFractions[fMapSymbolElement["H"]] = 0.01;
  fConcreteFractions[fMapSymbolElement["O"]] = 0.529;
  fConcreteFractions[fMapSymbolElement["Na"]] = 0.016;
  //  fConcreteFractions[fMapSymbolElement["Hg"]] = 0.002 ;
  fConcreteFractions[fMapSymbolElement["Al"]] = 0.034;
  fConcreteFractions[fMapSymbolElement["Si"]] = 0.337 + 0.002;
  fConcreteFractions[fMapSymbolElement["K"]] = 0.013;
  fConcreteFractions[fMapSymbolElement["Ca"]] = 0.044;
  fConcreteFractions[fMapSymbolElement["Fe"]] = 0.014;
  fConcreteFractions[fMapSymbolElement["C"]] = 0.001;

}

Tst33MaterialFactory::~Tst33MaterialFactory(){
}

void Tst33MaterialFactory::FillElementMap(const G4String &name, 
				   const G4String &symbol,
				   G4int Z,
				   G4double A) {
  Tst33MapSymbolElement::iterator it = fMapSymbolElement.find(symbol);
  if (it!=fMapSymbolElement.end()) {
    G4std::G4cout << "Tst33MaterialFactory::FillElementMap: symbol: " 
	   << symbol << ", already defined" << G4endl;
  }
  else {
    fMapSymbolElement[symbol] = new G4Element(name, symbol, Z, A);
    if (!fMapSymbolElement[symbol]) {
      G4std::G4Exception("Tst33MaterialFactory::FillElementMap: new failed to create G4Element!");
    }
  }
  return;
}

G4Material *Tst33MaterialFactory::CreateConcrete(){
  G4double density = 2.03*G4std::g/G4std::cm3;
  G4Material* Concrete = 0;
  Concrete = new G4Material("Concrete", density, 
			    fConcreteFractions.size());
  for (Tst33MapElementFraction::iterator it = fConcreteFractions.begin();
       it != fConcreteFractions.end(); ++it) {
    Concrete->AddElement(it->first , it->second);
  }

  return Concrete;

}

G4Material *Tst33MaterialFactory::CreateLightConcrete(){
  G4double density = 0.0203*G4std::g/G4std::cm3;
  G4Material* LightConcrete = 0;
  LightConcrete = new G4Material("LightConcrete", density, 
				 fConcreteFractions.size());
  for (Tst33MapElementFraction::iterator it = fConcreteFractions.begin();
       it != fConcreteFractions.end(); ++it) {
    LightConcrete->AddElement(it->first , it->second);
  }
  return LightConcrete;
}

G4Material *Tst33MaterialFactory::CreateGalactic(){
  G4double density = G4std::universe_mean_density;  //from PhysicalConstants.h
  G4double pressure    = 3.e-18*pascal;
  G4double temperature = 2.73*G4std::kelvin;
  return  new G4Material("Galactic", 1., 1.01*G4std::g/G4std::mole, density,
			 kStateGas,temperature,pressure);
}

