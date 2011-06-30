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

#include "Materials.hh"
#include "G4Material.hh"


Materials* Materials::instance = 0;


Materials* Materials::Instance () {

  if(instance == 0) {
     instance = new Materials;
  }
  return instance;
}


void Materials::Destroy() {

  if(!instance == 0) {

     delete instance;
     instance = 0;
  }
}


Materials::Materials() {

  G4String name, chemSymbol;         
  G4double moleMass, Z, density;  

  G4int ncomponents, isotopeN, natoms;
  G4double abundance, fractionmass;
  G4double temperature, pressure;

  moleMass = 1.01 * g/mole;
  G4Element* H  = new G4Element(name="Hydrogen",chemSymbol="H" , Z = 1., moleMass);

  moleMass = 14.01 * g/mole;
  G4Element* N  = new G4Element(name="Nitrogen",chemSymbol="N" , Z = 7., moleMass);

  moleMass = 16.00 * g/mole;
  G4Element* O  = new G4Element(name="Oxygen"  ,chemSymbol="O" , Z = 8., moleMass);

  moleMass = 28.09 * g/mole;
  G4Element* Si = new G4Element(name="Silicon",chemSymbol="Si" , Z = 14., moleMass);

  moleMass = 235.01 * g/mole;
  G4Isotope* U5 = new G4Isotope(name="U235", 92, isotopeN = 235, moleMass);

  moleMass = 238.03 * g/mole;
  G4Isotope* U8 = new G4Isotope(name="U238", 92, isotopeN = 238, moleMass);

  G4Element* U  = new G4Element(name="enriched Uranium",chemSymbol="U",ncomponents=2);
  U->AddIsotope(U5, abundance= 90.*perCent);
  U->AddIsotope(U8, abundance= 10.*perCent);

  moleMass = 1.01 * g/mole;
  density = 8.3748e-5 * g/cm3; 
  hydrogen = new G4Material(name="Hydrogen", Z = 1., moleMass, density);
  hydrogen -> GetIonisation()-> SetMeanExcitationEnergy(19.2 * eV);
 
  moleMass = 9.012 * g/mole;
  density = 1.848 * g/cm3; 
  beryllium = new G4Material(name="Beryllium", Z = 4., moleMass, density);
  beryllium -> GetIonisation()-> SetMeanExcitationEnergy(63.7 * eV);

  moleMass = 12.01 * g/mole;
  density = 1.7 * g/cm3; 
  graphite = new G4Material(name="Graphite", Z = 6., moleMass, density);
  graphite -> GetIonisation()->SetMeanExcitationEnergy(78. * eV);
 
  moleMass = 16.00 * g/mole;
  density = 0.00133 * g/cm3;
  ossigeno = new G4Material(name="Oxygen", Z = 8.,moleMass,density);
  ossigeno->GetIonisation()->SetMeanExcitationEnergy(95. * eV);

  moleMass = 24.312 * g/mole;
  density = 1.738 * g/cm3; 
  magnesium = new G4Material(name="Magnesium", Z = 12., moleMass, density);
  magnesium -> GetIonisation()-> SetMeanExcitationEnergy(156. * eV);

  moleMass = 26.981 * g/mole;
  density = 2.6989 * g/cm3; 
  aluminium = new G4Material(name="Aluminium", Z = 13., moleMass, density);
  aluminium->GetIonisation()->SetMeanExcitationEnergy(166.0 * eV);

  moleMass = 28.085 * g/mole;
  density = 2.336 * g/cm3; 
  silicon = new G4Material(name="Silicon", Z = 14., moleMass, density);
  silicon->GetIonisation()->SetMeanExcitationEnergy(173.0 * eV);

  moleMass = 39.95 * g/mole;
  density = 1.390 * g/cm3;
  liquidArgon = new G4Material(name="liquidArgon", Z = 18., moleMass, density);
  liquidArgon ->GetIonisation()->SetMeanExcitationEnergy(188.0 * eV);

  moleMass = 47.867 * g/mole;
  density = 4.54 * g/cm3;
  titanium = new G4Material(name="Titanium", Z = 22., moleMass, density);
  titanium -> GetIonisation()->SetMeanExcitationEnergy(233.0 * eV);
 
  moleMass = 55.845 * g/mole;
  density = 7.874 * g/cm3;
  iron = new G4Material(name="Iron", Z = 26., moleMass, density);
  iron -> GetIonisation()->SetMeanExcitationEnergy(286. * eV);
 
  moleMass = 58.933 * g/mole;
  density = 8.9 * g/cm3;
  cobalt = new G4Material(name="Cobalt", Z = 27., moleMass, density);
  cobalt -> GetIonisation()->SetMeanExcitationEnergy(297. * eV);

  moleMass = 58.69 * g/mole;
  density = 8.902 * g/cm3;
  nickel = new G4Material(name="Nickel", Z = 28., moleMass, density);
  nickel -> GetIonisation()->SetMeanExcitationEnergy(311. * eV);

  moleMass = 63.546 * g/mole;
  density = 8.96 * g/cm3;
  copper = new G4Material(name="Copper", Z = 29., moleMass, density);
  copper -> GetIonisation()->SetMeanExcitationEnergy(322. * eV);
 
  moleMass = 65.409 * g/mole;
  density = 7.133 * g/cm3;
  zinc = new G4Material(name="Zinc", Z = 30., moleMass, density);
  zinc -> GetIonisation()->SetMeanExcitationEnergy(330. * eV);
 
  moleMass = 69.723 * g/mole;
  density = 5.904 * g/cm3;
  gallium = new G4Material(name="Gallium", Z = 31., moleMass, density);
  gallium -> GetIonisation()->SetMeanExcitationEnergy(334. * eV);

  moleMass = 72.64 * g/mole;
  density = 5.323 * g/cm3;
  germanium = new G4Material(name="Germanium", Z = 32., moleMass, density);
  germanium -> GetIonisation()->SetMeanExcitationEnergy(350. * eV);

  moleMass = 91.224  * g/mole;
  density = 6.506 * g/cm3;
  zirconium = new G4Material(name="Zirconium", Z = 40., moleMass, density);
  zirconium -> GetIonisation()->SetMeanExcitationEnergy(393. * eV);
 
  moleMass = 95.94  * g/mole;
  density = 10.22 * g/cm3;
  molybdenium = new G4Material(name="Molybdenium", Z = 42., moleMass, density);
  molybdenium -> GetIonisation()->SetMeanExcitationEnergy(424. * eV);

  moleMass = 107.8682 * g/mole;
  density = 10.5 * g/cm3;
  silver = new G4Material(name="Silver", Z = 47., moleMass, density);
  silver -> GetIonisation()->SetMeanExcitationEnergy(470. * eV);

  moleMass = 112.411 * g/mole;
  density = 8.65 * g/cm3;
  cadmium = new G4Material(name="Cadmium", Z = 48., moleMass, density);
  cadmium -> GetIonisation()->SetMeanExcitationEnergy(469. * eV);

  moleMass = 114.818 * g/mole;
  density = 7.310 * g/cm3;
  indium = new G4Material(name="Indium", Z = 49., moleMass, density);
  indium -> GetIonisation()->SetMeanExcitationEnergy(488. * eV);

  moleMass = 118.71 * g/mole;
  density = 7.310 * g/cm3;
  tin = new G4Material(name="Tin", Z = 50., moleMass, density);
  tin -> GetIonisation()->SetMeanExcitationEnergy(488. * eV);

  moleMass = 132.90545 * g/mole;
  density = 1.873 * g/cm3;
  cesium = new G4Material(name="Cesium", Z = 55., moleMass, density);
  cesium -> GetIonisation()->SetMeanExcitationEnergy(488. * eV);

  moleMass = 150.36 * g/mole;
  density = 7.46 * g/cm3;
  samarium = new G4Material(name="Samarium", Z = 62., moleMass, density);
  samarium -> GetIonisation()->SetMeanExcitationEnergy(574. * eV);

  moleMass = 173.04 * g/mole;
  density = 6.73 * g/cm3;
  ytterbium = new G4Material(name="Ytterbium", Z = 70., moleMass, density);
  ytterbium -> GetIonisation()->SetMeanExcitationEnergy(684. * eV);

  moleMass = 180.9947 * g/mole;
  density = 16.65 * g/cm3;
  tantalum = new G4Material(name="Tantalum", Z = 73., moleMass, density);
  tantalum -> GetIonisation()->SetMeanExcitationEnergy(718. * eV);

  moleMass = 183.85 * g/mole;
  density = 19.3 * g/cm3;
  tungsten = new G4Material(name="Tungsten", Z = 74., moleMass, density);
  tungsten -> GetIonisation()->SetMeanExcitationEnergy(727. * eV);

  moleMass = 196.966 * g/mole;
  density = 19.32 * g/cm3;
  gold = new G4Material(name="Gold", Z = 79., moleMass, density);
  gold -> GetIonisation()->SetMeanExcitationEnergy(790. * eV);

  moleMass = 207.19 * g/mole;
  density = 11.35 * g/cm3;
  lead = new G4Material(name="Lead", Z = 82., moleMass, density);
  lead -> GetIonisation()->SetMeanExcitationEnergy(823. * eV);

  moleMass = 238. * g/mole;
  density = 18.95 * g/cm3;
  uranium = new G4Material(name="Uranium", Z = 92., moleMass, density);
  uranium -> GetIonisation()->SetMeanExcitationEnergy(890. * eV);
 
  density = 1.000 * g/cm3;
  water = new G4Material(name="Water", density, ncomponents=2);
  //water->SetChemicalFormula("H_2O");
  water->AddElement(H, natoms=2);
  water->AddElement(O, natoms=1);
  water->GetIonisation()->SetMeanExcitationEnergy(75.0 * eV);

  density = 2.200 * g/cm3;
  quartz = new G4Material(name="Quartz", density, ncomponents=2);
  quartz->AddElement(Si, natoms=1);
  quartz->AddElement(O , natoms=2);

  density = 1.290 * mg/cm3;
  air = new G4Material(name="Air", density, ncomponents=2);
  air->AddElement(N, fractionmass=0.7);
  air->AddElement(O, fractionmass=0.3);
  air->GetIonisation()->SetMeanExcitationEnergy(85.7 * eV);

  density = universe_mean_density;
  moleMass = 1.01 * g/mole;
  pressure = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  vacuum = new G4Material(name="Galactic", Z = 1., moleMass,
				      density,kStateGas,temperature,pressure);

  density = 1.16*mg/cm3;
  moleMass = 14. * g/mole;
  nytrogen = new G4Material(name="Nytrogen", Z = 7., moleMass, density);
  nytrogen->GetIonisation()->SetMeanExcitationEnergy(82.0 * eV);
}


Materials::~Materials() {

  delete nytrogen;
  delete vacuum;
  delete air;
  delete quartz;
  delete water;
  delete uranium; 
  delete lead;
  delete gold;
  delete tungsten;
  delete tantalum;
  delete ytterbium;
  delete samarium;
  delete cesium;
  delete tin;
  delete indium;
  delete cadmium;
  delete silver;
  delete molybdenium;
  delete zirconium;
  delete germanium;
  delete gallium;
  delete zinc;
  delete copper;
  delete nickel;
  delete cobalt;
  delete iron; 
  delete titanium;
  delete liquidArgon;
  delete silicon;
  delete aluminium;
  delete magnesium;
  delete graphite;
  delete beryllium;
  delete hydrogen;
}


G4Material* Materials::GetMaterial(G4String matName) {

  return G4Material::GetMaterial(matName);
}

