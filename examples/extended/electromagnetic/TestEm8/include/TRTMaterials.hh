#ifndef TRT_MATERIALS_HH
#define TRT_MATERIALS_HH

//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel, nComponents;

  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen", symbol="O", iz=8., a);
  a = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon", symbol="C", iz=6., a);
  a = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen", symbol="H", iz=1., a);
  a = 131.29*g/mole;
  G4Element* elXe = new G4Element(name="Xenon", symbol="Xe", iz=54., a);  
  a = 19.00*g/mole;
  G4Element* elF  = new G4Element(name="Fluorine", symbol="F", iz=9., a);

  density = 1.205*mg/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2,
					kStateGas,293.15*kelvin,1.*atmosphere);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);
  
  density = 1.39*g/cm3;
  G4Material* Kapton = new G4Material(name="Kapton", density, nel=3);
  Kapton->AddElement(elO,2);
  Kapton->AddElement(elC,5);
  Kapton->AddElement(elH,4);
  
  a = 63.55*g/mole;
  density = 8.96*g/cm3;
  G4Material* Copper = new G4Material(name="Copper", z=29.,a, density);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminum", z=13., a, density);

  a = 28.09*g/mole;
  density = 2.33*g/cm3;
  G4Material* Si = new G4Material(name="Silicon", z=14., a, density);

  density = 1.977*mg/cm3;
  G4Material* CO2 = new G4Material(name="CO2", density, nel=2,
				       kStateGas,273.15*kelvin,1.*atmosphere);
  CO2->AddElement(elC,1);
  CO2->AddElement(elO,2);

  G4double TRT_Xe_density = 5.485*mg/cm3;
  G4Material* TRT_Xe = new G4Material(name="TRT_Xe", TRT_Xe_density, nel=1,
				      kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_Xe->AddElement(elXe,1);

  G4double TRT_CO2_density = 1.842*mg/cm3;
  G4Material* TRT_CO2 = new G4Material(name="TRT_CO2", TRT_CO2_density, nel=2,
				       kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_CO2->AddElement(elC,1);
  TRT_CO2->AddElement(elO,2);

  G4double TRT_CF4_density = 3.9*mg/cm3;
  G4Material* TRT_CF4 = new G4Material(name="TRT_CF4", TRT_CF4_density, nel=2,
                                       kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_CF4->AddElement(elC,1);
  TRT_CF4->AddElement(elF,4);

  G4double XeCO2CF4_density = 4.76*mg/cm3;
  G4Material* XeCO2CF4 = new G4Material(name="XeCO2CF4", XeCO2CF4_density,
					nComponents=3,
					kStateGas,293.15*kelvin,1.*atmosphere);
  XeCO2CF4->AddMaterial(TRT_Xe,0.807);
  XeCO2CF4->AddMaterial(TRT_CO2,0.039);
  XeCO2CF4->AddMaterial(TRT_CF4,0.154);
      
  density = 0.935*g/cm3;
  G4Material* TRT_CH2 = new G4Material(name="TRT_CH2",density, nel=2);
  TRT_CH2->AddElement(elC,1);
  TRT_CH2->AddElement(elH,2);

  density = 0.059*g/cm3;
  G4Material* Radiator = new G4Material(name="Radiator",density, nel=2);
  Radiator->AddElement(elC,1);
  Radiator->AddElement(elH,2);

  density = 0.145*g/cm3;
  G4Material* CarbonFiber = new G4Material(name="CarbonFiber",density, nel=1);
  CarbonFiber->AddElement(elC,1);

#endif

