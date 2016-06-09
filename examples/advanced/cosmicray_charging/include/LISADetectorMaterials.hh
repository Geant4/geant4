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
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************


void LISADetectorConstruction::ConstructMaterials(void) {

  // Elements

  G4Element* elementH  = new G4Element("Hydrogen",  "H",  1.,  1.0079*g/mole);
  G4Element* elementC  = new G4Element("Carbon",    "C",  6.,  12.011*g/mole);
  G4Element* elementN  = new G4Element("Nitrogen",  "N",  7.,  14.007*g/mole);
  G4Element* elementO  = new G4Element("Oxygen",    "O",  8., 15.9994*g/mole);
  G4Element* elementMg = new G4Element("Magnesium", "Mg",12., 24.3050*g/mole);
  G4Element* elementAl = new G4Element("Aluminium", "Al",13., 26.9815*g/mole);
  G4Element* elementSi = new G4Element("Silicon",   "Si",14., 28.0855*g/mole);
  G4Element* elementTi = new G4Element("Titanium",  "Ti",22.,   47.90*g/mole);
  G4Element* elementV  = new G4Element("Vanadium",  "V", 23., 50.9415*g/mole);
  G4Element* elementFe = new G4Element("Iron",      "Fe",26.,  55.845*g/mole);
  G4Element* elementMo = new G4Element("Molybdenum","Mo",42.,   95.94*g/mole);
  G4Element* elementPt = new G4Element("Platinum",  "Pt",78.,  195.08*g/mole);
  G4Element* elementAu = new G4Element("Gold",      "Au",79.,  196.97*g/mole);


  // Materials

  vacuum = new G4Material("vacuum", 1, 1.00794*g/mole, 
     1.0E-25*g/cm3, kStateGas, 0.1*kelvin, 1.0E-19*pascal);

  // Aluminium alloy 6061-T6
  Al6061 = new G4Material("Al6061",          2.70*g/cm3, 4);
  Al6061->AddElement(elementAl, 0.980);
  Al6061->AddElement(elementMg, 0.010);
  Al6061->AddElement(elementSi, 0.006);
  Al6061->AddElement(elementFe, 0.004);

  // Aluminium Alloy honeycomb
  AlHoneycomb = new G4Material("AlHoneycomb",0.05*g/cm3, 4);
  AlHoneycomb->AddElement(elementAl, 0.980);
  AlHoneycomb->AddElement(elementMg, 0.010);
  AlHoneycomb->AddElement(elementSi, 0.006);
  AlHoneycomb->AddElement(elementFe, 0.004);

  // Solar Cells****
  Scell = new G4Material("Scell",            7.82*g/cm3, 1);
  Scell->AddElement(elementSi, 1.00);

  // MLI Blanket: mylar (composition from STAR Stopping Power database)
  MLImat = new G4Material("MLImat",          1.40*g/cm3, 3);
  MLImat->AddElement(elementH, 0.041959);
  MLImat->AddElement(elementC, 0.625017);
  MLImat->AddElement(elementO, 0.333025);

  // Molybdenum
  molybdenum = new G4Material("molybdenum", 10.22*g/cm3, 1);
  molybdenum->AddElement(elementMo, 1.00);

  // Titanium Alloy (Ti-6Al-4V)
  TiAlloy   = new G4Material("TiAlloy",    4.43*g/cm3, 3);
  TiAlloy->AddElement(elementTi, 0.90);
  TiAlloy->AddElement(elementAl, 0.06);
  TiAlloy->AddElement(elementV,  0.04);

  // Gold
  gold       = new G4Material("gold",       19.32*g/cm3, 1);
  gold->AddElement(elementAu, 1.00);

  // Gold-Platinum 70%/30%
  AuPt = new G4Material("AuPt",             19.92*g/cm3, 2);
  AuPt->AddElement(elementAu, 0.70);
  AuPt->AddElement(elementPt, 0.30);

  // CRFP (Carbon Fiber Reinforced Polymer): M55 Quasiisotropic Layup
  CFRP = new G4Material("CFRP",               1.66*g/cm3, 1);
  CFRP->AddElement(elementC,1);

  // Silicate glass
  G4Material* SiGlass = new G4Material("SiGlass", 2.20*g/cm3,  2);
  SiGlass->AddElement(elementO, 2);
  SiGlass->AddElement(elementSi,1);
  // Titanium glass
  G4Material* TiGlass = new G4Material("TiGlass", 4.25*g/cm3,  2);
  TiGlass->AddElement(elementO, 2);
  TiGlass->AddElement(elementTi,1);
  // Corning 7972 ULE Titanium Silicate Glass
  ULEglass = new G4Material("ULEglass",      2.21*g/cm3, 2);
  ULEglass->AddMaterial(SiGlass, 0.925);
  ULEglass->AddMaterial(TiGlass, 0.075);

  // SHAPAL-M (AlN)
  SHAPAL = new G4Material("SHAPAL-M",        2.90*g/cm3, 2);
  SHAPAL->AddElement(elementAl,1);
  SHAPAL->AddElement(elementN, 1);

  // Silicon carbide
  SiC = new G4Material("SiC",                 3.1*g/cm3, 2);
  SiC->AddElement(elementSi,1);
  SiC->AddElement(elementC, 1);

  // Foam: Polystyrene-based
  foam = new G4Material("foam",0.05*g/cm3, 2);
  foam->AddElement(elementC, 0.90);
  foam->AddElement(elementH, 0.10);


}
