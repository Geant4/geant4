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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4IsotopeList
//
// Author  V. Ivantchenko, 22 October 2020
//
 
#ifndef G4IsotopeList_h
#define G4IsotopeList_h 1

static const G4int amin[] = {
   0,
   1,   3,   6,   9,  10,  12,  14,  16,  19,  20, //1-10
  22,  24,  26,  27,  31,  32,  35,  36,  39,  40, //11-20
  44,  44,  48,  50,  52,  54,  56,  56,  63,  64, //21-30
  67,  70,  71,  74,  77,  78,  85,  83,  87,  88, //31-40
  91,  92,  96,  96,  99, 102, 106, 106, 113, 112, //41-50
 121, 120, 126, 123, 133, 130, 137, 136, 141, 142, //51-60
 145, 144, 151, 148, 158, 156, 163, 162, 169, 168, //61-70
 173, 174, 179, 180, 185, 184, 190, 190, 197, 196, //71-80
 202, 204, 208, 208, 210, 222, 223, 223, 225, 227, //81-90
 231, 232, 235, 236, 241, 240, 247, 249, 253, 255};//91-100

static const G4int amax[] = {
   0,
   3,   4,   7,   9,  11,  14,  15,  18,  19,  22, //1-10
  23,  27,  27,  32,  33,  36,  37,  41,  41,  48, //11-20
  48,  50,  51,  54,  55,  60,  62,  66,  67,  70, //21-30
  71,  76,  77,  82,  82,  86,  88,  90,  91,  96, //31-40
  95, 100,  99, 106, 105, 110, 111, 116, 115, 126, //41-50
 127, 132, 135, 136, 137, 140, 140, 144, 143, 150, //51-60
 151, 154, 157, 161, 160, 165, 166, 172, 171, 176, //61-70
 177, 182, 182, 188, 188, 193, 193, 198, 197, 204, //71-80
 205, 208, 210, 209, 210, 222, 223, 226, 227, 234, //81-90
 233, 241, 239, 246, 244, 250, 250, 254, 255, 255};//91-100

static const G4double aeff[] = {
 0.,
 1.00794, 4.00264, 6.94003, 9.01218,  10.811, 12.0107, 14.0068, 15.9994, 18.9984,   20.18,  //1-10
 22.9898,  24.305, 26.9815, 28.0854, 30.9738, 32.0661, 35.4526, 39.9477, 39.0983,  40.078,  //11-20
 44.9559, 47.8667, 50.9415, 51.9961,  54.938, 55.8451, 58.9332, 58.6933, 63.5456, 65.3955,  //21-30
 69.7231, 72.6128, 74.9216, 78.9594, 79.9035, 83.7993, 85.4677, 87.6166, 88.9058, 91.2236,  //31-40
 92.9064, 95.9313, 97.9072, 101.065, 102.906, 106.415, 107.868, 112.411, 114.818,  118.71,  //41-50
  121.76, 127.603, 126.904, 131.292, 132.905, 137.327, 138.905, 140.115, 140.908, 144.236,  //51-60
 144.913, 150.366, 151.964, 157.252, 158.925, 162.497,  164.93, 167.256, 168.934, 173.038,  //61-70
 174.967, 178.485, 180.948, 183.842, 186.207, 190.225, 192.216, 195.078, 196.967, 200.599,  //71-80
 204.383, 207.217,  208.98, 208.982, 209.987, 222.018,  223.02, 226.025, 227.028, 232.038,  //81-90
 231.036, 238.029, 237.048, 244.064, 243.061,  247.07,  247.07,  251.08, 252.083, 257.095}; //91-100

static const G4String elementName[] = {
  "",
  "Hydrogen",     "Helium",      "Lithium",     "Berylium",   "Boron",        "Carbon",
  "Nitrogen",     "Oxygen",      "Fluorine",    "Neon",       "Sodium",       "Magnesium",
  "Aluminum",     "Silicon",     "Phosphorous", "Sulfur",     "Chlorine",     "Argon",
  "Potassium",    "Calcium",     "Scandium",    "Titanium",   "Vanadium",     "Chromium",
  "Manganese",    "Iron",        "Cobalt",      "Nickel",     "Copper",       "Zinc",
  "Gallium",      "Germanium",   "Arsenic",     "Selenium",   "Bromine",      "Krypton",
  "Rubidium",     "Strontium",   "Yttrium",     "Zirconium",  "Niobium",      "Molybdenum",
  "Technetium",   "Ruthenium",   "Rhodium",     "Palladium",  "Silver",       "Cadmium",
  "Indium",       "Tin",         "Antimony",    "Tellurium",  "Iodine",       "Xenon",
  "Cesium",       "Barium",      "Lanthanum",   "Cerium",     "Praseodymium", "Neodymium",
  "Promethium",   "Samarium",    "Europium",    "Gadolinium", "Terbium",      "Dysprosium",
  "Holmium",      "Erbium",      "Thulium",     "Ytterbium",  "Lutetium",     "Hafnium",
  "Tantalum",     "Tungsten",    "Rhenium",     "Osmium",     "Iridium",      "Platinium",
  "Gold",         "Mercury",     "Thallium",    "Lead",       "Bismuth",      "Polonium",
  "Astatine",     "Radon",       "Francium",    "Radium",     "Actinium",     "Thorium",
  "Protactinium", "Uranium",     "Neptunium",   "Plutonium",  "Americium",    "Curium",
  "Berkelium",    "Californium", "Einsteinium", "Fermium"};

#endif
