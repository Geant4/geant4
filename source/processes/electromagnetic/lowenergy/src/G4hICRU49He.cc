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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4hICRU49He
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
//
// Creation date: 20 July 2000
//
// Modifications:
// 20/07/2000  V.Ivanchenko First implementation
// 18/09/2000  V.Ivanchenko clean up - all variable are the same as in ICRU
// 03/10/2000  V.Ivanchenko clean up accoding to CodeWizard
// 10/05/2001  V.Ivanchenko Clean up againist Linux compilation with -Wall
// 26/08/2004  V.Ivanchenko Fix a problem of effective charge
//
// Class Description:
//
// Electronic stopping power parametrised according to
// ICRU Report N49, 1993. J.F. Ziegler model for He ion.
//
// Class Description: End
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hICRU49He.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hICRU49He::G4hICRU49He():G4VhElectronicStoppingPower(),
  rateMass(4.0026/1.007276),
  iMolecula(0)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hICRU49He::~G4hICRU49He()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4hICRU49He::HasMaterial(const G4Material* material)
{
  G4String chFormula = material->GetChemicalFormula() ;
  G4String myFormula = G4String(" ");

  if (myFormula == chFormula ) {
    if(1 == (material->GetNumberOfElements())) {return true;}
    return false ;
  }

  // ICRU Report N49, 1993. Power's model for He.
  const G4int numberOfMolecula = 30 ;
  static const G4String name[numberOfMolecula] = {
    "H_2", "Be-Solid", "C-Solid", "Graphite", "N_2",
    "O_2", "Al-Solid", "Si-Solid", "Ar-Solid", "Cu-Solid",
    "Ge", "W-Solid", "Au-Solid", "Pb-Solid", "C_2H_2",
    "CO_2", "Cellulose-Nitrat", "C_2H_4", "LiF",
    "CH_4", "Nylon", "Polycarbonate", "(CH_2)_N-Polyetilene", "PMMA",
    "(C_8H_8)_N", "SiO_2", "CsI", "H_2O", "H_2O-Gas"};

  // Special treatment for water in gas state

  myFormula = G4String("H_2O") ;
  const G4State theState = material->GetState() ;
  if( theState == kStateGas && myFormula == chFormula) {
    chFormula = G4String("H_2O-Gas");
  }

  // Search for the material in the table
  for (G4int i=0; i<numberOfMolecula; ++i) {
      if (chFormula == name[i]) {
        SetMoleculaNumber(i) ;
	return true ;
      }
  }
  return false ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hICRU49He::StoppingPower(const G4Material* material,
                                          G4double kineticEnergy)
{
  G4double ionloss = 0.0 ;

  // pure material (normally not the case for this function)
  if(1 == (material->GetNumberOfElements())) {
    G4double z = material->GetZ() ;
    ionloss = ElectronicStoppingPower( z, kineticEnergy ) ;

  // The data and the fit from:
  // ICRU Report N49, 1993. Power's model for He.
  } else if ( iMolecula < 30 ) {

    // Reduced kinetic energy
    // in internal units of parametrisation formula (MeV)
    G4double T = kineticEnergy*rateMass/MeV ;

    static const G4double c[30][7] = {
      {8.0080,  3.6287,  23.0700,  14.9900,  0.8507, 0.60, 2.0
   },{ 13.3100,  3.7432,  39.4130,  12.1990,  1.0950, 0.38, 1.4
   },{ 22.7240,  3.6040,  47.1810,  17.5490,  0.9040, 0.40, 1.4
   },{ 24.4040,  2.4032,  48.9440,  27.9730,  1.2933, 0.40, 1.6
   },{ 58.4719,  1.5115,  77.6421,  102.490,  1.5811, 0.50, 2.0
   },{ 60.5408,  1.6297,  91.7601,  94.1260,  1.3662, 0.50, 2.0
   },{ 48.4480,  6.4323,  59.2890,  18.3810,  0.4937, 0.48, 1.6
   },{ 59.0346,  5.1305,  47.0866,  30.0857,  0.3500, 0.60, 2.0
   },{ 71.8691,  2.8250,  51.1658,  57.1235,  0.4477, 0.60, 2.0
   },{ 78.3520,  4.0961,  136.731,  28.4470,  1.0621, 0.52, 1.2
   },{ 120.553,  1.5374,  49.8740,  82.2980,  0.8733, 0.45, 1.6
   },{ 249.896,  0.6996,  -37.274,  248.592,  1.1052, 0.50, 1.5
   },{ 246.698,  0.6219,  -58.391,  292.921,  0.8186, 0.56, 1.8
   },{ 248.563,  0.6235,  -36.8968, 306.960,  1.3214, 0.50, 2.0
   },{ 25.5860,  1.7125,  154.723,  118.620,  2.2580, 0.50, 2.0
   },{ 138.294,  25.6413, 231.873,  17.3780,  0.3218, 0.58, 1.3
   },{ 83.2091,  1.1294,  135.7457, 190.865,  2.3461, 0.50, 2.0
   },{ 263.542,  1.4754,  1541.446, 781.898,  1.9209, 0.40, 2.0
   },{ 59.5545,  1.5354,  132.1523, 153.3537, 2.0262, 0.50, 2.0
   },{ 31.7380,  19.820,  125.2100, 6.8910,   0.7242, 0.50, 1.1
   },{ 31.7549,  1.5682,  97.4777,  106.0774, 2.3204, 0.50, 2.0
   },{ 230.465,  4.8967,  1845.320, 358.641,  1.0774, 0.46, 1.2
   },{ 423.444,  5.3761,  1189.114, 319.030,  0.7652, 0.48, 1.5
   },{ 86.3410,  3.3322,  91.0433,  73.1091,  0.4650, 0.50, 2.0
   },{ 146.105,  9.4344,  515.1500, 82.8860,  0.6239, 0.55, 1.5
   },{ 238.050,  5.6901,  372.3575, 146.1835, 0.3992, 0.50, 2.0
   },{ 124.2338, 2.6730,  133.8175, 99.4109,  0.7776, 0.50, 2.0
   },{ 221.723,  1.5415,  87.7315,  192.5266, 1.0742, 0.50, 2.0
   },{ 26.7537,  1.3717,  90.8007,  77.1587,  2.3264, 0.50, 2.0
   },{ 37.6121,  1.8052,  73.0250,  66.2070,  1.4038, 0.50, 2.0} };

    G4double a1,a2 ;

  // Free electron gas model
    if ( T < 0.001 ) {
      G4double T0 = 0.001 ;
      a1 = 1.0 - G4Exp(-c[iMolecula][1]*std::pow(T0,-2.0+c[iMolecula][5])) ;
      a2 = (c[iMolecula][0]*std::log(T0)/T0 + c[iMolecula][2]/T0) *
            G4Exp(-c[iMolecula][4]*std::pow(T0,-c[iMolecula][6])) +
            c[iMolecula][3]/(T0*T0) ;

      ionloss = a1*a2*std::sqrt(T/T0) ;

  // Main parametrisation
    } else {
      a1 = 1.0 - G4Exp(-c[iMolecula][1]*std::pow(T,-2.0+c[iMolecula][5])) ;
      a2 = (c[iMolecula][0]*std::log(T)/T + c[iMolecula][2]/T) *
            G4Exp(-c[iMolecula][4]*std::pow(T,-c[iMolecula][6])) +
            c[iMolecula][3]/(T*T) ;
      ionloss = a1*a2;
    }

  // He effective charge
    G4double z = (material->GetTotNbOfElectPerVolume()) /
                 (material->GetTotNbOfAtomsPerVolume()) ;

    ionloss /= HeEffChargeSquare(z, kineticEnergy*rateMass) ;

    if ( ionloss < 0.0) ionloss = 0.0 ;
  }

  return ionloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hICRU49He::ElectronicStoppingPower(G4double z,
                                              G4double kineticEnergy) const
{
  G4double ionloss ;
  G4int i = G4int(z)-1 ;  // index of atom
  if(i < 0)  i = 0 ;
  if(i > 91) i = 91 ;

  // The data and the fit from:
  // ICRU Report 49, 1993. Ziegler's type of parametrisations
  // Reduced kinetic energy

  // He energy in internal units of parametrisation formula (MeV)
  G4double T = kineticEnergy*rateMass/MeV ;

  static const G4double a[92][5] = {
    {0.35485, 0.6456, 6.01525,  20.8933, 4.3515
   },{ 0.58,    0.59,   6.3,	 130.0,   44.07
   },{ 1.42,    0.49,   12.25,    32.0,    9.161
   },{ 2.1895,  0.47183,7.2362,   134.30,  197.96
   },{ 3.691,   0.4128, 18.48,    50.72,   9.0
   },{ 3.83523, 0.42993,12.6125,  227.41,  188.97
   },{ 1.9259,  0.5550, 27.15125, 26.0665, 6.2768
   },{ 2.81015, 0.4759, 50.0253,  10.556,  1.0382
   },{ 1.533,   0.531,  40.44,    18.41,   2.718
   },{ 2.303,   0.4861, 37.01,    37.96,   5.092
   },{ 9.894,   0.3081, 23.65,    0.384,   92.93
   },{ 4.3,     0.47,   34.3,     3.3,     12.74
   },{ 2.5,     0.625,  45.7,     0.1,     4.359
   },{ 2.1,     0.65,   49.34,    1.788,   4.133
   },{ 1.729,   0.6562, 53.41,    2.405,   3.845
   },{ 1.402,   0.6791, 58.98,    3.528,   3.211
   },{ 1.117,   0.7044, 69.69,    3.705,    2.156
   },{ 2.291,   0.6284, 73.88,    4.478,    2.066
   },{ 8.554,   0.3817, 83.61,    11.84,    1.875
   },{ 6.297,   0.4622, 65.39,    10.14,    5.036
   },{ 5.307,   0.4918, 61.74,    12.4,	   6.665
   },{ 4.71,    0.5087, 65.28,    8.806,    5.948
   },{ 6.151,   0.4524, 83.0,	 18.31,    2.71
   },{ 6.57,    0.4322, 84.76,    15.53,    2.779
   },{ 5.738,   0.4492, 84.6,	 14.18,    3.101
   },{ 5.013,   0.4707, 85.8,	 16.55,    3.211
   },{ 4.32,    0.4947, 76.14,    10.85,    5.441
   },{ 4.652,   0.4571, 80.73,    22.0,	   4.952
   },{ 3.114,   0.5236, 76.67,    7.62,	   6.385
   },{ 3.114,   0.5236, 76.67,    7.62,	   7.502
   },{ 3.114,   0.5236, 76.67,    7.62,	   8.514
   },{ 5.746,   0.4662, 79.24,    1.185,    7.993
   },{ 2.792,   0.6346, 106.1,    0.2986,   2.331
   },{ 4.667,   0.5095, 124.3,    2.102,    1.667
   },{ 2.44,    0.6346, 105.0,    0.83,	   2.851
   },{ 1.413,   0.7377, 147.9,    1.466,    1.016
   },{ 11.72,   0.3826, 102.8,    9.231,    4.371
   },{ 7.126,   0.4804, 119.3,    5.784,    2.454
   },{ 11.61,   0.3955, 146.7,    7.031,    1.423
   },{ 10.99,   0.41,   163.9,	 7.1,	   1.052
   },{ 9.241,   0.4275, 163.1,    7.954,    1.102
   },{ 9.276,   0.418,  157.1,	 8.038,    1.29
   },{ 3.999,   0.6152, 97.6,	 1.297,    5.792
   },{ 4.306,   0.5658, 97.99,    5.514,    5.754
   },{ 3.615,   0.6197, 86.26,    0.333,    8.689
   },{ 5.8,     0.49,   147.2,	 6.903,    1.289
   },{ 5.6,     0.49,   130.0,	 10.0,	   2.844
   },{ 3.55,    0.6068, 124.7,    1.112,    3.119
   },{ 3.6,     0.62,   105.8,	 0.1692,   6.026
   },{ 5.4,     0.53,   103.1,	 3.931,    7.767
   },{ 3.97,    0.6459, 131.8,    0.2233,   2.723
   },{ 3.65,    0.64,   126.8,	 0.6834,   3.411
   },{ 3.118,   0.6519, 164.9,    1.208,    1.51
   },{ 3.949,   0.6209, 200.5,    1.878,    0.9126
   },{ 14.4,    0.3923, 152.5,    8.354,    2.597
   },{ 10.99,   0.4599, 138.4,    4.811,    3.726
   },{ 16.6,    0.3773, 224.1,    6.28,	   0.9121
   },{ 10.54,   0.4533, 159.3,	 4.832,	   2.529
   },{ 10.33,   0.4502, 162.0,	 5.132,	   2.444
   },{ 10.15,   0.4471, 165.6,	 5.378,	   2.328
   },{ 9.976,   0.4439, 168.0,	 5.721,	   2.258
   },{ 9.804,   0.4408, 176.2,	 5.675,	   1.997
   },{ 14.22,   0.363,  228.4,	 7.024,	   1.016
   },{ 9.952,   0.4318, 233.5,	 5.065,	   0.9244
   },{ 9.272,   0.4345, 210.0,	 4.911,	   1.258
   },{ 10.13,   0.4146, 225.7,	 5.525,	   1.055
   },{ 8.949,   0.4304, 213.3,	 5.071,	   1.221
   },{ 11.94,   0.3783, 247.2,	 6.655,	   0.849
   },{ 8.472,   0.4405, 195.5,	 4.051,	   1.604
   },{ 8.301,   0.4399, 203.7,	 3.667,	   1.459
   },{ 6.567,   0.4858, 193.0,	 2.65,	   1.66
   },{ 5.951,   0.5016, 196.1,	 2.662,	   1.589
   },{ 7.495,   0.4523, 251.4,	 3.433,	   0.8619
   },{ 6.335,   0.4825, 255.1,	 2.834,	   0.8228
   },{ 4.314,   0.5558, 214.8,	 2.354,	   1.263
   },{ 4.02,    0.5681, 219.9,	 2.402,	   1.191
   },{ 3.836,   0.5765, 210.2,	 2.742,	   1.305
   },{ 4.68,    0.5247, 244.7,	 2.749,	   0.8962
   },{ 3.223,   0.5883, 232.7,	 2.954,	   1.05
   },{ 2.892,   0.6204, 208.6,	 2.415,	   1.416
   },{ 4.728,   0.5522, 217.0,	 3.091,	   1.386
   },{ 6.18,    0.52,   170.0,	 4.0,	   3.224
   },{ 9.0,     0.47,   198.0,	 3.8,	   2.032
   },{ 2.324,   0.6997, 216.0,	 1.599,	   1.399
   },{ 1.961,   0.7286, 223.0,	 1.621,	   1.296
   },{ 1.75,    0.7427, 350.1,	 0.9789,   0.5507
   },{ 10.31,   0.4613, 261.2,	 4.738,	   0.9899
   },{ 7.962,   0.519,  235.7,	 4.347,	   1.313
   },{ 6.227,   0.5645, 231.9,	 3.961,	   1.379
   },{ 5.246,   0.5947, 228.6,	 4.027,	   1.432
   },{ 5.408,   0.5811, 235.7,	 3.961,	   1.358
   },{ 5.218,   0.5828, 245.0,	 3.838,	   1.25}
  };

  // Free electron gas model
  if ( T < 0.001 ) {
    G4double slow  = a[i][0] ;
    G4double shigh = std::log( 1.0 + a[i][3]*1000.0 + a[i][4]*0.001 )
                   * a[i][2]*1000.0 ;
    ionloss  = slow*shigh / (slow + shigh) ;
    ionloss *= std::sqrt(T*1000.0) ;

  // Main parametrisation
  } else {
    G4double slow  = a[i][0] * std::pow((T*1000.0), a[i][1]) ;
    G4double shigh = std::log( 1.0 + a[i][3]/T + a[i][4]*T ) * a[i][2]/T ;
    ionloss = slow*shigh / (slow + shigh) ;

  }
  if ( ionloss < 0.0) ionloss = 0.0 ;

  // He effective charge
  ionloss /= HeEffChargeSquare(z, kineticEnergy*rateMass) ;

  return ionloss;
}
