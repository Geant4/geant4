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

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4XNDeltaTable.hh"
#include "G4PhysicsFreeVector.hh"

// Energies (GeV) corresponding to the cross section table
// Units are assigned when filling the PhysicsVector

const G4double G4XNDeltaTable::energyTable[121] =
{
  0.0, 
  2.014,  2.014,  2.016,  2.018,  2.022,  2.026,  2.031,  2.037,  2.044,  2.052,   
  2.061,  2.071,  2.082,  2.094,  2.107,  2.121,  2.135,  2.151,  2.168,  2.185,   
  2.204,  2.223,  2.244,  2.265,  2.287,  2.311,  2.335,  2.360,  2.386,  2.413,   
  2.441,  2.470,  2.500,  2.531,  2.562,  2.595,  2.629,  2.664,  2.699,  2.736,
  2.773,  2.812,  2.851,  2.891,  2.933,  2.975,  3.018,  3.062,  3.107,  3.153,   
  3.200,  3.248,  3.297,  3.347,  3.397,  3.449,  3.502,  3.555,  3.610,  3.666,   
  3.722,  3.779,  3.838,  3.897,  3.957,  4.018,  4.081,  4.144,  4.208,  4.273,   
  4.339,  4.406,  4.473,  4.542,  4.612,  4.683,  4.754,  4.827,  4.900,  4.975,
  5.000,  6.134,  7.269,  8.403,  9.538, 10.672, 11.807, 12.941, 14.076, 15.210,  
 16.345, 17.479, 18.613, 19.748, 20.882, 22.017, 23.151, 24.286, 25.420, 26.555, 
 27.689, 28.824, 29.958, 31.092, 32.227, 33.361, 34.496, 35.630, 36.765, 37.899,  
 39.034, 40.168, 41.303, 42.437, 43.571, 44.706, 45.840, 46.975, 48.109, 49.244
};

const G4double G4XNDeltaTable::sigmaND1232[121] = 
{
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.86797, 4.66734, 7.69136, 10.6795, 13.4186, 16.1709, 18.589, 20.8729, 22.7591, 24.3304, 25.6082, 26.4774, 27.0037, 27.2032, 27.102, 26.7335, 26.1361, 25.35, 24.4149, 23.4022, 22.2798, 21.1099, 19.9181, 18.7585, 17.5815, 16.4636, 15.3542, 14.3173, 13.3284, 12.3683, 11.4842, 10.6528, 9.87343, 9.14494, 8.46562, 7.83342, 7.24607, 6.70118, 6.19628, 5.73764, 5.30454, 4.90425, 4.54101, 4.19898, 3.88329, 3.59684, 3.33204, 3.08328, 2.85746, 2.64878, 2.45594, 2.27499, 2.11055, 1.95855, 1.81805, 1.68815, 1.56803, 1.4585, 1.3556, 1.26038, 1.17227, 1.09178, 1.01616, 0.94702, 0.882049, 0.861647, 0.338175, 0.159097, 0.084589, 0.0490573, 0.030404, 0.0198328, 0.0134898, 0.00949051, 0.00687213, 0.00509677, 0.00386061, 0.00297765, 0.0023329, 0.00185401, 0.00149168, 0.00121394, 0.000997782, 0.000827808, 0.000692435, 0.000583729, 0.000495476, 0.000423349, 0.000363883, 0.000314463, 0.00027319, 0.000238444, 0.000209075, 0.000184071, 0.000162714, 0.000144352, 0.000128522, 0.000114792, 0.000102858, 9.24351e-05, 8.32923e-05, 7.52591e-05, 6.81657e-05, 6.18939e-05, 5.63225e-05, 
};



G4XNDeltaTable::G4XNDeltaTable() : size(121)
{ }


G4XNDeltaTable::~G4XNDeltaTable()
{ }


G4PhysicsVector* G4XNDeltaTable::CrossSectionTable() const
{
  G4PhysicsFreeVector* sigma = new G4PhysicsFreeVector(size);

  G4int i;
  for (i=0; i<size; i++)
    {
      G4double value = 0.5*sigmaND1232[i] * millibarn;
      G4double energy = energyTable[i] * GeV;
      sigma->PutValue(i,energy,value);
    }
  return sigma;
}
