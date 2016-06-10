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
// p p -> Delta Delta cross section tables
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4XDeltaDeltaTable.hh"
#include "G4PhysicsFreeVector.hh"

// Energies (GeV) corresponding to the cross section table
// Units are assigned when filling the PhysicsVector

const G4double G4XDeltaDeltaTable::energyTable[121] =
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

// Cross-sections in mb, from S.A. Bass et al., Prog.Part.Nucl.Phys.41:225-370,1998 
// Units are assigned when filling the PhysicsVector

const G4double G4XDeltaDeltaTable::sigmaDD1232[121] = 
{
  0.0, 
  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
  0.001, 0.000, 0.000, 0.000, 0.002, 0.008, 0.029, 0.078,
  0.159, 0.262, 0.374, 0.488, 0.599, 0.706, 0.806, 0.899,
  0.985, 1.064, 1.135, 1.200, 1.257, 1.309, 1.354, 1.394,
  1.429, 1.458, 1.482, 1.503, 1.519, 1.531, 1.540, 1.545,
  1.548, 1.548, 1.545, 1.540, 1.533, 1.523, 1.513, 1.500,
  1.486, 1.471, 1.455, 1.437, 1.419, 1.400, 1.381, 1.361,
  1.340, 1.319, 1.297, 1.275, 1.254, 1.231, 1.209, 1.187,
  1.180, 0.889, 0.681, 0.534, 0.430, 0.351, 0.292, 0.247,
  0.211, 0.183, 0.160, 0.141, 0.125, 0.111, 0.100, 0.090,
  0.082, 0.075, 0.068, 0.063, 0.058, 0.054, 0.050, 0.046,
  0.043, 0.040, 0.038, 0.035, 0.033, 0.031, 0.030, 0.028,
  0.027, 0.025, 0.024, 0.023, 0.022, 0.021, 0.020, 0.019
};



G4XDeltaDeltaTable::G4XDeltaDeltaTable()  : size(121)
{ }


G4XDeltaDeltaTable::~G4XDeltaDeltaTable()
{ }


G4PhysicsVector* G4XDeltaDeltaTable::CrossSectionTable() const
{
  G4PhysicsFreeVector*sigma = new G4PhysicsFreeVector(size);

  G4int i;
  for (i=0; i<size; i++)
    {
      G4double value = 0.5*sigmaDD1232[i] * millibarn;
      G4double energy = energyTable[i] * GeV;
      sigma->PutValue(i,energy,value);
    }
  return sigma;
}











