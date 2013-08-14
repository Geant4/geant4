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
// $Id: G4GamP2NPipAngDst.cc 67433 2013-02-20 21:27:31Z mkelsey $
// Author:  Dennis Wright (SLAC)
// Date:    28 January 2013
//
// Description: implementation of numerically integrated angular distribution
//              class for gamma p -> n pi+ reaction
//
// 20130219	Inherit from templated base, move arrays to namespace statics

#include "G4GamP2NPipAngDst.hh"

namespace {
  static const G4double eBins[15] =
    { 0.1515, 0.185, 0.260, 0.350, 0.400, 0.450, 0.603,
      0.698, 0.793, 0.902, 1.056, 1.162, 1.269, 1.480,
      1.770 };

  static const G4double angleBins[19] = 
    { -1.000, -0.985, -0.940, -0.866, -0.766, -0.643, -0.500, -0.342, -0.174,
       0.000,  0.174,  0.342,  0.500,  0.643,  0.766,  0.866,  0.940,  0.985,
       1.000 };

  static const G4double integralTable[15][19] = {
    {0.0000, 0.0075, 0.0300, 0.0670, 0.1170, 0.1785, 0.2500, 0.3290, 0.4130,
     0.5000, 0.5870, 0.6710, 0.7500, 0.8215, 0.8830, 0.9330, 0.9700, 0.9925,
     1.0000}, 
    {0.0000, 0.0090, 0.0357, 0.0792, 0.1372, 0.2071, 0.2866, 0.3721, 0.4604,
     0.5488, 0.6344, 0.7141, 0.7866, 0.8501, 0.9032, 0.9452, 0.9757, 0.9940,
     1.0000},
    {0.0000, 0.0079, 0.0313, 0.0694, 0.1215, 0.1881, 0.2705, 0.3680, 0.4767,
     0.5898, 0.6974, 0.7901, 0.8633, 0.9162, 0.9516, 0.9746, 0.9891, 0.9973,
     1.0000},
    {0.0000, 0.0061, 0.0245, 0.0550, 0.0974, 0.1516, 0.2172, 0.2927, 0.3767,
     0.4688, 0.5660, 0.6624, 0.7510, 0.8257, 0.8849, 0.9314, 0.9674, 0.9915,
     1.0000}, 
    {0.0000, 0.0044, 0.0171, 0.0379, 0.0675, 0.1081, 0.1608, 0.2246, 0.2982,
     0.3824, 0.4770, 0.5778, 0.6776, 0.7674, 0.8428, 0.9048, 0.9544, 0.9880,
     1.0000},
    {0.0000, 0.0041, 0.0162, 0.0362, 0.0647, 0.1033, 0.1536, 0.2148, 0.2859,
     0.3677, 0.4605, 0.5614, 0.6634, 0.7567, 0.8351, 0.8996, 0.9514, 0.9871,
     1.0000}, 
    {0.0000, 0.0027, 0.0118, 0.0282, 0.0512, 0.0809, 0.1191, 0.1678, 0.2268,
     0.2967, 0.3810, 0.4821, 0.5950, 0.7048, 0.7981, 0.8741, 0.9370, 0.9827,
     1.0000},
    {0.0000, 0.0021, 0.0101, 0.0262, 0.0499, 0.0808, 0.1226, 0.1789, 0.2489,
     0.3289, 0.4186, 0.5208, 0.6328, 0.7403, 0.8285, 0.8958, 0.9483, 0.9858,
     1.0000}, 
    {0.0000, 0.0036, 0.0151, 0.0348, 0.0622, 0.0996, 0.1524, 0.2224, 0.3024,
     0.3827, 0.4607, 0.5428, 0.6338, 0.7270, 0.8107, 0.8814, 0.9406, 0.9837,
     1.0000},
    {0.0000, 0.0025, 0.0156, 0.0411, 0.0711, 0.1105, 0.1676, 0.2283, 0.2776,
    0.3269, 0.3870, 0.4513, 0.5284, 0.6381, 0.7618, 0.8593, 0.9260, 0.9774,
     1.0000},
    {0.0000, 0.0015, 0.0114, 0.0326, 0.0603, 0.0994, 0.1528, 0.2003, 0.2283,
     0.2591, 0.3160, 0.3948, 0.4960, 0.6302, 0.7754, 0.8863, 0.9502, 0.9862,
     1.0000},
    {0.0000, 0.0026, 0.0136, 0.0368, 0.0688, 0.1006, 0.1276, 0.1536, 0.1845,
     0.2234, 0.2723, 0.3394, 0.4359, 0.5635, 0.7050, 0.8342, 0.9296, 0.9833,
     1.0000},
    {0.0000, 0.0047, 0.0194, 0.0444, 0.0761, 0.1076, 0.1326, 0.1521, 0.1743,
     0.2063, 0.2490, 0.3036, 0.3820, 0.4985, 0.6479, 0.7999, 0.9170, 0.9812,
     1.0000},
    {0.0000, 0.0008, 0.0024, 0.0250, 0.0631, 0.0940, 0.1133, 0.1345, 0.1630,
     0.1868, 0.2026, 0.2346, 0.3133, 0.4410, 0.5910, 0.7383, 0.8688, 0.9642,
     1.0000},
    {0.0000, 0.0007, 0.0066, 0.0238, 0.0507, 0.0778, 0.0970, 0.1107, 0.1288,
     0.1552, 0.1813, 0.1968, 0.2088, 0.2487, 0.3531, 0.5327, 0.7509, 0.9307,
     1.0000}
  };
}

// Constructor passes arrays to templated base class

G4GamP2NPipAngDst::G4GamP2NPipAngDst(G4int verbose)
  : G4NumIntTwoBodyAngDst<15,19>("G4GamP2NPipAngDist", eBins, angleBins,
				 integralTable, 3.0, verbose) {;}
