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
// $Id: $
// Author:  Dennis Wright (SLAC)
// Date:    6 January 2014
//
// Description: implementation of numerically integrated angular distribution
//              class for pi+ p and pi- n elastic scattering
//

#include "G4PipP2PipPAngDst.hh"

namespace {
  static const G4double eBins[11] =
    { 0.000, 0.031, 0.097, 0.240, 0.550, 0.873, 1.31, 1.53, 2.60, 3.86, 10.66};

  static const G4double angleBins[19] = 
    {-1.000, -0.985, -0.940, -0.866, -0.766, -0.643, -0.500, -0.342, -0.174,
      0.000,  0.174,  0.342,  0.500,  0.643,  0.766,  0.866,  0.940,  0.985,
      1.000};

  static const G4double integralTable[11][19] = {
    {0.0000, 0.0110, 0.0438, 0.0964, 0.1655, 0.2470, 0.3370, 0.4308, 0.5239,
     0.6133, 0.6959, 0.7695, 0.8329, 0.8856, 0.9278, 0.9600, 0.9825, 0.9957, 1.0000},
    {0.0000, 0.0197, 0.0771, 0.1667, 0.2785, 0.4016, 0.5255, 0.6411, 0.7413,
     0.8228, 0.8849, 0.9292, 0.9588, 0.9773, 0.9883, 0.9945, 0.9978, 0.9995, 1.0000},
    {0.0000, 0.0180, 0.0701, 0.1498, 0.2461, 0.3484, 0.4473, 0.5355, 0.6099,
     0.6711, 0.7221, 0.7668, 0.8096, 0.8526, 0.8955, 0.9360, 0.9698, 0.9922, 1.0000},
    {0.0000, 0.0113, 0.0434, 0.0910, 0.1461, 0.2011, 0.2510, 0.2935, 0.3309,
     0.3689, 0.4145, 0.4738, 0.5502, 0.6418, 0.7413, 0.8383, 0.9243, 0.9811, 1.0000},
    {0.0000, 0.0007, 0.0034, 0.0086, 0.0154, 0.0219, 0.0264, 0.0291, 0.0356,
     0.0566, 0.1044, 0.1880, 0.3089, 0.4594, 0.6218, 0.7746, 0.8971, 0.9741, 1.0000},
    {0.0000, 0.0053, 0.0204, 0.0428, 0.0694, 0.0999, 0.1335, 0.1653, 0.1911,
     0.2147, 0.2483, 0.3057, 0.3956, 0.5161, 0.6535, 0.7884, 0.9009, 0.9746, 1.0000},
    {0.0000, 0.0112, 0.0372, 0.0640, 0.0878, 0.1198, 0.1709, 0.2346, 0.2914,
     0.3302, 0.3530, 0.3677, 0.3800, 0.4012, 0.4583, 0.5797, 0.7591, 0.9294, 1.0000},
    {0.0000, 0.0043, 0.0141, 0.0240, 0.0349, 0.0555, 0.0930, 0.1414, 0.1879,
     0.2252, 0.2548, 0.2770, 0.2904, 0.3077, 0.3637, 0.4978, 0.7079, 0.9136, 1.0000},
    {0.0000, 0.0035, 0.0095, 0.0143, 0.0206, 0.0281, 0.0351, 0.0408, 0.0441,
     0.0466, 0.0533, 0.0655, 0.0806, 0.0957, 0.1238, 0.2268, 0.4855, 0.8311, 1.0000},
    {0.0000, 0.0009, 0.0025, 0.0041, 0.0053, 0.0062, 0.0068, 0.0073, 0.0076,
     0.0078, 0.0079, 0.0085, 0.0109, 0.0173, 0.0327, 0.0869, 0.3238, 0.7572, 1.0000},
    {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0002, 0.0269, 0.4555, 1.0000} 
  }; 
}

// Constructor passes arrays to templated base class

G4PipP2PipPAngDst::G4PipP2PipPAngDst(G4int verbose)
 : G4NumIntTwoBodyAngDst<11,19>("G4PipP2PipPAngDst", eBins, angleBins,
                                integralTable, 7.43, verbose) 
{}

