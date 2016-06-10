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
//              class for pi- p and pi+ n elastic scattering
//

#include "G4PimP2PimPAngDst.hh"

namespace {
  static const G4double eBins[11] =
    { 0.000, 0.031, 0.097, 0.240, 0.550, 0.873, 1.31, 1.53, 2.60, 3.86, 10.66};

  static const G4double angleBins[19] = 
    {-1.000, -0.985, -0.940, -0.866, -0.766, -0.643, -0.500, -0.342, -0.174,
      0.000,  0.174,  0.342,  0.500,  0.643,  0.766,  0.866,  0.940,  0.985,
      1.000};

  static const G4double integralTable[11][19] = {
    {0.0000, 0.0065, 0.0262, 0.0590, 0.1037, 0.1595, 0.2258, 0.3007, 0.3820,
     0.4682, 0.5562, 0.6428, 0.7260, 0.8028, 0.8699, 0.9251, 0.9664, 0.9916, 1.0000}, 
    {0.0000, 0.0023, 0.0098, 0.0237, 0.0457, 0.0776, 0.1215, 0.1784, 0.2481,
     0.3301, 0.4224, 0.5209, 0.6223, 0.7213, 0.8122, 0.8900, 0.9499, 0.9874, 1.0000}, 
    {0.0000, 0.0045, 0.0179, 0.0395, 0.0678, 0.1017, 0.1408, 0.1849, 0.2350,
     0.2933, 0.3621, 0.4426, 0.5356, 0.6381, 0.7437, 0.8434, 0.9263, 0.9810, 1.0000}, 
    {0.0000, 0.0143, 0.0550, 0.1149, 0.1834, 0.2512, 0.3119, 0.3630, 0.4075,
     0.4516, 0.5019, 0.5630, 0.6365, 0.7192, 0.8037, 0.8819, 0.9453, 0.9860, 1.0000}, 
    {0.0000, 0.0097, 0.0368, 0.0745, 0.1137, 0.1464, 0.1682, 0.1799, 0.1879,
     0.2032, 0.2366, 0.2953, 0.3812, 0.4907, 0.6158, 0.7460, 0.8668, 0.9654, 1.0000},
    {0.0000, 0.0015, 0.0123, 0.0435, 0.0982, 0.1650, 0.2257, 0.2681, 0.2928,
     0.3080, 0.3195, 0.3272, 0.3353, 0.3627, 0.4392, 0.5820, 0.7691, 0.9343, 1.0000}, 
    {0.0000, 0.0029, 0.0108, 0.0221, 0.0369, 0.0571, 0.0840, 0.1158, 0.1492,
     0.1791, 0.1997, 0.2099, 0.2191, 0.2506, 0.3352, 0.4943, 0.7130, 0.9153, 1.0000}, 
    {0.0000, 0.0018, 0.0057, 0.0094, 0.0133, 0.0222, 0.0390, 0.0620, 0.0874,
     0.1142, 0.1395, 0.1578, 0.1711, 0.1984, 0.2731, 0.4292, 0.6653, 0.8989, 1.0000}, 
    {0.0000, 0.0005, 0.0016, 0.0030, 0.0043, 0.0051, 0.0058, 0.0069, 0.0085,
     0.0123, 0.0206, 0.0329, 0.0471, 0.0607, 0.0866, 0.1828, 0.4420, 0.8093, 1.0000}, 
    {0.0000, 0.0009, 0.0025, 0.0041, 0.0053, 0.0062, 0.0068, 0.0073, 0.0076,
     0.0078, 0.0079, 0.0085, 0.0109, 0.0173, 0.0327, 0.0869, 0.3238, 0.7572, 1.0000}, 
    {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0002, 0.0269, 0.4555, 1.0000}
  }; 
}

// Constructor passes arrays to templated base class

G4PimP2PimPAngDst::G4PimP2PimPAngDst(G4int verbose)
 : G4NumIntTwoBodyAngDst<11,19>("G4PimP2PimPAngDst", eBins, angleBins,
                                integralTable, 7.43, verbose) 
{}

