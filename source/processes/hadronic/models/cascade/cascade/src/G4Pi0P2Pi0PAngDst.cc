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
//              class for pi0 p and pi0 n elastic scattering
//

#include "G4Pi0P2Pi0PAngDst.hh"

namespace {
  static const G4double eBins[11] =
    { 0.000, 0.031, 0.097, 0.240, 0.550, 0.873, 1.31, 1.53, 2.60, 3.86, 10.66};

  static const G4double angleBins[19] = 
    {-1.000, -0.985, -0.940, -0.866, -0.766, -0.643, -0.500, -0.342, -0.174,
      0.000,  0.174,  0.342,  0.500,  0.643,  0.766,  0.866,  0.940,  0.985,
      1.000};

  static const G4double integralTable[11][19] = {
    {0.0000, 0.0088, 0.0350, 0.0777, 0.1346, 0.2032, 0.2814, 0.3658, 0.4530,
     0.5408, 0.6260, 0.7062, 0.7794, 0.8442, 0.8988, 0.9426, 0.9744, 0.9936, 1.0000},
    {0.0000, 0.0110, 0.0434, 0.0952, 0.1621, 0.2396, 0.3235, 0.4098, 0.4947,
     0.5764, 0.6536, 0.7250, 0.7906, 0.8493, 0.9002, 0.9422, 0.9738, 0.9934, 1.0000},
    {0.0000, 0.0112, 0.0440, 0.0946, 0.1570, 0.2250, 0.2940, 0.3602, 0.4224,
     0.4822, 0.5421, 0.6047, 0.6726, 0.7454, 0.8196, 0.8897, 0.9480, 0.9866, 1.0000},
    {0.0000, 0.0128, 0.0492, 0.1030, 0.1648, 0.2262, 0.2814, 0.3282, 0.3692,
     0.4102, 0.4582, 0.5184, 0.5934, 0.6805, 0.7725, 0.8601, 0.9348, 0.9836, 1.0000},
    {0.0000, 0.0052, 0.0201, 0.0416, 0.0646, 0.0842, 0.0973, 0.1045, 0.1118,
     0.1299, 0.1705, 0.2416, 0.3450, 0.4750, 0.6188, 0.7603, 0.8820, 0.9698, 1.0000},
    {0.0000, 0.0034, 0.0164, 0.0432, 0.0838, 0.1324, 0.1796, 0.2167, 0.2420,
     0.2614, 0.2839, 0.3164, 0.3654, 0.4394, 0.5464, 0.6852, 0.8350, 0.9544, 1.0000},
    {0.0000, 0.0070, 0.0240, 0.0430, 0.0624, 0.0884, 0.1274, 0.1752, 0.2203,
     0.2546, 0.2764, 0.2888, 0.2996, 0.3259, 0.3968, 0.5370, 0.7360, 0.9224, 1.0000},
    {0.0000, 0.0030, 0.0099, 0.0167, 0.0241, 0.0388, 0.0660, 0.1017, 0.1376,
     0.1697, 0.1972, 0.2174, 0.2308, 0.2530, 0.3184, 0.4635, 0.6866, 0.9062, 1.0000},
    {0.0000, 0.0020, 0.0056, 0.0086, 0.0124, 0.0166, 0.0204, 0.0238, 0.0263,
     0.0294, 0.0370, 0.0492, 0.0638, 0.0782, 0.1052, 0.2048, 0.4638, 0.8202, 1.0000},
    {0.0000, 0.0009, 0.0025, 0.0041, 0.0053, 0.0062, 0.0068, 0.0073, 0.0076,
     0.0078, 0.0079, 0.0085, 0.0109, 0.0173, 0.0327, 0.0869, 0.3238, 0.7572, 1.0000},
    {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0002, 0.0269, 0.4555, 1.0000}
  }; 
}

// Constructor passes arrays to templated base class

G4Pi0P2Pi0PAngDst::G4Pi0P2Pi0PAngDst(G4int verbose)
 : G4NumIntTwoBodyAngDst<11,19>("G4Pi0P2Pi0PAngDst", eBins, angleBins,
                                integralTable, 7.43, verbose) 
{}

