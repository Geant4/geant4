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
//              class for pi- p charge exchange scattering
//

#include "G4PimP2Pi0NAngDst.hh"

namespace {
  static const G4double eBins[11] =
    { 0.000, 0.031, 0.097, 0.240, 0.550, 0.873, 1.31, 1.53, 2.60, 3.86, 10.66};

  static const G4double angleBins[19] = 
    {-1.000, -0.985, -0.940, -0.866, -0.766, -0.643, -0.500, -0.342, -0.174,
      0.000,  0.174,  0.342,  0.500,  0.643,  0.766,  0.866,  0.940,  0.985,
      1.000};

  static const G4double integralTable[11][19] = {
    {0.0000, 0.0097, 0.0384, 0.0851, 0.1467, 0.2204, 0.3033, 0.3914, 0.4810,
     0.5695, 0.6539, 0.7316, 0.8011, 0.8612, 0.9108, 0.9498, 0.9778, 0.9945, 1.0000}, 
    {0.0000, 0.0175, 0.0687, 0.1489, 0.2501, 0.3634, 0.4804, 0.5925, 0.6934,
     0.7798, 0.8493, 0.9021, 0.9402, 0.9659, 0.9822, 0.9917, 0.9969, 0.9993, 1.0000}, 
    {0.0000, 0.0217, 0.0846, 0.1805, 0.2960, 0.4178, 0.5345, 0.6369, 0.7202,
     0.7846, 0.8329, 0.8692, 0.8983, 0.9237, 0.9469, 0.9678, 0.9849, 0.9961, 1.0000}, 
    {0.0000, 0.0098, 0.0379, 0.0798, 0.1287, 0.1781, 0.2236, 0.2632, 0.2988,
     0.3357, 0.3806, 0.4398, 0.5178, 0.6134, 0.7196, 0.8257, 0.9170, 0.9785, 1.0000}, 
    {0.0000, 0.0054, 0.0196, 0.0377, 0.0539, 0.0656, 0.0741, 0.0837, 0.0995,
     0.1254, 0.1644, 0.2207, 0.3010, 0.4113, 0.5504, 0.7067, 0.8550, 0.9614, 1.0000}, 
    {0.0000, 0.0058, 0.0353, 0.1089, 0.2242, 0.3480, 0.4383, 0.4785, 0.4902,
     0.5098, 0.5536, 0.6067, 0.6457, 0.6688, 0.7015, 0.7706, 0.8714, 0.9631, 1.0000}, 
    {0.0000, 0.0204, 0.0697, 0.1265, 0.1852, 0.2546, 0.3326, 0.4007, 0.4506,
     0.4997, 0.5632, 0.6281, 0.6754, 0.7154, 0.7784, 0.8691, 0.9513, 0.9918, 1.0000}, 
    {0.0000, 0.0093, 0.0328, 0.0640, 0.1041, 0.1580, 0.2202, 0.2780, 0.3369,
     0.4192, 0.5215, 0.6020, 0.6368, 0.6540, 0.6963, 0.7805, 0.8874, 0.9707, 1.0000}, 
    {0.0000, 0.0086, 0.0305, 0.0549, 0.0735, 0.0905, 0.1111, 0.1319, 0.1481,
     0.1669, 0.2076, 0.2882, 0.4029, 0.5084, 0.5852, 0.6656, 0.7921, 0.9352, 1.0000}, 
    {0.0000, 0.0010, 0.0027, 0.0043, 0.0056, 0.0065, 0.0072, 0.0077, 0.0081,
     0.0090, 0.0108, 0.0138, 0.0224, 0.0534, 0.1173, 0.1715, 0.3494, 0.7745, 1.0000}, 
    {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0002, 0.0269, 0.4555, 1.0000}
  }; 
}

// Constructor passes arrays to templated base class

G4PimP2Pi0NAngDst::G4PimP2Pi0NAngDst(G4int verbose)
 : G4NumIntTwoBodyAngDst<11,19>("G4PimP2Pi0NAngDst", eBins, angleBins,
                                integralTable, 7.43, verbose) 
{}

