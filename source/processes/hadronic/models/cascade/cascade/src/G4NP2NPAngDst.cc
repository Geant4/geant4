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
// Date:    19 March 2013
//
// Description: implementation of numerically integrated angular distribution
//              class for n p -> n p reaction
//

#include "G4NP2NPAngDst.hh"

namespace {
  static const G4double eBins[11] =
    { 0.000, 0.090, 0.200, 0.300, 0.542, 0.802, 1.240,
      2.250, 4.250, 5.900, 10.00};

  static const G4double angleBins[19] = 
    { -1.000, -0.985, -0.940, -0.866, -0.766, -0.643, -0.500, -0.342, -0.174,
       0.000,  0.174,  0.342,  0.500,  0.643,  0.766,  0.866,  0.940,  0.985,
       1.000 };

  static const G4double integralTable[11][19] = {
    {0.0000, 0.0075, 0.0300, 0.0670, 0.1170, 0.1785, 0.2500, 0.3290, 0.4130,
     0.5000, 0.5870, 0.6710, 0.7500, 0.8215, 0.8830, 0.9330, 0.9700, 0.9925,
     1.0000}, 
    {0.0000, 0.0149, 0.0569, 0.1182, 0.1889, 0.2613, 0.3320, 0.3995, 0.4642,
     0.5264, 0.5858, 0.6428, 0.6998, 0.7596, 0.8229, 0.8872, 0.9450, 0.9855,
     1.0000}, 
    {0.0000, 0.0235, 0.0876, 0.1746, 0.2638, 0.3428, 0.4101, 0.4702, 0.5288,
     0.5873, 0.6421, 0.6897, 0.7313, 0.7731, 0.8219, 0.8795, 0.9384, 0.9833,
     1.0000}, 
    {0.0000, 0.0193, 0.0722, 0.1447, 0.2200, 0.2874, 0.3448, 0.3965, 0.4488,
     0.5062, 0.5685, 0.6331, 0.6983, 0.7637, 0.8290, 0.8923, 0.9478, 0.9863,
     1.0000}, 
    {0.0000, 0.0295, 0.1003, 0.1749, 0.2291, 0.2670, 0.3030, 0.3400, 0.3710,
     0.3959, 0.4280, 0.4811, 0.5543, 0.6352, 0.7165, 0.8028, 0.8942, 0.9700,
     1.0000}, 
    {0.0000, 0.0212, 0.0733, 0.1306, 0.1745, 0.2050, 0.2312, 0.2566, 0.2785,
     0.2958, 0.3152, 0.3473, 0.4003, 0.4782, 0.5832, 0.7129, 0.8505, 0.9588,
     1.0000}, 
    {0.0000, 0.0245, 0.0785, 0.1255, 0.1514, 0.1692, 0.1889, 0.2050, 0.2130,
     0.2204, 0.2342, 0.2512, 0.2745, 0.3292, 0.4451, 0.6201, 0.8090, 0.9493,
     1.0000}, 
    {0.0000, 0.0015, 0.0052, 0.0098, 0.0140, 0.0175, 0.0202, 0.0224, 0.0244,
     0.0266, 0.0295, 0.0337, 0.0416, 0.0568, 0.0920, 0.1773, 0.3946, 0.7789,
     1.0000}, 
    {0.0000, 0.0002, 0.0008, 0.0014, 0.0021, 0.0026, 0.0030, 0.0033, 0.0036,
     0.0040, 0.0044, 0.0050, 0.0061, 0.0082, 0.0135, 0.0372, 0.1874, 0.6556,
     1.0000}, 
    {0.0000, 0.000052, 0.00019, 0.00038, 0.00059, 0.00078, 0.00095, 0.0011,
     0.0013, 0.0014,   0.0016,  0.0019,  0.0022,  0.0028,  0.0045,  0.0147,
     0.1074, 0.5550,   1.0000}, 
    {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     0.0000, 0.0000, 2.5e-22, 5.0e-17, 2.3e-12, 2.3e-08, 0.000043, 0.0111,
     0.3243, 1.0000}
  }; 
}

// Constructor passes arrays to templated base class

G4NP2NPAngDst::G4NP2NPAngDst(G4int verbose)
 : G4NumIntTwoBodyAngDst<11,19>("G4NP2NPAngDst", eBins, angleBins,
                                integralTable, 8.0, verbose) 
{}

