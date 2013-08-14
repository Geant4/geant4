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
// $Id: G4GamP2PPi0AngDst.cc 67433 2013-02-20 21:27:31Z mkelsey $
// Author:  Dennis Wright (SLAC)
// Date:    28 January 2013
//
// Description: implementation of numerically integrated angular distribution
//              class for gamma p -> p pi0 reaction
//
// 20130219	Inherit from templated base, move arrays to namespace statics

#include "G4GamP2PPi0AngDst.hh"

namespace {
  static const G4double eBins[15] =
    { 0.145, 0.169, 0.240, 0.300, 0.360, 0.420, 0.480,
      0.700, 1.050, 1.238, 1.400, 1.575, 1.825, 2.450,
      2.900 };

  static const G4double angleBins[19] = 
    { -1.000, -0.985, -0.940, -0.866, -0.766, -0.643, -0.500, -0.342, -0.174,
       0.000,  0.174,  0.342,  0.500,  0.643,  0.766,  0.866,  0.940,  0.985,
       1.000 };

  static const G4double integralTable[15][19] = {
    {0.0000, 0.0075, 0.0300, 0.0670, 0.1170, 0.1785, 0.2500, 0.3290, 0.4130,
     0.5000, 0.5870, 0.6710, 0.7500, 0.8215, 0.8830, 0.9330, 0.9700, 0.9925,
     1.0000},
    {0.0000, 0.0123, 0.0482, 0.1041, 0.1739, 0.2514, 0.3315, 0.4098, 0.4840,
     0.5543, 0.6214, 0.6862, 0.7500, 0.8124, 0.8711, 0.9231, 0.9644, 0.9909,
     1.0000},
    {0.0000, 0.0059, 0.0242, 0.0564, 0.1039, 0.1676, 0.2475, 0.3409, 0.4434,
     0.5499, 0.6534, 0.7473, 0.8274, 0.8910, 0.9374, 0.9686, 0.9875, 0.9971,
     1.0000},
    {0.0000, 0.0039, 0.0164, 0.0399, 0.0769, 0.1300, 0.2003, 0.2865, 0.3851,
     0.4913, 0.5980, 0.6981, 0.7866, 0.8598, 0.9158, 0.9558, 0.9816, 0.9956,
     1.0000},
    {0.0000, 0.0025, 0.0113, 0.0294, 0.0605, 0.1082, 0.1749, 0.2598, 0.3594,
     0.4686, 0.5798, 0.6848, 0.7780, 0.8550, 0.9136, 0.9549, 0.9814, 0.9956,
     1.0000},
    {0.0000, 0.0015, 0.0074, 0.0212, 0.0476, 0.0910, 0.1544, 0.2377, 0.3374,
     0.4485, 0.5627, 0.6714, 0.7684, 0.8486, 0.9098, 0.9530, 0.9806, 0.9954,
     1.0000},
    {0.0000, 0.0009, 0.0051, 0.0160, 0.0383, 0.0765, 0.1342, 0.2116, 0.3062,
     0.4135, 0.5260, 0.6356, 0.7360, 0.8220, 0.8902, 0.9408, 0.9748, 0.9939,
     1.0000},
    {0.0000, 0.0022, 0.0105, 0.0290, 0.0630, 0.1170, 0.1937, 0.2914, 0.4047,
     0.5261, 0.6454, 0.7529, 0.8421, 0.9092, 0.9542, 0.9807, 0.9938, 0.9988,
     1.0000},
    {0.0000, 0.0080, 0.0333, 0.0785, 0.1424, 0.2194, 0.3003, 0.3759, 0.4433,
     0.5071, 0.5757, 0.6547, 0.7432, 0.8313, 0.9058, 0.9578, 0.9863, 0.9975,
     1.0000},
    {0.0000, 0.0064, 0.0288, 0.0717, 0.1333, 0.2036, 0.2718, 0.3330, 0.3919,
     0.4566, 0.5307, 0.6114, 0.6937, 0.7726, 0.8448, 0.9078, 0.9575, 0.9892,
     1.0000},
    {0.0000, 0.0129, 0.0517, 0.1132, 0.1849, 0.2517, 0.3063, 0.3558, 0.4142,
     0.4878, 0.5668, 0.6343, 0.6844, 0.7265, 0.7772, 0.8444, 0.9189, 0.9778,
     1.0000},
    {0.0000, 0.0027, 0.0197, 0.0638, 0.1277, 0.1895, 0.2402, 0.2972, 0.3824,
     0.4892, 0.5796, 0.6246, 0.6388, 0.6663, 0.7363, 0.8357, 0.9271, 0.9827,
     1.0000},
    {0.0000, 0.0025, 0.0142, 0.0401, 0.0761, 0.1140, 0.1572, 0.2217, 0.3174,
     0.4282, 0.5179, 0.5649, 0.5863, 0.6227, 0.6997, 0.8068, 0.9095, 0.9773,
     1.0000},
    {0.0000, 0.0010, 0.0036, 0.0062, 0.0083, 0.0111, 0.0174, 0.0302, 0.0518,
     0.0844, 0.1307, 0.1957, 0.2864, 0.4073, 0.5541, 0.7132, 0.8598, 0.9630,
     1.0000},
    {0.0000, 0.0013, 0.0044, 0.0079, 0.0111, 0.0149, 0.0219, 0.0351, 0.0569,
     0.0895, 0.1359, 0.2010, 0.2917, 0.4121, 0.5582, 0.7160, 0.8612, 0.9634,
     1.0000}
  };
}

// Constructor passes arrays to templated base class

G4GamP2PPi0AngDst::G4GamP2PPi0AngDst(G4int verbose)
  : G4NumIntTwoBodyAngDst<15,19>("G4GamP2PPi0AngDist", eBins, angleBins,
				 integralTable, 1.5, verbose) {;}
