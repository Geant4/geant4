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
// Date:    19 May 2013
//
// Description: implementation of numerically integrated angular distribution
//              class for pp and nn elastic reactions
//

#include "G4PP2PPAngDst.hh"

namespace {
  static const G4double eBins[11] =
    { 0.000, 0.130, 0.300, 0.542, 0.817, 1.24, 1.74, 2.25, 2.75, 4.40, 6.15};

  static const G4double angleBins[19] = 
    {-1.000, -0.985, -0.940, -0.866, -0.766, -0.643, -0.500, -0.342, -0.174,
      0.000,  0.174,  0.342,  0.500,  0.643,  0.766,  0.866,  0.940,  0.985,
      1.000};

  static const G4double integralTable[11][19] = {
    {0.0000, 0.0075, 0.0300, 0.0670, 0.1170, 0.1785, 0.2500, 0.3290, 0.4130,
     0.5000, 0.5870, 0.6710, 0.7500, 0.8215, 0.8830, 0.9330, 0.9700, 0.9925, 1.0000},
    {0.0000, 0.0104, 0.0388, 0.0808, 0.1334, 0.1954, 0.2652, 0.3405, 0.4193,
     0.5000, 0.5807, 0.6595, 0.7348, 0.8046, 0.8666, 0.9192, 0.9611, 0.9896, 1.0000},
    {0.0000, 0.0099, 0.0353, 0.0730, 0.1227, 0.1837, 0.2549, 0.3331, 0.4154,
     0.5000, 0.5846, 0.6669, 0.7451, 0.8163, 0.8773, 0.9270, 0.9647, 0.9901, 1.0000},
    {0.0000, 0.0136, 0.0483, 0.0975, 0.1569, 0.2229, 0.2922, 0.3620, 0.4311,
     0.5000, 0.5689, 0.6380, 0.7078, 0.7771, 0.8431, 0.9025, 0.9517, 0.9864, 1.0000},
    {0.0000, 0.0283, 0.0993, 0.1898, 0.2786, 0.3529, 0.4086, 0.4481, 0.4767,
     0.5000, 0.5235, 0.5521, 0.5915, 0.6472, 0.7216, 0.8104, 0.9009, 0.9719, 1.0000},
    {0.0000, 0.0498, 0.1615, 0.2784, 0.3677, 0.4237, 0.4557, 0.4749, 0.4885,
     0.5000, 0.5115, 0.5251, 0.5443, 0.5763, 0.6323, 0.7216, 0.8385, 0.9502, 1.0000},
    {0.0000, 0.0761, 0.2271, 0.3524, 0.4233, 0.4566, 0.4732, 0.4838, 0.4923,
     0.5000, 0.5077, 0.5162, 0.5268, 0.5434, 0.5767, 0.6476, 0.7729, 0.9239, 1.0000},
    {0.0000, 0.0988, 0.2802, 0.4044, 0.4565, 0.4752, 0.4841, 0.4904, 0.4955,
     0.5000, 0.5045, 0.5096, 0.5159, 0.5248, 0.5435, 0.5956, 0.7198, 0.9012, 1.0000},
    {0.0000, 0.1132, 0.3125, 0.4319, 0.4713, 0.4839, 0.4904, 0.4948, 0.4977,
     0.5000, 0.5023, 0.5052, 0.5096, 0.5161, 0.5287, 0.5681, 0.6875, 0.8868, 1.0000},
    {0.0000, 0.1746, 0.3626, 0.4549, 0.4817, 0.4911, 0.4956, 0.4980, 0.4992,
     0.5000, 0.5008, 0.5020, 0.5044, 0.5089, 0.5183, 0.5451, 0.6374, 0.8254, 1.0000},
    {0.0000, 0.2290, 0.4423, 0.4863, 0.4949, 0.4982, 0.4995, 0.4999, 0.5000,
     0.5000, 0.5000, 0.5001, 0.5005, 0.5018, 0.5051, 0.5137, 0.5577, 0.7710, 1.0000}
  }; 
}

// Constructor passes arrays to templated base class

G4PP2PPAngDst::G4PP2PPAngDst(G4int verbose)
 : G4NumIntTwoBodyAngDst<11,19>("G4PP2PPAngDst", eBins, angleBins,
                                integralTable, 7.94, verbose) 
{}

