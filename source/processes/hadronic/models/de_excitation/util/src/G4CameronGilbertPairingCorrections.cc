//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4CameronGilbertPairingCorrections.cc,v 1.2 2005/06/04 13:29:20 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4CameronGilbertPairingCorrections.hh"


// Data comes from:
// A. Gilbert and A.G.W. Cameron, Can. J. Phys., 43, 1446(1965)

// P(Z)
const G4double G4CameronGilbertPairingCorrections::PairingZTable
[G4CameronGilbertPairingCorrections::ZTableSize] = { // 88 values from Z = 11 to Z = 98
    0.00,  2.46,  0.00,  2.09,  0.00,  1.62,  0.00,  1.62,  0.00,  1.83,
    0.00,  1.73,  0.00,  1.35,  0.00,  1.54,  0.00,  1.20,  0.00,  1.06,
    0.00,  1.36,  0.00,  1.43,  0.00,  1.17,  0.00,  1.24,  0.00,  1.20,  
    0.00,  1.28,  0.00,  1.28,  0.00,  1.35,  0.00,  1.36,  0.00,  1.19,
    0.00,  1.14,  0.00,  1.12,  0.00,  1.58,  0.00,  1.17,  0.00,  1.18,
    0.00,  1.22,  0.00,  0.97,  0.00,  0.92,  0.00,  0.62,  0.00,  0.68,
    0.00,  0.64,  0.00,  0.72,  0.00,  0.75,  0.00,  0.71,  0.00,  0.87,
    0.00,  0.83,  0.00,  0.89,  0.00,  0.79,  0.00,  0.89,  0.00,  0.78,
    0.00,  0.69,  0.00,  0.61,  0.00,  0.72,  0.00,  0.77
};
// P(N)
const G4double G4CameronGilbertPairingCorrections::PairingNTable
[G4CameronGilbertPairingCorrections::NTableSize] = { // 140 values from Z = 11 to Z = 150
    0.00,  2.67,  0.00,  1.80,  0.00,  1.67,  0.00,  1.86,  0.00,  2.04,
    0.00,  1.64,  0.00,  1.44,  0.00,  1.54,  0.00,  1.30,  0.00,  1.27,
    0.00,  1.29,  0.00,  1.41,  0.00,  1.50,  0.00,  1.50,  0.00,  1.43,
    0.00,  1.88,  0.00,  1.47,  0.00,  1.57,  0.00,  1.46,  0.00,  0.93,
    0.00,  0.72,  0.00,  1.12,  0.00,  1.29,  0.00,  0.94,  0.00,  1.24,
    0.00,  1.25,  0.00,  1.14,  0.00,  1.32,  0.00,  1.15,  0.00,  1.24,
    0.00,  1.43,  0.00,  1.09,  0.00,  1.20,  0.00,  1.04,  0.00,  0.70,
    0.00,  0.85,  0.00,  0.76,  0.00,  0.92,  0.00,  0.99,  0.00,  1.10,
    0.00,  0.92,  0.00,  0.73,  0.00,  0.70,  0.00,  0.87,  0.00,  0.61,
    0.00,  0.69,  0.00,  0.55,  0.00,  0.40,  0.00,  0.73,  0.00,  0.58,
    0.00,  0.86,  0.00,  1.13,  0.00,  0.84,  0.00,  0.79,  0.00,  0.82,
    0.00,  0.71,  0.00,  0.41,  0.00,  0.38,  0.00,  0.67,  0.00,  0.61,
    0.00,  0.78,  0.00,  0.67,  0.00,  0.67,  0.00,  0.79,  0.00,  0.60,
    0.00,  0.57,  0.00,  0.49,  0.00,  0.43,  0.00,  0.50,  0.00,  0.39
};

G4CameronGilbertPairingCorrections* G4CameronGilbertPairingCorrections::theInstance = 0;

G4CameronGilbertPairingCorrections::G4CameronGilbertPairingCorrections()
{;}


G4CameronGilbertPairingCorrections* G4CameronGilbertPairingCorrections::GetInstance()
{
  if (!theInstance) theInstance = new G4CameronGilbertPairingCorrections();
  return theInstance;
}
