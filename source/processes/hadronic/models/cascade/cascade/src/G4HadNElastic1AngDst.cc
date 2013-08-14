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
// $Id: G4HadNElastic1AngDst.cc 67633 2013-02-27 20:38:04Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    20 February 2013
//
// Description: class containing parametrized angular distributions
//              in the CM for elastic scattering of pi+p, pi0p, gammap,
//		k+p, k0bp, pi-n, pi0n, gamman, k-n, or k0n

#include "G4HadNElastic1AngDst.hh"

namespace {
  static const G4double hn1ke[10] =
    {0.0, 0.062, 0.12, 0.217, 0.533, 0.873, 1.34, 2.86, 5.86, 10.0};

  static const G4double hn1A[10] =
    {0.0,  0.0,   27.58,  19.83,   6.46,   4.59,   6.47,   6.68,   6.43,   6.7};
  static const G4double hn1C[10] =
    {0.0, -26.4, -30.55, -19.42,  -5.05,  -5.24,  -1.00,   2.14,   2.9,    6.7};
  static const G4double hn1Cos[10] =
    {1.0,  1.0,    0.174, -0.174, -0.7,   -0.295,  0.5,    0.732,  0.837,  0.89};
  static const G4double hn1Frac[10] =
    {0.0,  0.0,    0.2980, 0.7196, 0.9812, 0.8363, 0.5602, 0.9601, 0.9901, 1.0};
}

// Constructor passes arrays to templated base class

G4HadNElastic1AngDst::G4HadNElastic1AngDst(G4int verbose)
  : G4ParamExpTwoBodyAngDst<10>("G4HadNElastic1AngDist", hn1ke, hn1Frac,
			     hn1A, hn1C, hn1Cos, verbose) {;}
