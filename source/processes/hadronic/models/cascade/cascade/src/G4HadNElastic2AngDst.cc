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
// $Id: G4HadNElastic2AngDst.cc 67633 2013-02-27 20:38:04Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    20 February 2013
//
// Description: class containing parametrized angular distributions
//              in the CM for elastic scattering of pi-p, pi+n, k-p,
//		k0bn, k+n, or k0p

#include "G4HadNElastic2AngDst.hh"

namespace {
  static const G4double hn2ke[10] = 
    {0.0, 0.062, 0.12, 0.217, 0.533, 0.873, 1.34, 2.86, 5.86, 10.0};

  static const G4double hn2A[10] =
    {0.0, 27.08, 19.32,   9.92,   7.74,   9.86,   5.51,   7.25,   7.23,   7.3};
  static const G4double hn2C[10] =
    {0.0,  0.0, -19.49, -15.78,  -9.78,  -2.74,  -1.16,   2.31,   2.96,   7.3};
  static const G4double hn2Cos[10] =
    {-1.0, -1.0,  -0.235, -0.259, -0.276,  0.336,  0.250,  0.732,  0.875,  0.9};
  static const G4double hn2Frac[10] =
    {1.0,  1.0,   0.6918, 0.6419, 0.7821, 0.6542, 0.8382, 0.9722, 0.9784, 1.0};
}

// Constructor passes arrays to templated base class

G4HadNElastic2AngDst::G4HadNElastic2AngDst(G4int verbose)
  : G4ParamExpTwoBodyAngDst<10>("G4HadNElastic2AngDist", hn2ke, hn2Frac,
			     hn2A, hn2C, hn2Cos, verbose) {;}
