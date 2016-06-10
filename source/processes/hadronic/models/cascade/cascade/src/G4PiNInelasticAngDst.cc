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
// $Id: G4PiNInelasticAngDst.cc 67633 2013-02-27 20:38:04Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    20 February 2013
//
// Description: class containing parametrized angular distributions
//              in the CM for pion - nucleon inelastic two-body scattering

#include "G4PiNInelasticAngDst.hh"

namespace {
  static const G4double qxke[10] =
    {0.0, 0.062, 0.12, 0.217, 0.533, 0.873, 1.34, 2.86, 5.86, 10.0};

  static const G4double qxFrac[10] =
    {0.0,   0.0,    0.1156, 0.5832, 0.8125, 0.3357, 0.3269, 0.7765, 0.8633, 1.0};
  static const G4double qxA[10] =
    {0.0,   0.0,    2.48,   7.93,  10.0,    9.78,   5.08,   8.13,   8.13,   8.13};
  static const G4double qxC[10] =
    {0.0, -39.58, -12.55,  -4.38,   1.81,  -1.99,  -0.33,   1.2,    1.43,   8.13};
  static const G4double qxCos[10] =
    {1.0,   1.0,    0.604, -0.033,  0.25,   0.55,   0.65,   0.80,   0.916,  0.916};
}

// Constructor passes arrays to templated base class

G4PiNInelasticAngDst::G4PiNInelasticAngDst(G4int verbose)
  : G4ParamExpTwoBodyAngDst<10>("G4PiNInelasticAngDist", qxke, qxFrac,
			     qxA, qxC, qxCos, verbose) {;}
