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
// $Id: G4NuclNuclAngDst.cc 67633 2013-02-27 20:38:04Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    20 February 2013
//
// Description: class containing parametrized angular distributions
//              in the CM for nucleon/hyperon - nucleon two-body scattering

#include "G4NuclNuclAngDst.hh"

namespace {
  static const G4double nnke[9] =
    { 0.0, 0.44, 0.59, 0.80, 1.00, 2.24, 4.40, 6.15, 10.00};
  static const G4double nnFrac[9] =
    {1.0,   1.0, 0.4898, 0.7243, 0.7990, 0.8892, 0.8493, 0.9583,  1.0};
  static const G4double nnA[9] =
    { 0.0,   0.34, 2.51,   4.59,   4.2,    5.61,   6.38,   7.93,   8.7};
  static const G4double nnC[9] =
    { 0.0,   0.0,  1.21,   1.54,   1.88,   1.24,   1.91,   4.04,   8.7};
  static const G4double nnCos[9] =
    {-1.0,  -1.0, 0.4226, 0.4226, 0.4384, 0.7193, 0.8788, 0.9164,  0.95};
}

// Constructor passes arrays to templated base class

G4NuclNuclAngDst::G4NuclNuclAngDst(G4int verbose)
  : G4ParamExpTwoBodyAngDst<9>("G4NuclNuclAngDist", nnke, nnFrac,
			     nnA, nnC, nnCos, verbose) {;}
