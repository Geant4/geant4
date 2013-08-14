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
// $Id: G4NuclNucl3BodyMomDst.cc 67796 2013-03-08 06:18:39Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    7 March 2013
//
// Description: class containing parametrized momentum distributions
//              in the CM for nucleon/nucleon 3-body final states

#include "G4NuclNucl3BodyMomDst.hh"

namespace {
  // Powers of Ekin^0..3, blocks of S^0..3 for PQ,PR: outgoing N; outgoing h,K,Y
  static const G4double pqprC[2][4][4] = {
    { { 0.5028,  0.9348, -0.0967,  -0.025},
      { 3.1442,  -10.59,  4.7335, -0.6248},
      {-7.8172,  29.227, -14.298,  2.0282},
      { 8.1667,  -34.55,  17.685, -2.5895} },
    { { 1.1965,   0.287, -0.2449,  0.0373},
      {-0.8289, -4.9065,  2.9191,  -0.422},
      { 1.0426,  16.264, -9.5776,  1.3883},
      { -1.909, -19.904,  11.938, -1.7476} }
  };

  // Powers of Ekin^0..2 for PS: outgoing N; outgoing h,K,Y
  static const G4double psC[2][3] = {
    { 0.1451,  0.4652,  -0.033 }, { 0.1538,  0.2744, -0.0146 }
  };
}

// Constructor passes arrays to templated base class

G4NuclNucl3BodyMomDst::G4NuclNucl3BodyMomDst(G4int verbose)
  : G4InuclParamMomDst("G4NuclNucl3BodyMomDist", pqprC, psC, verbose) {;}
