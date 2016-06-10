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
// $Id: G4HadNucl4BodyMomDst.cc 67874 2013-03-12 05:37:09Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    7 March 2013
//
// Description: class containing parametrized momentum distributions
//              in the CM for hadron/nucleon >= 4-body final states

#include "G4HadNucl4BodyMomDst.hh"

namespace {
  // Powers of Ekin^0..3, blocks of S^0..3 for PQ,PR: outgoing N; outgoing h,K,Y
  static const G4double pqprC[2][4][4] = {
    { { 1.9439, -0.3464,  0.0271, -0.0007},
      {-4.6268,  1.1093, -0.1164,  0.0051},
      { 9.7879, -1.9313,  0.2697,  -0.015},
      {-9.6074,  1.7064, -0.3185,  0.0196} },
    { { 1.8693, -0.4996,  0.0462, -0.0013},
      {-5.5678,  1.7874, -0.1854,  0.0058},
      { 14.795,  -4.133,  0.4531, -0.0146},
      {-16.903,  3.8393, -0.4627,  0.0156} }
  };

  // Powers of Ekin^0..2 for PS: outgoing N; outgoing h,K,Y
  static const G4double psC[2][3] = {
    { 0.1491, 0.385, -0.0128 }, { 0.1802, 0.3302, -0.0094 }
  };
}

// Constructor passes arrays to templated base class

G4HadNucl4BodyMomDst::G4HadNucl4BodyMomDst(G4int verbose)
  : G4InuclParamMomDst("G4HadNucl4BodyMomDist", pqprC, psC, verbose) {;}
