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
// $Id$
// Author:  Michael Kelsey (SLAC)
// Date:    22 April 2013
//
// Description: class containing parametrized angular distributions
//              in the CM for hadron/nucleon 3-body final states

#include "G4HadNucl3BodyAngDst.hh"

namespace {
  // Powers of Ekin^0..3, blocks of S^0..3 for AB
  static const G4double abC[2][4][4] = {
    // -------- Initial state (h,K,Y,g)-nucleon, outgoing N --------
    { { 0.1729, -0.145, 0.0454,-0.0048 }, {  7.108,-13.032, 8.3515,-1.4095 },
      {-17.961, 41.781, -30.26, 5.3505 }, { 16.403,-40.799, 32.882,-6.0946 } 
    },
    // -------- Initial state (h,K,Y,g)-nucleon, outgoing h,K,Y --------
    { { 0.0376, 0.2383,-0.1541,  0.025 }, { 1.4331, 1.8253,-1.5201, 0.3059 },
      { -3.135, 1.7648,-1.5692, 0.3252 }, { 6.4864,-16.735, 17.185,-3.5277 } 
    }
  };
}

// Constructor passes arrays to templated base class

G4HadNucl3BodyAngDst::G4HadNucl3BodyAngDst(G4int verbose)
  : G4InuclParamAngDst("G4HadNucl3BodyAngDist", abC, verbose) {;}
