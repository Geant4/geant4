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
// $Id: G4HadNucl3BodyMomDst.cc 67796 2013-03-08 06:18:39Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    7 March 2013
//
// Description: class containing parametrized momentum distributions
//              in the CM for hadron/nucleon 3-body final states

#include "G4HadNucl3BodyMomDst.hh"

namespace {
  // Powers of Ekin^0..3, blocks of S^0..3 for PQ,PR: outgoing N; outgoing h,K,Y
  static const G4double pqprC[2][4][4] = {
    { { 0.6305,  2.1801, -1.2886,  0.2091},
      {-3.7333,  1.5163,  -2.457,  0.5228},
      { 13.464,  -16.38,  15.129, -2.8687},
      {-18.594,  27.944, -23.295,  4.2688} },
    { { 0.9336,  1.7811, -1.5264,  0.2713},
      {-1.8181, -8.2927,  6.8433, -1.1944},
      { 5.5157,  20.607, -16.067,  2.7495},
      {-8.5216, -20.827,  16.845, -2.9045} }
  };

  // Powers of Ekin^0..2 for PS: outgoing N; outgoing h,K,Y
  static const G4double psC[2][3] = {
    { 0.0929,  0.5389, -0.0545 }, { 0.1303,  0.4071, -0.0288 }
  };
}

// Constructor passes arrays to templated base class

G4HadNucl3BodyMomDst::G4HadNucl3BodyMomDst(G4int verbose)
  : G4InuclParamMomDst("G4HadNucl3BodyMomDist", pqprC, psC, verbose) {;}
