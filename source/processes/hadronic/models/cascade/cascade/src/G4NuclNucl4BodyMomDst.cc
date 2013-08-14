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
// $Id: G4NuclNucl4BodyMomDst.cc 67796 2013-03-08 06:18:39Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    7 March 2013
//
// Description: class containing parametrized momentum distributions
//              in the CM for nucleon/nucleon >= 4-body final states

#include "G4NuclNucl4BodyMomDst.hh"

namespace {
  // Powers of Ekin^0..3, blocks of S^0..3 for PQ,PR: outgoing N; outgoing h,K,Y
  static const G4double pqprC[2][4][4] = {
    { { 1.6208, -0.2009,  0.0126, -0.0002},
      {-4.3139,  1.3641, -0.0835,  0.0014},
      { 12.291,  -3.403,   0.186, -0.0024},
      {-15.288,  3.8559, -0.2004,  0.0022} },
    { {1.2419,  -0.244,  0.0157, -0.0003},
      {-4.3633,  1.3158, -0.0826,  0.0014},
      {13.743, -3.5691,  0.2143, -0.0034},
      {-18.592,  4.3867, -0.2585,  0.0039} }
  };

  // Powers of Ekin^0..2 for PS: outgoing N; outgoing h,K,Y
  static const G4double psC[2][3] = {
    { 0.6296, 0.1787, -0.0026 }, { 0.8381, 0.0086, 0.0033 }
  };
}

// Constructor passes arrays to templated base class

G4NuclNucl4BodyMomDst::G4NuclNucl4BodyMomDst(G4int verbose)
  : G4InuclParamMomDst("G4NuclNucl4BodyMomDist", pqprC, psC, verbose) {;}
