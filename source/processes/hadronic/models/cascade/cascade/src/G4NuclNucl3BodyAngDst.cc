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
//              in the CM for nucleon/nucleon 3-body final states

#include "G4NuclNucl3BodyAngDst.hh"

namespace {
  // Powers of Ekin^0..3, blocks of S^0..3 for AB
  static const G4double abC[2][4][4] = {
    // -------- Initial state nucleon-nucleon, outgoing N --------
    { { 0.0856, 0.0543,-0.0511, 0.0075 }, {  5.039,-9.2324, 4.6003,-0.6253 },
      {-13.782, 36.397,-20.534, 2.9159 }, { 14.661,-42.962, 27.731,-4.1101 } 
    },
    // -------- Initial state nucleon-nucleon, outgoing h,K,Y --------
    { { 0.0716, 0.0926,-0.0515, 0.0058 }, {  3.096,-3.2186, 0.8989,-0.0017 },
      {-11.125, 20.273,-7.5084, 0.7022 }, {  18.13,-33.245, 13.188,-1.4854 } 
    }
  };
}

// Constructor passes arrays to templated base class

G4NuclNucl3BodyAngDst::G4NuclNucl3BodyAngDst(G4int verbose)
  : G4InuclParamAngDst("G4NuclNucl3BodyAngDist", abC, verbose) {;}
