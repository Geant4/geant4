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
// $Id: G4KaonSampler.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 20150323  Split KaonHypSampler to allow different energy binning
// 20150713  Change binning arguments to match new K+p etc. tables

#include "G4KaonSampler.hh"

namespace {
  static const G4double bins[30] = 
    { 0.0,  0.01, 0.02, 0.04, 0.06, 0.08, 0.11, 0.14, 0.21, 0.28,
      0.36, 0.45, 0.53, 0.67, 0.80, 0.90, 0.99, 1.23, 1.47, 1.9,
      2.4,  3.2,  4.2,  5.6,  7.5,  10.0, 13.0, 18.0, 24.0, 32.0 };
}

// Define constructor for initialization
G4KaonSampler::G4KaonSampler() : G4CascadeSampler<30,8>(bins) {}
