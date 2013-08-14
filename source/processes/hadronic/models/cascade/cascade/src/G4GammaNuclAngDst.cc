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
// $Id: G4GammaNuclAngDst.cc 67633 2013-02-27 20:38:04Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    20 February 2013
//
// Description: class containing parametrized angular distributions
//              in the CM for gamma - nucleon two-body scattering

#include "G4GammaNuclAngDst.hh"

namespace {
  static const G4double gnke[10] =
    {0.0, 0.11, 0.22, 0.26, 0.30, 0.34, 0.42, 0.59, 0.79, 10.0};

  static const G4double gnA[10] = 
    {0.0,   0.0,   5.16,   5.55,  5.33,  7.40, 13.55,  13.44,  13.31,   7.3};
  static const G4double gnC[10] =
    {0.0, -10.33, -5.44,  -5.92, -4.27, -0.66,  1.37,   1.07,   0.52,   7.3};
  static const G4double gnCos[10] =
    {1.0,   1.0,   0.906,  0.940, 0.940, 0.906, 0.906,  0.91,   0.91,   0.94};
  static const G4double gnFrac[10] =
    {0.0,   0.0,   0.028,  0.012, 0.014, 0.044, 0.087,  0.122,  0.16,   1.0};
}

// Constructor passes arrays to templated base class

G4GammaNuclAngDst::G4GammaNuclAngDst(G4int verbose)
  : G4ParamExpTwoBodyAngDst<10>("G4GammaNuclAngDist", gnke, gnFrac,
			     gnA, gnC, gnCos, verbose) {;}
