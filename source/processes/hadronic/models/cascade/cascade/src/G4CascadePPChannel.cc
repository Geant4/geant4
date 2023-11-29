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

#include "G4CascadePPChannel.hh"
#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

namespace {
  // p p : Outgoing particle types of a given multiplicity
  static const G4int pp2bfs[1][2] =
  {{pro,pro}};

  static const G4int pp3bfs[6][3] =
  {{pro,pro,pi0}, {pro,neu,pip}, {pro,kpl,lam},
   {pro,kpl,s0},  {pro,k0,sp},   {neu,kpl,sp}};

  static const G4int pp4bfs[18][4] =
  {{pro,pro,pip,pim}, {pro,neu,pip,pi0}, {pro,pro,pi0,pi0},
   {neu,neu,pip,pip}, {pro,pi0,kpl,lam}, {pro,pip,k0,lam},
   {neu,pip,kpl,lam}, {neu,pip,kpl,s0},  {pro,pi0,kpl,s0},
   {pro,pip,k0,s0},   {pro,pip,kpl,sm},  {pro,pi0,k0,sp},
   {neu,pip,k0,sp},   {pro,pim,kpl,sp},  {neu,pi0,kpl,sp},
   {pro,pro,kpl,kmi}, {pro,pro,k0,k0b},  {pro,neu,kpl,k0b}};

  static const G4int pp5bfs[32][5] =
  {{pro,pro,pip,pim,pi0}, {pro,pro,pi0,pi0,pi0}, {pro,neu,pip,pip,pim},
   {pro,neu,pip,pi0,pi0}, {neu,neu,pip,pip,pi0}, {pro,pip,pim,kpl,lam},
   {pro,pi0,pi0,kpl,lam}, {pro,pip,pi0,k0,lam},  {neu,pip,pi0,kpl,lam},
   {neu,pip,pip,k0,lam},  {pro,pip,pim,kpl,s0},  {pro,pi0,pi0,kpl,s0},
   {pro,pip,pi0,k0,s0},   {neu,pip,pi0,kpl,s0},  {neu,pip,pip,k0,s0},
   {pro,pip,pim,k0,sp},   {pro,pi0,pi0,k0,sp},   {pro,pim,pi0,kpl,sp},
   {neu,pip,pi0,k0,sp},   {neu,sp,kpl,pip,pim},  {neu,pi0,pi0,kpl,sp},  
   {pro,pip,pi0,kpl,sm},  {pro,pip,pip,k0,sm},   {neu,pip,pip,kpl,sm},
   {pro,pro,pip,k0,kmi},  {pro,pro,pim,kpl,k0b}, {pro,pro,pi0,k0,k0b},
   {pro,pro,pi0,kpl,kmi}, {pro,neu,pip,k0,k0b},  {pro,neu,pip,kpl,kmi},
   {pro,neu,pi0,kpl,k0b}, {neu,neu,pip,kpl,k0b}};

  static const G4int pp6bfs[48][6] =
  {{pro,pro,pip,pip,pim,pim}, {pro,pro,pip,pim,pi0,pi0},
   {pro,pro,pi0,pi0,pi0,pi0}, {pro,neu,pip,pip,pim,pi0},
   {pro,neu,pip,pi0,pi0,pi0}, {neu,neu,pip,pip,pip,pim},
   {neu,neu,pip,pip,pi0,pi0}, {pro,pip,pim,pi0,kpl,lam},
   {pro,pi0,pi0,pi0,kpl,lam}, {pro,pip,pip,pim,k0,lam},
   {pro,pip,pi0,pi0,k0,lam},  {neu,pip,pip,pim,kpl,lam},
   {neu,pip,pi0,pi0,kpl,lam}, {neu,pip,pip,pi0,k0,lam}, 
   {pro,pip,pim,pi0,kpl,s0},  {pro,pi0,pi0,pi0,kpl,s0}, 
   {pro,pip,pip,pim,k0,s0},   {pro,pip,pi0,pi0,k0,s0}, 
   {neu,pip,pip,pim,kpl,s0},  {neu,pip,pi0,pi0,kpl,s0}, 
   {neu,pip,pip,pi0,k0,s0},   {pro,pim,pi0,pi0,kpl,sp},
   {pro,pip,pim,pim,kpl,sp},  {pro,pi0,pi0,pi0,k0,sp},
   {pro,pip,pim,pi0,k0,sp},   {neu,pi0,pi0,pi0,kpl,sp},
   {neu,pip,pim,pi0,kpl,sp},  {neu,pip,pi0,pi0,k0,sp},
   {neu,pip,pip,pim,k0,sp},   {pro,pip,pi0,pi0,kpl,sm},
   {pro,pip,pip,pim,kpl,sm},  {pro,pip,pip,pi0,k0,sm},
   {neu,pip,pip,pi0,kpl,sm},  {neu,pip,pip,pip,k0,sm}, 
   {pro,pro,pim,pi0,kpl,k0b}, {pro,pro,pi0,pi0,kpl,kmi},
   {pro,pro,pip,pim,kpl,kmi}, {pro,pro,pi0,pi0,k0,k0b}, 
   {pro,pro,pip,pim,k0,k0b},  {pro,pro,pip,pi0,kmi,k0},
   {pro,neu,pi0,pi0,kpl,k0b}, {pro,neu,pip,pim,kpl,k0b},
   {pro,neu,pip,pi0,kpl,kmi}, {pro,neu,pip,pi0,k0,k0b},
   {pro,neu,pip,pip,kmi,k0},  {neu,neu,pip,pi0,kpl,k0b},
   {neu,neu,pip,pip,kpl,kmi}, {neu,neu,pip,pip,k0,k0b}};

  static const G4int pp7bfs[63][7] =
  {{pro,pro,pip,pip,pim,pim,pi0}, {pro,pro,pip,pim,pi0,pi0,pi0},
   {pro,pro,pi0,pi0,pi0,pi0,pi0}, {pro,neu,pip,pip,pip,pim,pim},
   {pro,neu,pip,pip,pim,pi0,pi0}, {pro,neu,pip,pi0,pi0,pi0,pi0},
   {neu,neu,pip,pip,pip,pim,pi0}, {neu,neu,pip,pip,pi0,pi0,pi0},
   {pro,pip,pip,pim,pim,kpl,lam}, {pro,pip,pim,pi0,pi0,kpl,lam},
   {pro,pi0,pi0,pi0,pi0,kpl,lam}, {pro,pip,pip,pim,pi0,k0,lam},
   {pro,pip,pi0,pi0,pi0,k0,lam},  {neu,pip,pip,pim,pi0,kpl,lam},
   {neu,pip,pi0,pi0,pi0,kpl,lam}, {neu,pip,pip,pi0,pi0,k0,lam},
   {neu,pip,pip,pip,pim,k0,lam},  {pro,pip,pip,pim,pim,kpl,s0},
   {pro,pip,pim,pi0,pi0,kpl,s0},  {pro,pi0,pi0,pi0,pi0,kpl,s0},
   {pro,pip,pip,pim,pi0,k0,s0},   {pro,pip,pi0,pi0,pi0,k0,s0},
   {neu,pip,pip,pim,pi0,kpl,s0},  {neu,pip,pi0,pi0,pi0,kpl,s0},
   {neu,pip,pip,pi0,pi0,k0,s0},   {neu,pip,pip,pip,pim,k0,s0},
   {pro,pip,pim,pim,pi0,kpl,sp},  {pro,pim,pi0,pi0,pi0,kpl,sp},
   {pro,pip,pip,pim,pim,k0,sp},   {pro,pip,pim,pi0,pi0,k0,sp},
   {pro,pi0,pi0,pi0,pi0,k0,sp},   {neu,pip,pip,pim,pim,kpl,sp},
   {neu,pip,pim,pi0,pi0,kpl,sp},  {neu,pi0,pi0,pi0,pi0,kpl,sp},
   {neu,pip,pip,pim,pi0,k0,sp},   {neu,pip,pi0,pi0,pi0,k0,sp},
   {pro,pip,pip,pim,pi0,kpl,sm},  {pro,pip,pi0,pi0,pi0,kpl,sm},
   {pro,pip,pip,pip,pim,k0,sm},   {pro,pip,pip,pi0,pi0,k0,sm},
   {neu,pip,pip,pip,pim,kpl,sm},  {neu,pip,pip,pi0,pi0,kpl,sm},
   {neu,pip,pip,pip,pi0,k0,sm},   {pro,pro,pip,pim,pim,kpl,k0b},
   {pro,pro,pim,pi0,pi0,kpl,k0b}, {pro,pro,pip,pim,pi0,kpl,kmi},
   {pro,pro,pi0,pi0,pi0,kpl,kmi}, {pro,pro,pip,pim,pi0,k0,k0b},
   {pro,pro,pi0,pi0,pi0,k0,k0b},  {pro,pro,pip,pip,pim,kmi,k0}, 
   {pro,pro,pip,pi0,pi0,kmi,k0},  {pro,neu,pip,pim,pi0,kpl,k0b},
   {pro,neu,pi0,pi0,pi0,kpl,k0b}, {pro,neu,pip,pip,pim,kpl,kmi},
   {pro,neu,pip,pi0,pi0,kpl,kmi}, {pro,neu,pip,pip,pim,k0,k0b},
   {pro,neu,pip,pi0,pi0,k0,k0b},  {pro,neu,pip,pip,pi0,kmi,k0},
   {neu,neu,pip,pip,pim,kpl,k0b}, {neu,neu,pip,pi0,pi0,kpl,k0b},
   {neu,neu,pip,pip,pi0,kpl,kmi}, {neu,neu,pip,pip,pi0,k0,k0b},
   {neu,neu,pip,pip,pip,kmi,k0}};
 
  static const G4int pp8bfs[73][8] =
  {{pro,pro,pip,pip,pip,pim,pim,pim}, {pro,pro,pip,pip,pim,pim,pi0,pi0},
   {pro,pro,pip,pim,pi0,pi0,pi0,pi0}, {pro,pro,pi0,pi0,pi0,pi0,pi0,pi0},
   {pro,neu,pip,pip,pip,pim,pim,pi0}, {pro,neu,pip,pip,pim,pi0,pi0,pi0},
   {pro,neu,pip,pi0,pi0,pi0,pi0,pi0}, {neu,neu,pip,pip,pip,pip,pim,pim},
   {neu,neu,pip,pip,pip,pim,pi0,pi0}, {neu,neu,pip,pip,pi0,pi0,pi0,pi0},
   {pro,pip,pip,pim,pim,pi0,kpl,lam}, {pro,pip,pim,pi0,pi0,pi0,kpl,lam},
   {pro,pip,pip,pip,pim,pim,k0,lam},  {pro,pip,pip,pim,pi0,pi0,k0,lam},
   {pro,pip,pi0,pi0,pi0,pi0,k0,lam},  {neu,pip,pip,pip,pim,pim,kpl,lam},
   {neu,pip,pip,pim,pi0,pi0,kpl,lam}, {neu,pip,pi0,pi0,pi0,pi0,kpl,lam}, 
   {neu,pip,pip,pip,pim,pi0,k0,lam},  {neu,pip,pip,pi0,pi0,pi0,k0,lam}, 
   {pro,pip,pip,pim,pim,pi0,kpl,s0},  {pro,pip,pim,pi0,pi0,pi0,kpl,s0},
   {pro,pip,pip,pip,pim,pim,k0,s0},   {pro,pip,pip,pim,pi0,pi0,k0,s0}, 
   {pro,pip,pi0,pi0,pi0,pi0,k0,s0},   {neu,pip,pip,pip,pim,pim,kpl,s0},
   {neu,pip,pip,pim,pi0,pi0,kpl,s0},  {neu,pip,pi0,pi0,pi0,pi0,kpl,s0},
   {neu,pip,pip,pip,pim,pi0,k0,s0},   {neu,pip,pip,pi0,pi0,pi0,k0,s0},
   {pro,pip,pip,pim,pim,pim,kpl,sp},  {pro,pip,pim,pim,pi0,pi0,kpl,sp},
   {pro,pim,pi0,pi0,pi0,pi0,kpl,sp},  {pro,pip,pip,pim,pim,pi0,k0,sp},
   {pro,pip,pim,pi0,pi0,pi0,k0,sp},   {neu,pip,pip,pim,pim,pi0,kpl,sp},
   {neu,pip,pim,pi0,pi0,pi0,kpl,sp},  {neu,pip,pip,pip,pim,pim,k0,sp},
   {neu,pip,pip,pim,pi0,pi0,k0,sp},   {neu,pip,pi0,pi0,pi0,pi0,k0,sp},
   {pro,pip,pip,pip,pim,pim,kpl,sm},  {pro,pip,pip,pim,pi0,pi0,kpl,sm},
   {pro,pip,pi0,pi0,pi0,pi0,kpl,sm},  {pro,pip,pip,pip,pim,pi0,k0,sm},
   {pro,pip,pip,pi0,pi0,pi0,k0,sm},   {neu,pip,pip,pip,pim,pi0,kpl,sm},
   {neu,pip,pip,pi0,pi0,pi0,kpl,sm},  {neu,pip,pip,pip,pi0,pi0,k0,sm},
   {neu,pip,pip,pip,pip,pim,k0,sm},   {pro,pro,pip,pim,pim,pi0,kpl,k0b},
   {pro,pro,pim,pi0,pi0,pi0,kpl,k0b}, {pro,pro,pip,pim,pi0,pi0,kpl,kmi}, 
   {pro,pro,pip,pip,pim,pim,kpl,kmi}, {pro,pro,pip,pim,pi0,pi0,k0,k0b},
   {pro,pro,pip,pip,pim,pim,k0,k0b},  {pro,pro,pip,pip,pim,pi0,kmi,k0},
   {pro,pro,pip,pi0,pi0,pi0,kmi,k0},  {pro,neu,pip,pim,pi0,pi0,kpl,k0b},
   {pro,neu,pip,pip,pim,pim,kpl,k0b}, {pro,neu,pi0,pi0,pi0,pi0,kpl,k0b},
   {pro,neu,pip,pip,pim,pi0,kpl,kmi}, {pro,neu,pip,pi0,pi0,pi0,kpl,kmi},
   {pro,neu,pip,pip,pim,pi0,k0,k0b},  {pro,neu,pip,pi0,pi0,pi0,k0,k0b},
   {pro,neu,pip,pip,pi0,pi0,kmi,k0},  {pro,neu,pip,pip,pip,pim,kmi,k0},
   {neu,neu,pip,pip,pim,pi0,kpl,k0b}, {neu,neu,pip,pi0,pi0,pi0,kpl,k0b},
   {neu,neu,pip,pip,pi0,pi0,kpl,kmi}, {neu,neu,pip,pip,pip,pim,kpl,kmi},
   {neu,neu,pip,pip,pi0,pi0,k0,k0b},  {neu,neu,pip,pip,pip,pim,k0,k0b},
   {neu,neu,pip,pip,pip,pi0,kmi,k0}};

  static const G4int pp9bfs[79][9] =
  {{pro,pro,pip,pip,pip,pim,pim,pim,pi0}, {pro,pro,pip,pip,pim,pim,pi0,pi0,pi0},
   {pro,pro,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pro,neu,pip,pip,pip,pip,pim,pim,pim},
   {pro,neu,pip,pip,pip,pim,pim,pi0,pi0}, {pro,neu,pip,pip,pim,pi0,pi0,pi0,pi0},
   {pro,neu,pip,pi0,pi0,pi0,pi0,pi0,pi0}, {neu,neu,pip,pip,pip,pip,pim,pim,pi0},
   {neu,neu,pip,pip,pip,pim,pi0,pi0,pi0}, {neu,neu,pip,pip,pi0,pi0,pi0,pi0,pi0},
   {pro,pip,pip,pip,pim,pim,pim,kpl,lam}, {pro,pip,pip,pim,pim,pi0,pi0,kpl,lam},
   {pro,pip,pim,pi0,pi0,pi0,pi0,kpl,lam}, {pro,pip,pip,pip,pim,pim,pi0,k0,lam}, 
   {pro,pip,pip,pim,pi0,pi0,pi0,k0,lam},  {neu,pip,pip,pip,pim,pim,pi0,kpl,lam},
   {neu,pip,pip,pim,pi0,pi0,pi0,kpl,lam}, {neu,pip,pip,pip,pip,pim,pim,k0,lam},
   {neu,pip,pip,pip,pim,pi0,pi0,k0,lam},  {neu,pip,pip,pi0,pi0,pi0,pi0,k0,lam},
   {pro,pip,pip,pip,pim,pim,pim,kpl,s0},  {pro,pip,pip,pim,pim,pi0,pi0,kpl,s0},
   {pro,pip,pim,pi0,pi0,pi0,pi0,kpl,s0},  {pro,pip,pip,pip,pim,pim,pi0,k0,s0},
   {pro,pip,pip,pim,pi0,pi0,pi0,k0,s0},   {neu,pip,pip,pip,pim,pim,pi0,kpl,s0},
   {neu,pip,pip,pim,pi0,pi0,pi0,kpl,s0},  {neu,pip,pip,pip,pip,pim,pim,k0,s0},
   {neu,pip,pip,pip,pim,pi0,pi0,k0,s0},   {neu,pip,pip,pi0,pi0,pi0,pi0,k0,s0},
   {pro,pip,pip,pim,pim,pim,pi0,kpl,sp},  {pro,pip,pim,pim,pi0,pi0,pi0,kpl,sp},
   {pro,pip,pip,pip,pim,pim,pim,k0,sp},   {pro,pip,pip,pim,pim,pi0,pi0,k0,sp},
   {pro,pip,pim,pi0,pi0,pi0,pi0,k0,sp},   {neu,pip,pip,pip,pim,pim,pim,kpl,sp},
   {neu,pip,pip,pim,pim,pi0,pi0,kpl,sp},  {neu,pip,pim,pi0,pi0,pi0,pi0,kpl,sp},
   {neu,pip,pip,pip,pim,pim,pi0,k0,sp},   {neu,pip,pip,pim,pi0,pi0,pi0,k0,sp},
   {pro,pip,pip,pip,pim,pim,pi0,kpl,sm},  {pro,pip,pip,pim,pi0,pi0,pi0,kpl,sm},
   {pro,pip,pip,pip,pip,pim,pim,k0,sm},   {pro,pip,pip,pip,pim,pi0,pi0,k0,sm}, 
   {pro,pip,pip,pi0,pi0,pi0,pi0,k0,sm},   {neu,pip,pip,pip,pip,pim,pim,kpl,sm},
   {neu,pip,pip,pip,pim,pi0,pi0,kpl,sm},  {neu,pip,pip,pi0,pi0,pi0,pi0,kpl,sm},
   {neu,pip,pip,pip,pip,pim,pi0,k0,sm},   {neu,pip,pip,pip,pi0,pi0,pi0,k0,sm},
   {pro,pro,pip,pip,pim,pim,pim,kpl,k0b}, {pro,pro,pip,pim,pim,pi0,pi0,kpl,k0b},
   {pro,pro,pim,pi0,pi0,pi0,pi0,kpl,k0b}, {pro,pro,pip,pip,pim,pim,pi0,kpl,kmi},
   {pro,pro,pip,pim,pi0,pi0,pi0,kpl,kmi}, {pro,pro,pip,pip,pim,pim,pi0,k0,k0b},
   {pro,pro,pip,pim,pi0,pi0,pi0,k0,k0b},  {pro,pro,pip,pip,pip,pim,pim,kmi,k0}, 
   {pro,pro,pip,pip,pim,pi0,pi0,kmi,k0},  {pro,pro,pip,pi0,pi0,pi0,pi0,kmi,k0},
   {pro,neu,pip,pip,pim,pim,pi0,kpl,k0b}, {pro,neu,pip,pim,pi0,pi0,pi0,kpl,k0b},
   {pro,neu,pip,pip,pip,pim,pim,kpl,kmi}, {pro,neu,pip,pip,pim,pi0,pi0,kpl,kmi},
   {pro,neu,pip,pi0,pi0,pi0,pi0,kpl,kmi}, {pro,neu,pip,pip,pip,pim,pim,k0,k0b},
   {pro,neu,pip,pip,pim,pi0,pi0,k0,k0b},  {pro,neu,pip,pi0,pi0,pi0,pi0,k0,k0b},
   {pro,neu,pip,pip,pip,pim,pi0,kmi,k0},  {pro,neu,pip,pip,pi0,pi0,pi0,kmi,k0},
   {neu,neu,pip,pip,pip,pim,pim,kpl,k0b}, {neu,neu,pip,pip,pim,pi0,pi0,kpl,k0b},
   {neu,neu,pip,pi0,pi0,pi0,pi0,kpl,k0b}, {neu,neu,pip,pip,pip,pim,pi0,kpl,kmi},
   {neu,neu,pip,pip,pi0,pi0,pi0,kpl,kmi}, {neu,neu,pip,pip,pip,pim,pi0,k0,k0b},
   {neu,neu,pip,pip,pi0,pi0,pi0,k0,k0b},  {neu,neu,pip,pip,pip,pip,pim,kmi,k0},
   {neu,neu,pip,pip,pip,pi0,pi0,kmi,k0}};

}

namespace {
  // Total p p cross sections as a function of kinetic energy
  static const G4double ppTotXSec[30] = 
  // Stepanov cross sections below 130 MeV
   {17613.0,  863.3,  674.6,  495.2, 376.0, 285.4, 205.8,  135.7,   93.7,   69.1,
       56.0,   46.0,   40.0,   35.6,  33.0,  34.9,  44.515, 46.855, 44.868, 46.0,
       44.012, 41.838, 41.177, 40.65, 40.22, 40.0,  39.26,  38.36,  38.37,  38.41};


  static const G4double ppCrossSections[320][30] = {
  //
  // multiplicity 2  index 0   (1 channel)
  //
  //  p p
  // Stepanov cross sections below 130 MeV 
   {17613.0, 863.3, 674.6, 495.2, 376.0, 285.4, 205.8, 135.7, 93.7, 69.1,
       56.0,  46.0,  40.0,  35.6,  32.25, 28.7,  26.0,  23.2, 20.7, 18.0,
       15.7,  14.0,  12.5,  11.2,  10.1,   9.4,   8.9,   8.4,  8.0,  7.7},
  //
  // multiplicity 3  index 1 - 6  (6 channels)
  //
  //  p p pi0 (n n pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.15, 1.1, 3.5,  4.3,  4.03, 3.9,
     3.6, 3.1, 2.8, 2.3, 1.8,  1.5, 1.25, 1.05, 0.9,  0.8},

  //  p n pi+ (p n pi-)
   { 0.0,  0.0,  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0, 0.0, 0.6, 5.1, 15.0, 19.1, 18.0, 16.5,
     13.3, 10.0, 7.8, 5.9, 4.4, 3.4,  2.6,  2.0,  1.65, 1.35},

  //  p L K+ (n L K0)
   { 0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.03, 0.05, 0.058, 0.055, 0.05, 0.044, 0.04, 0.036, 0.033, 0.03},

  //  p S0 K+ (n S0 K0)
   { 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.015, 0.025, 0.029, 0.027, 0.025, 0.022, 0.02, 0.018, 0.016, 0.015},

  //  p S+ K0 (n S- K+)
   { 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.015, 0.025, 0.029, 0.027, 0.025, 0.022, 0.02, 0.018, 0.016, 0.014},

  //  n S+ K+ (p S- K0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.03, 0.06, 0.07, 0.06, 0.048, 0.04, 0.034, 0.028, 0.024, 0.02},
  //
  // multiplicity 4  index 7 - 24  (18 channels)
  //
  //  p p pi+ pi- (n n pi+ pi-) (391)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.005, 0.085, 0.6, 1.9,
     2.8, 3.0, 3.0, 2.8, 2.5, 2.1, 1.9,   1.6,   1.4, 1.2},

  //  p n pi+ pi0 (p n pi- pi0)  (332)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,  0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.006, 0.12, 1.0, 3.5,
     4.2, 3.9, 3.5, 3.1, 2.8, 2.4, 2.2,   1.9,  1.7, 1.5},

  //  p p pi0 pi0 (n n pi0 pi0)  (416)
   { 0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0,   0.0,   0.0,  0.0,
     0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.002, 0.025, 0.24, 0.76,
     1.12, 1.2, 1.2, 1.12, 1.0, 0.87, 0.76,  0.65,  0.56, 0.48},

  //  n n pi+ pi+ (p p pi- pi-)
   { 0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0,   0.0,   0.0,  0.0,
     0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.002, 0.025, 0.24, 0.76,
     1.12, 1.2, 1.2, 1.12, 1.0, 0.87, 0.76,  0.65,  0.56, 0.48},

  //  L K+ p pi0 (L K0 n pi0)  267
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.033, 0.055, 0.06, 0.055, 0.05, 0.045, 0.04, 0.035},

  //  L K0 p pi+ (L K+ n pi-)    263
   { 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,
     0.0, 0.015, 0.05, 0.08, 0.085, 0.08, 0.07, 0.06, 0.052, 0.045},

  //  L K+ n pi+ (L K0 p pi-)   271
   { 0.0, 0.0,  0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,
     0.0, 0.0,  0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,
     0.0, 0.01, 0.04, 0.075, 0.075, 0.07, 0.06, 0.05, 0.042, 0.035},

  //  S0 K+ n pi+ (S0 K0 p pi-)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.005, 0.02, 0.037, 0.037, 0.035, 0.03, 0.025, 0.02, 0.017},

  //  S0 K+ p pi0 (S0 K0 n pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.005, 0.016, 0.027, 0.03, 0.027, 0.025, 0.022, 0.02, 0.017},

  //  S0 K0 p pi+ (S0 K+ n pi-)  303
   { 0.0, 0.0, 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.0, 0.01, 0.033, 0.04, 0.035, 0.03, 0.025, 0.02, 0.015},

  //  S- K+ p pi+ (S+ K0 n pi-)   314
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.003, 0.015, 0.035, 0.04, 0.035, 0.03, 0.024, 0.019, 0.015},

  //  S+ K0 p pi0 (S- K+ n pi0)  291
   { 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.001, 0.01, 0.02, 0.025, 0.025, 0.02, 0.015, 0.012, 0.009},

  //  S+ K0 n pi+ (S- K+ p pi-)  295
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.001, 0.007, 0.025, 0.035, 0.035, 0.03, 0.024, 0.019, 0.015},

  //  S+ K+ p pi- (S- K0 n pi+)   292
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.002, 0.012, 0.035, 0.04, 0.035, 0.03, 0.024, 0.019, 0.015},

  //  S+ K+ n pi0 (S- K0 p pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.002, 0.012, 0.035, 0.04, 0.035, 0.03, 0.024, 0.019, 0.015},

  //  p p K+ K- (n n K0 K0b)
   { 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.003, 0.01, 0.02, 0.025, 0.023, 0.019, 0.016, 0.013},

  //  p p K0 K0b (n n K+ K-)   (424)
   { 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.003, 0.01, 0.02, 0.025, 0.023, 0.019, 0.016, 0.013},

  //  p n K+ K0b (p n K0 K-)    (344)
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.001, 0.012, 0.019, 0.02, 0.018, 0.014, 0.01, 0.008, 0.006},
  //
  // multiplicity 5  index 25 - 56  (32 channels)
  //
  //  p p pi+ pi- pi0 (n n pi+ pi- pi0) (385)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.007, 0.1,
     0.4,  1.1,  1.8,  2.2,  2.2,  2.2,  2.0,  1.7,  1.5,  1.3 },

  //  p p pi0 pi0 pi0 (n n pi0 pi0 pi0)
   { 0.0,   0.0,   0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.01,
     0.067, 0.183, 0.3, 0.367, 0.367, 0.367, 0.332, 0.283, 0.255, 0.217},

  //  p n pi+ pi+ pi- (p n pi+ pi- pi-)  (336)
   { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0,   0.0, 0.0,   0.0,
     0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0,   0.0, 0.024, 0.26,
     0.7, 1.6, 2.19, 2.3, 2.3, 2.06, 1.882, 1.6, 1.3,   1.0},

  //  p n pi+ pi0 pi0 (p n pi- pi0 pi0)
   { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0,   0.0, 0.0,   0.0,
     0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0,   0.0, 0.024, 0.26,
     0.7, 1.6, 2.19, 2.3, 2.3, 2.06, 1.882, 1.6, 1.3,   1.0},

  //  n n pi+ pi+ pi0 (p p pi- pi- pi0) 
   { 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,   0.0,
     0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.003, 0.05,
     0.2, 0.55, 0.9, 1.2, 1.2, 1.1, 1.0, 0.85, 0.758, 0.65},

  //  p L K+ pi+ pi- (n L K0 pi+ pi-)  (262)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.005, 0.025, 0.048, 0.05, 0.047, 0.04, 0.035, 0.03},

  //  p L K+ pi0 pi0 (n L K0 pi0 pi0) 
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.002, 0.012, 0.023, 0.025, 0.023, 0.02, 0.017, 0.015},

  //  p L K0 pi+ pi0 (n L K+ pi- pi0)  (261) 
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.005, 0.025, 0.048, 0.05, 0.047, 0.04, 0.035, 0.03},

  //  n L K+ pi+ pi0 (p L K0 pi- pi0) 
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.005, 0.025, 0.048, 0.05, 0.047, 0.04, 0.035, 0.03},

  //  n L K0 pi+ pi+ (p L K+ pi- pi-)  (274)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.001, 0.012, 0.033, 0.039, 0.035, 0.03, 0.025, 0.02},

  // p S0 K+ pi+ pi- (n S0 K0 pi+ pi-)  (302)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.002, 0.011, 0.023, 0.025, 0.023, 0.02, 0.017, 0.015},

  //  p S0 K+ pi0 pi0 (n S0 K0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.001, 0.006, 0.011, 0.012, 0.011, 0.01, 0.009, 0.008},

  //  p S0 K0 pi+ pi0 (n S0 K+ pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.002, 0.011, 0.023, 0.025, 0.023, 0.02, 0.017, 0.015},

  //  n S0 K+ pi+ pi0 (p S0 K0 pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.002, 0.011, 0.023, 0.025, 0.023, 0.02, 0.017, 0.015},

  //  n S0 K0 pi+ pi+ (p S0 K+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.006, 0.016, 0.02, 0.017, 0.015, 0.012, 0.01},

  //  p S+ K0 pi+ pi- (n S- K+ pi+ pi-)  (289)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.004, 0.02, 0.033, 0.037, 0.035, 0.03, 0.025},

  //  p S+ K0 pi0 pi0 (n S- K+ pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.016, 0.019, 0.017, 0.015, 0.013},

  //  p S+ K+ pi- pi0 (n S- K0 pi+ pi0) (290)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.004, 0.02, 0.033, 0.037, 0.035, 0.03, 0.025},

  //  n S+ K0 pi+ pi0 (p S- K+ pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.008, 0.04, 0.066, 0.074, 0.07, 0.06, 0.05},

  //  n S+ K+ pi+ pi- (p S- K0 pi+ pi-)  (294)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.008, 0.04, 0.066, 0.074, 0.07, 0.06, 0.05},

  //  n S+ K+ pi0 pi0 (p S- K0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.004, 0.02, 0.033, 0.037, 0.035, 0.03, 0.025},

  //  p S- K+ pi+ pi0 (n S+ K0 pi- pi0)  (313)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.005, 0.02, 0.03, 0.03, 0.026, 0.022, 0.018},

  //  p S- K0 pi+ pi+ (n S+ K+ pi- pi-)  (316)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.005, 0.02, 0.03, 0.03, 0.026, 0.022, 0.018},

  //  n S- K+ pi+ pi+ (p S+ K0 pi- pi-)  (317)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.009, 0.015, 0.015, 0.013, 0.011, 0.009},

  //  p p pi+ K0 K- (n n pi- K+ K0b) (399)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.005, 0.02, 0.027, 0.025, 0.022, 0.02, 0.017},

  //  p p pi- K+ K0b (n n pi+ K0 K-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.005, 0.02, 0.027, 0.025, 0.022, 0.02, 0.017},

  //  p p pi0 K0 K0b (n n pi0 K+ K-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.005, 0.02, 0.027, 0.025, 0.022, 0.02, 0.017},

  //  p p pi0 K+ K- (n n pi0 K0 K0bar)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.005, 0.02, 0.027, 0.025, 0.022, 0.02, 0.017},

  //  p n pi+ K0 K0b (p n pi- K+ K-)  (335)
   { 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.01, 0.04, 0.055, 0.05, 0.045, 0.039, 0.035},

  //  p n pi+ K+ K- (p n pi- K0 K0b)
   { 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.01, 0.04, 0.055, 0.05, 0.045, 0.039, 0.035},

  //  p n pi0 K+ K0b (p n pi0 K0 K-)
   { 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.01, 0.04, 0.055, 0.05, 0.045, 0.039, 0.035},

  // n n pi+ K+ K0b (p p pi- K0 K-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.005, 0.02, 0.027, 0.025, 0.022, 0.02, 0.017},
  //
  // multiplicity 6  index 57 - 104   (48 channels)
  //
  //  p p pi+ pi+ pi- pi- (n n pi+ pi+ pi- pi-)  (403)
   { 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.002, 0.018, 0.13, 0.33, 0.423, 0.461, 0.461, 0.461, 0.461, 0.461},

  //  p p pi+ pi- pi0 pi0 (n n pi+ pi- pi0 pi0)
   { 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.003, 0.024, 0.173, 0.45, 0.552, 0.614, 0.614, 0.614, 0.614, 0.614},

  //  p p pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.003, 0.022, 0.05, 0.066, 0.073, 0.073, 0.073, 0.073, 0.073},

  //  p n pi+ pi+ pi- pi0 (p n pi+ pi- pi- pi0)
   { 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.004, 0.036, 0.23, 0.55, 0.705, 0.833, 0.833, 0.833, 0.833, 0.833},

  //  p n pi+ pi0 pi0 pi0 (p n pi- pi0 pi0 pi0)
   { 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.003, 0.024, 0.173, 0.43, 0.552, 0.614, 0.614, 0.614, 0.614, 0.614},

  //  n n pi+ pi+ pi+ pi- (p p pi+ pi- pi- pi-)
   { 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.001, 0.012, 0.087, 0.227, 0.275, 0.307, 0.307, 0.307, 0.307, 0.307},

  //  n n pi+ pi+ pi0 pi0 (p p pi- pi- pi0 pi0)
   { 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.002, 0.018, 0.13, 0.34, 0.423, 0.461, 0.461, 0.461, 0.461, 0.461},

  //  p L K+ pi+ pi- pi0 (n L K0 pi+ pi- pi0)  (260)
   { 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.004, 0.02, 0.04, 0.046, 0.044, 0.044, 0.044, 0.044},

  //  p L K+ pi0 pi0 pi0 (n L K0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.008, 0.008, 0.008, 0.008},

  //  p L K0 pi+ pi+ pi- (n L K+ pi+ pi- pi-)  (266)
   { 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.01, 0.02, 0.023, 0.023, 0.023, 0.023, 0.023},

  //  p L K0 pi+ pi0 pi0 (n L K+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.01, 0.02, 0.023, 0.023, 0.023, 0.023, 0.023},

  //  n L K+ pi+ pi+ pi- (p L K0 pi+ pi- pi-)  (273)
   { 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.01, 0.02, 0.023, 0.023, 0.023, 0.023, 0.023},

  //  n L K+ pi+ pi0 pi0 (p L K0 pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.01, 0.02, 0.023, 0.023, 0.023, 0.023, 0.023},

  //  n L K0 pi+ pi+ pi0 (p L K+ pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.01, 0.02, 0.023, 0.023, 0.023, 0.023, 0.023},

  //  p S0 K+ pi+ pi- pi0 (n S0 K0 pi+ pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.005, 0.015, 0.022, 0.023, 0.023, 0.023, 0.023},

  //  p S0 K+ pi0 pi0 pi0 (n S0 K0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.005, 0.005, 0.005},

  //  p S0 K0 pi+ pi+ pi- (n S0 K+ pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.011, 0.011, 0.011, 0.011},

  //  p S0 K0 pi+ pi0 pi0 (n S0 K+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.011, 0.011, 0.011, 0.011},

  //  n S0 K+ pi+ pi+ pi- (p S0 K0 pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.011, 0.011, 0.011, 0.011},

  //  n S0 K+ pi+ pi0 pi0 (p S0 K0 pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.011, 0.011, 0.011, 0.011},

  //  n S0 K0 pi+ pi+ pi0 (p S0 K+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.011, 0.011, 0.011, 0.011},

  //  p S+ K+ pi- pi0 pi0 (n S- K0 pi+ pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.011, 0.011, 0.011, 0.011},

  //  p S+ K+ pi+ pi- pi- (n S- K0 pi+ pi+ pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.011, 0.011, 0.011, 0.011},

  //  p S+ K0 pi0 pi0 pi0 (n S- K+ pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  //  p S+ K0 pi+ pi- pi0 (n S- K+ pi+ pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.005, 0.015, 0.022, 0.023, 0.023, 0.023, 0.023},

  //  n S+ K+ pi0 pi0 pi0 (p S- K0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  //  n S+ K+ pi+ pi- pi0 (p S- K0 pi+ pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.005, 0.015, 0.022, 0.023, 0.023, 0.023, 0.023},

  //  n S+ K0 pi+ pi0 pi0 (p S- K+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.011, 0.011, 0.011, 0.011},

  //  n S+ K0 pi+ pi+ pi- (p S- K+ pi+ pi- pi-)  (296)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.011, 0.011, 0.011, 0.011},

  //  p S- K+ pi+ pi0 pi0 (n S+ K0 pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.013, 0.013, 0.013, 0.013, 0.013},

  //  p S- K+ pi+ pi+ pi- (n S+ K0 pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.013, 0.013, 0.013, 0.013, 0.013},

  //  p S- K0 pi+ pi+ pi0 (n S+ K+ pi- pi- pi0)  (315)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.013, 0.013, 0.013, 0.013, 0.013},

  //  n S- K+ pi+ pi+ pi0 (p S+ K0 pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.013, 0.013, 0.013, 0.013, 0.013},

  //  n S- K0 pi+ pi+ pi+ (p S+ K+ pi- pi- pi-)  (318)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  //  p p K+ K0b pi- pi0 (n n K- K0 pi+ pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.006, 0.01, 0.01, 0.01, 0.01, 0.01},

  //  p p K+ K- pi0 pi0 (n n K0 K0b pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.005, 0.005, 0.005, 0.005, 0.005},

  //  p p K+ K- pi+ pi- (n n K0 K0b pi+ pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.006, 0.01, 0.01, 0.01, 0.01, 0.01},

  //  p p K0 K0b pi0 pi0 (n n K+ K- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.005, 0.005, 0.005, 0.005, 0.005},

  //  p p K0 K0b pi+ pi- (n n K+ K- pi+ pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.006, 0.01, 0.01, 0.01, 0.01, 0.01},

  //  p p K0 K- pi+ pi0 (n n K+ K0b pi- pi0)  (389)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.006, 0.01, 0.01, 0.01, 0.01, 0.01},

  //  p n K+ K0b pi0 pi0 (p n K- K0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.007, 0.007, 0.007, 0.007, 0.007},

  //  p n K+ K0b pi+ pi- (p n K- K0 pi+ pi-)  (334)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.013, 0.013, 0.013, 0.013, 0.013},

  //  p n K+ K- pi+ pi0 (p n K0 K0b pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.013, 0.013, 0.013, 0.013, 0.013},

  //  p n K0 K0b pi+ pi0 (p n K+ K- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.013, 0.013, 0.013, 0.013, 0.013},

  //  p n K0 K- pi+ pi+ (p n K+ K0b pi- pi-) (339)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.007, 0.007, 0.007, 0.007, 0.007},

  //  n n K+ K0b pi+ pi0 (p p K- K0 pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.006, 0.007, 0.007, 0.007, 0.007},

  //  n n K+ K- pi+ pi+ (p p K0 K0b pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.003, 0.003, 0.003, 0.003, 0.003},

  //  n n K0 K0b pi+ pi+ (p p K+ K- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.003, 0.003, 0.003, 0.003, 0.003},
  //
  // multiplicity 7  index 105 - 167 (63 channels)
  //
  //  p p pi+ pi+ pi- pi- pi0 (n n pi+ pi+ pi- pi- pi0)  (400)
   { 0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0, 0.0,
     0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0, 0.0,
     0.0, 0.006, 0.045, 0.14, 0.32, 0.53, 0.681, 0.84, 1.0, 1.13},

  //  p p pi+ pi- pi0 pi0 pi0 (n n pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.004, 0.03, 0.093, 0.213, 0.352, 0.447, 0.564, 0.672, 0.761},

  //  p p pi0 pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.005, 0.01, 0.017, 0.022, 0.028, 0.033, 0.037},

  // p n pi+ pi+ pi+ pi- pi- (p n pi+ pi+ pi- pi- pi-) (340)
   { 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.005, 0.04, 0.12, 0.254, 0.367, 0.447, 0.54, 0.644, 0.729},

  //  p n pi+ pi+ pi- pi0 pi0 (p n pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0, 0.0,
     0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0, 0.0,
     0.0, 0.006, 0.045, 0.14, 0.32, 0.52, 0.681, 0.84, 1.0, 1.13},

  //  p n pi+ pi0 pi0 pi0 pi0 (p n pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.002, 0.02, 0.06, 0.127, 0.185, 0.224, 0.27, 0.322, 0.364},

  //  n n pi+ pi+ pi+ pi- pi0 (p p pi+ pi- pi- pi- pi0)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.004, 0.03, 0.093, 0.213, 0.35, 0.447, 0.54, 0.644, 0.729},

  //  n n pi+ pi+ pi0 pi0 pi0 (p p pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.002, 0.02, 0.06, 0.127, 0.185, 0.224, 0.27, 0.322, 0.364},

  //  p L K+ pi+ pi+ pi- pi- (n L K0 pi+ pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.014, 0.02, 0.028, 0.036, 0.044},

  //  p L K+ pi+ pi- pi0 pi0 (n L K0 pi+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.005, 0.015, 0.029, 0.039, 0.055, 0.072, 0.089},

  //  p L K+ pi0 pi0 pi0 pi0 (n L K0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p L K0 pi+ pi+ pi- pi0 (n L K+ pi+ pi- pi- pi0)  (265)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.005, 0.015, 0.029, 0.039, 0.055, 0.072, 0.089},

  //  p L K0 pi+ pi0 pi0 pi0 (n L K+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.005, 0.01, 0.014, 0.019, 0.024, 0.03},

  //  n L K+ pi+ pi+ pi- pi0 (p L K0 pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.006, 0.021, 0.038, 0.05, 0.065, 0.082, 0.097},

  //  n L K+ pi+ pi0 pi0 pi0 (p L K0 pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.013, 0.017, 0.022, 0.027, 0.032},

  //  n L K0 pi+ pi+ pi0 pi0 (p L K+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.01, 0.019, 0.025, 0.033, 0.041, 0.048},

  //  n L K0 pi+ pi+ pi+ pi- (p L K+ pi+ pi- pi- pi-)  (275)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.013, 0.017, 0.022, 0.027, 0.032},

  //  p S0 K+ pi+ pi+ pi- pi- (n S0 K0 pi+ pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.007, 0.01, 0.014, 0.018, 0.022},

  //  p S0 K+ pi+ pi- pi0 pi0 (n S0 K0 pi+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.014, 0.02, 0.027, 0.036, 0.044},

  //  p S0 K+ pi0 pi0 pi0 pi0 (n S0 K0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005},

  //  p S0 K0 pi+ pi+ pi- pi0 (n S0 K+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.014, 0.02, 0.027, 0.036, 0.044},

  //  p S0 K0 pi+ pi0 pi0 pi0 (n S0 K+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.006, 0.008, 0.009, 0.01},

  //  n S0 K+ pi+ pi+ pi- pi0 (p S0 K0 pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.011, 0.019, 0.024, 0.033, 0.041, 0.048},

  //  n S0 K+ pi+ pi0 pi0 pi0 (p S0 K0 pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.009, 0.01, 0.011},

  //  n S0 K0 pi+ pi+ pi0 pi0 (p S0 K+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.005, 0.01, 0.012, 0.016, 0.02, 0.024},

  //  n S0 K0 pi+ pi+ pi+ pi- (p S0 K+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.009, 0.01, 0.011},

  //  p S+ K+ pi+ pi- pi- pi0 (n S- K0 pi+ pi+ pi- pi0)   2
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.015, 0.02, 0.023, 0.025, 0.027},

  //  p S+ K+ pi- pi0 pi0 pi0 (n S- K0 pi+ pi0 pi0 pi0)   6
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.009, 0.01, 0.01},

  //  p S+ K0 pi+ pi+ pi- pi- (n S- K+ pi+ pi+ pi- pi-)  4
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.007, 0.01, 0.012, 0.013, 0.014},

  //  p S+ K0 pi+ pi- pi0 pi0 (n S- K+ pi+ pi- pi0 pi0)  2
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.015, 0.02, 0.023, 0.025, 0.027},

  //  p S+ K0 pi0 pi0 pi0 pi0 (n S- K+ pi0 pi0 pi0 pi0)   24
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003},

  //  n S+ K+ pi+ pi+ pi- pi- (p S- K0 pi+ pi+ pi- pi-)   4
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.007, 0.01, 0.012, 0.014, 0.014},

  //  n S+ K+ pi+ pi- pi0 pi0 (p S- K0 pi+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.014, 0.02, 0.023, 0.025, 0.027},

  //  n S+ K+ pi0 pi0 pi0 pi0 (p S- K0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003},

  //  n S+ K0 pi+ pi+ pi- pi0 (p S- K+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.014, 0.02, 0.023, 0.025, 0.027},

  //  n S+ K0 pi+ pi0 pi0 pi0 (p S- K+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.009, 0.01, 0.01},

  //  p S- K+ pi+ pi+ pi- pi0 (n S+ K0 pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.014, 0.02, 0.023, 0.025, 0.027},

  //  p S- K+ pi+ pi0 pi0 pi0 (n S+ K0 pi- pi0 pi0 pi0) 
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.009, 0.01, 0.01},

  //  p S- K0 pi+ pi+ pi+ pi- (n S+ K+ pi+ pi- pi- pi-) 
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.009, 0.01, 0.01},

  //  p S- K0 pi+ pi+ pi0 pi0 (n S+ K+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.007, 0.01, 0.012, 0.013, 0.014},

  //  n S- K+ pi+ pi+ pi+ pi- (p S+ K0 pi+ pi- pi- pi-) 
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.009, 0.01, 0.01},

  //  n S- K+ pi+ pi+ pi0 pi0 (p S+ K0 pi- pi- pi0 pi0) 
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.007, 0.01, 0.012, 0.013, 0.014},

  //  n S- K0 pi+ pi+ pi+ pi0 (p S+ K+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.009, 0.01, 0.01},

  //  p p K+ K0b pi+ pi- pi- (n n K- K0 pi+ pi+ pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.005, 0.005},

  //  p p K+ K0b pi- pi0 pi0 (n n K- K0 pi+ pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.005, 0.005},

  //  p p K+ K- pi+ pi- pi0 (n n K0 K0b pi+ pi- pi0) 
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.007, 0.009, 0.01, 0.01, 0.01},

  //  p p K+ K- pi0 pi0 pi0 (n n K0 K0b pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002},

  //  p p K0 K0b pi+ pi- pi0 (n n K+ K- pi+ pi- pi0) (388)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.007, 0.009, 0.01, 0.01, 0.01},

  //  p p K0 K0b pi0 pi0 pi0 (n n K+ K- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002},

  //  p p K0 K- pi+ pi+ pi- (n n K+ K0b pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.005, 0.005},

  //  p p K0 K- pi+ pi0 pi0 (n n K+ K0b pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.005, 0.005},

  //  p n K+ K0b pi+ pi- pi0 (p n K- K0 pi+ pi- pi0) 
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.006, 0.024, 0.05, 0.064, 0.072, 0.08, 0.08},

  //  p n K+ K0b pi0 pi0 pi0 (p n K- K0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.008, 0.01, 0.012, 0.013, 0.013},

  //  p n K+ K- pi+ pi+ pi- (p n K0 K0b pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.003, 0.012, 0.025, 0.032, 0.036, 0.04, 0.04},

  //  p n K+ K- pi+ pi0 pi0 (p n K0 K0b pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.003, 0.012, 0.025, 0.032, 0.036, 0.04, 0.04},

  //  p n K0 K0b pi+ pi+ pi- (p n K+ K- pi+ pi- pi-) (338)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.003, 0.012, 0.025, 0.032, 0.036, 0.04, 0.04},

  //  p n K0 K0b pi+ pi0 pi0 (p n K+ K- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.003, 0.012, 0.025, 0.032, 0.036, 0.04, 0.04},

  //  p n K0 K- pi+ pi+ pi0 (p n K+ K0b pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.006, 0.013, 0.02, 0.024, 0.026, 0.028},

  //  n n K+ K0b pi+ pi+ pi- (p p K- K0 pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.007, 0.01, 0.012, 0.013, 0.014},

  //  n n K+ K0b pi+ pi0 pi0 (p p K- K0 pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.007, 0.01, 0.012, 0.013, 0.014},

  //  n n K+ K- pi+ pi+ pi0 (p p K0 K0b pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.007, 0.01, 0.012, 0.013, 0.014},

  //  n n K0 K0b pi+ pi+ pi0 (p p K+ K- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.007, 0.01, 0.012, 0.013, 0.014},

  //  n n K0 K- pi+ pi+ pi+ (p p K+ K0b pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004}, 
  //
  // multiplicity 8  index 168 - 240   (73 channels)
  //
  //  p p pi+ pi+ pi+ pi- pi- pi- (n n pi+ pi+ pi+ pi- pi- pi-) (409)
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.004, 0.012, 0.03, 0.05, 0.075, 0.098, 0.122, 0.146},

  //  p p pi+ pi+ pi- pi- pi0 pi0 (n n pi+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.004, 0.018, 0.058, 0.135, 0.24, 0.338, 0.443, 0.546, 0.65},

  //  p p pi+ pi- pi0 pi0 pi0 pi0 (n n pi+ pi- pi0 pi0 pi0 pi0) 
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.002, 0.006, 0.019, 0.045, 0.082, 0.112, 0.148, 0.19, 0.219},

  //  p p pi0 pi0 pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007},

  //  p n pi+ pi+ pi+ pi- pi- pi0 (p n pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.004, 0.024, 0.076, 0.18, 0.32, 0.448, 0.59, 0.732, 0.874},

  //  p n pi+ pi+ pi- pi0 pi0 pi0 (p n pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.004, 0.024, 0.076, 0.18, 0.328, 0.448, 0.59, 0.732, 0.875},

  //  p n pi+ pi0 pi0 pi0 pi0 pi0 (p n pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.007, 0.018, 0.033, 0.045, 0.059, 0.073, 0.088},

  //  n n pi+ pi+ pi+ pi+ pi- pi- (p p pi+ pi+ pi- pi- pi- pi-)
   { 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.003, 0.01, 0.022, 0.041, 0.056, 0.074, 0.092, 0.11},

  //  n n pi+ pi+ pi+ pi- pi0 pi0 (p p pi+ pi- pi- pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.002, 0.012, 0.038, 0.09, 0.164, 0.224, 0.295, 0.366, 0.437},

  //  n n pi+ pi+ pi0 pi0 pi0 pi0 (p p pi- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.003, 0.01, 0.022, 0.041, 0.056, 0.074, 0.092, 0.11},

  //  p L K+ pi+ pi+ pi- pi- pi0 (n L K0 pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.016, 0.024, 0.034, 0.043, 0.053},

  //  p L K+ pi+ pi- pi0 pi0 pi0 (n L K0 pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.005, 0.01, 0.016, 0.022, 0.029, 0.035},

  //  p L K+ pi0 pi0 pi0 pi0 pi0 (n L K0 pi0 pi0 pi0 pi0 pi0) negligible

  //  p L K0 pi+ pi+ pi+ pi- pi- (n L K+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.005, 0.008, 0.011, 0.014, 0.017},
  
  //  p L K0 pi+ pi+ pi- pi0 pi0 (n L K+ pi+ pi- pi- pi0 pi0) 
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.016, 0.024, 0.034, 0.043, 0.053},

  //  p L K0 pi+ pi0 pi0 pi0 pi0 (n L K+ pi- pi0 pi0 pi0 pi0) 
   { 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.006, 0.007, 0.008},

  //  n L K+ pi+ pi+ pi+ pi- pi- (p L K0 pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.005, 0.008, 0.011, 0.014, 0.017},

  //  n L K+ pi+ pi+ pi- pi0 pi0 (p L K0 pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.016, 0.024, 0.034, 0.043, 0.053},

  //  n L K+ pi+ pi0 pi0 pi0 pi0 (p L K0 pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.006, 0.007, 0.008},

  //  n L K0 pi+ pi+ pi+ pi- pi0 (p L K+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.005, 0.01, 0.016, 0.022, 0.029, 0.035},

  //  n L K0 pi+ pi+ pi0 pi0 pi0 (p L K+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.005, 0.008, 0.011, 0.014, 0.017},

  //  p S0 K+ pi+ pi+ pi- pi- pi0 (n S0 K0 pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.008, 0.012, 0.017, 0.022, 0.026},

  //  p S0 K+ pi+ pi- pi0 pi0 pi0 (n S0 K0 pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.005, 0.008, 0.011, 0.015, 0.018},

  //  p S0 K+ pi0 pi0 pi0 pi0 pi0 (n S0 K0 pi0 pi0 pi0 pi0 pi0) negligible

  //  p S0 K0 pi+ pi+ pi+ pi- pi- (n S0 K+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p S0 K0 pi+ pi+ pi- pi0 pi0 (n S0 K+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.008, 0.012, 0.017, 0.022, 0.026},

  //  p S0 K0 pi+ pi0 pi0 pi0 pi0 (n S0 K+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005},

   // n S0 K+ pi+ pi+ pi+ pi- pi- (p S0 K0 pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  n S0 K+ pi+ pi+ pi- pi0 pi0 (p S0 K0 pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.008, 0.012, 0.017, 0.022, 0.026},

  //  n S0 K+ pi+ pi0 pi0 pi0 pi0 (p S0 K0 pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005},

  //  n S0 K0 pi+ pi+ pi+ pi- pi0 (p S0 K+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.005, 0.008, 0.011, 0.015, 0.018},

  //  n S0 K0 pi+ pi+ pi0 pi0 pi0 (p S0 K+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p S+ K+ pi+ pi+ pi- pi- pi- (n S- K0 pi+ pi+ pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  //  p S+ K+ pi+ pi- pi- pi0 pi0 (n S- K0 pi+ pi+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.005, 0.01, 0.015, 0.022, 0.029, 0.035},

  //  p S+ K+ pi- pi0 pi0 pi0 pi0 (n S- K0 pi+ pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p S+ K0 pi+ pi+ pi- pi- pi0 (n S- K+ pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.005, 0.01, 0.015, 0.022, 0.029, 0.035},

  //  p S+ K0 pi+ pi- pi0 pi0 pi0 (n S- K+ pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.007, 0.01, 0.013, 0.015, 0.017},

  //  p S+ K0 pi0 pi0 pi0 pi0 pi0 (n S- K+ pi0 pi0 pi0 pi0 pi0) negligible

  // n S+ K+ pi+ pi+ pi- pi- pi0 (p S- K0 pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.005, 0.01, 0.015, 0.022, 0.029, 0.035},

  //  n S+ K+ pi+ pi- pi0 pi0 pi0 (p S- K0 pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.007, 0.01, 0.013, 0.015, 0.017},

  //  n S+ K+ pi0 pi0 pi0 pi0 pi0 (p S- K0 pi0 pi0 pi0 pi0 pi0) negligible

  //  n S+ K0 pi+ pi+ pi+ pi- pi- (p S- K+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  // n S+ K0 pi+ pi+ pi- pi0 pi0 (p S- K+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.005, 0.01, 0.015, 0.022, 0.029, 0.035},

  //  n S+ K0 pi+ pi0 pi0 pi0 pi0 (p S- K+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p S- K+ pi+ pi+ pi+ pi- pi- (n S+ K0 pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  //  p S- K+ pi+ pi+ pi- pi0 pi0 (n S+ K0 pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.005, 0.01, 0.015, 0.022, 0.029, 0.035},

  //  p S- K+ pi+ pi0 pi0 pi0 pi0 (n S+ K0 pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p S- K0 pi+ pi+ pi+ pi- pi0 (n S+ K+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.007, 0.01, 0.013, 0.015, 0.017},

  //  p S- K0 pi+ pi+ pi0 pi0 pi0 (n S+ K+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  //  n S- K+ pi+ pi+ pi+ pi- pi0 (p S+ K0 pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.007, 0.01, 0.013, 0.015, 0.017},

  //  n S- K+ pi+ pi+ pi0 pi0 pi0 (p S+ K0 pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  //  n S- K0 pi+ pi+ pi+ pi0 pi0 (p S+ K+ pi- pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  //  n S- K0 pi+ pi+ pi+ pi+ pi- (p S+ K+ pi+ pi- pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p p K+ K0b pi+ pi- pi- pi0 (n n K- K0 pi+ pi+ pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.017, 0.021, 0.026},

  //  p p K+ K0b pi- pi0 pi0 pi0 (n n K- K0 pi+ pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.006, 0.007},

  //  p p K+ K- pi+ pi- pi0 pi0 (n n K0 K0b pi+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.017, 0.021, 0.026},

  //  p p K+ K- pi+ pi+ pi- pi- (n n K0 K0b pi+ pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012},

  //  p p K+ K- pi0 pi0 pi0 pi0 (n n K0 K0b pi0 pi0 pi0 pi0) negligible

  //  p p K0 K0b pi+ pi- pi0 pi0 (n n K+ K- pi+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.017, 0.021, 0.026},

  //  p p K0 K0b pi+ pi+ pi- pi- (n n K+ K- pi+ pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012},

  //  p p K0 K0b pi0 pi0 pi0 pi0 (n n K+ K- pi0 pi0 pi0 pi0) negligible

  //  p p K0 K- pi+ pi+ pi- pi0 (n n K+ K0b pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.017, 0.021, 0.026},

  //  p p K0 K- pi+ pi0 pi0 pi0 (n n K+ K0b pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.006, 0.007},

  //  p n K+ K0b pi+ pi- pi0 pi0 (p n K- K0 pi+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.018, 0.025, 0.032, 0.043, 0.052},

  //  p n K+ K0b pi+ pi+ pi- pi- (p n K- K0 pi+ pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.017, 0.021, 0.026},

  //  p n K+ K0b pi0 pi0 pi0 pi0 (p n K- K0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p n K+ K- pi+ pi+ pi- pi0 (p n K0 K0b pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.018, 0.025, 0.032, 0.043, 0.052},

  //  p n K+ K- pi+ pi0 pi0 pi0 (p n K0 K0b pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.01, 0.012, 0.013},

  //  p n K0 K0b pi+ pi+ pi- pi0 (p n K+ K- pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.017, 0.025, 0.032, 0.043, 0.052},

  //  p n K0 K0b pi+ pi0 pi0 pi0 (p n K+ K- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.01, 0.012, 0.013},

  //  p n K0 K- pi+ pi+ pi0 pi0 (p n K+ K0b pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.017, 0.021, 0.026},

  //  p n K0 K- pi+ pi+ pi+ pi- (p n K+ K0b pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.01, 0.012, 0.014},

  //  n n K+ K0b pi+ pi+ pi- pi0 (p p K- K0 pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.017, 0.021, 0.026},

  //  n n K+ K0b pi+ pi0 pi0 pi0 (p p K- K0 pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.006, 0.007},

  //  n n K+ K- pi+ pi+ pi0 pi0 (p p K0 K0b pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  //  n n K+ K- pi+ pi+ pi+ pi- (p p K0 K0b pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.006, 0.007},

  //  n n K0 K0b pi+ pi+ pi0 pi0 (p p K+ K- pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  //  n n K0 K0b pi+ pi+ pi+ pi- (p p K+ K- pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.006, 0.007},

  //  n n K0 K- pi+ pi+ pi+ pi0 (p p K+ K0b pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.006, 0.007},
  //
  // multiplicity 9   index 241 - 319  (79 channels)
  //
  //  p p pi+ pi+ pi+ pi- pi- pi- pi0 (n n pi+ pi+ pi+ pi- pi- pi- pi0) (407)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0, 0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0, 0.0,   0.0,
     0.0, 0.001, 0.005, 0.019, 0.055, 0.105, 0.15, 0.2, 0.257, 0.32},

  //  p p pi+ pi+ pi- pi- pi0 pi0 pi0 (n n pi+ pi+ pi- pi- pi0 pi0 pi0)
   { 0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
     0.0, 0.002, 0.009, 0.03, 0.087, 0.165, 0.23, 0.31, 0.39, 0.48},

  //  p p pi+ pi- pi0 pi0 pi0 pi0 pi0 (n n pi+ pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.002, 0.006, 0.018, 0.033, 0.044, 0.06, 0.077, 0.097},

  //  p p pi0 pi0 pi0 pi0 pi0 pi0 pi0 (n n pi0 pi0 pi0 pi0 pi0 pi0 pi0) negligible

  //  p n pi+ pi+ pi+ pi+ pi- pi- pi- (p n pi+ pi+ pi+ pi- pi- pi- pi-) (341)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0, 0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0, 0.0,
     0.0, 0.0, 0.002, 0.008, 0.025, 0.055, 0.09, 0.144, 0.2, 0.243},

  //  p n pi+ pi+ pi+ pi- pi- pi0 pi0 (p n pi+ pi+ pi- pi- pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,
     0.0, 0.003, 0.015, 0.06, 0.165, 0.31, 0.44, 0.576, 0.772, 0.96},

  //  p n pi+ pi+ pi- pi0 pi0 pi0 pi0 (p n pi+ pi- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0,  0.0,
     0.0, 0.001, 0.006, 0.024, 0.075, 0.16, 0.26, 0.4, 0.55, 0.7},

  //  p n pi+ pi0 pi0 pi0 pi0 pi0 pi0 (p n pi- pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.005, 0.009, 0.014, 0.02, 0.026, 0.032},

  //  n n pi+ pi+ pi+ pi+ pi- pi- pi0 (p p pi+ pi+ pi- pi- pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.001, 0.004, 0.015, 0.041, 0.082, 0.115, 0.15, 0.193, 0.243},

  //  n n pi+ pi+ pi+ pi- pi0 pi0 pi0 (p p pi+ pi- pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.001, 0.005, 0.02, 0.055, 0.105, 0.15, 0.206, 0.258, 0.324},

  //  n n pi+ pi+ pi0 pi0 pi0 pi0 pi0 (p p pi- pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.001, 0.006, 0.016, 0.033, 0.045, 0.06, 0.077, 0.097},

  //  p L K+ pi+ pi+ pi+ pi- pi- pi- (n L K0 pi+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p L K+ pi+ pi+ pi- pi- pi0 pi0 (n L K0 pi+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.009, 0.014, 0.021, 0.03, 0.038},

  //  p L K+ pi+ pi- pi0 pi0 pi0 pi0 (n L K0 pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.01, 0.013},

  //  p L K+ pi0 pi0 pi0 pi0 pi0 pi0 (n L K0 pi0 pi0 pi0 pi0 pi0 pi0) negligible

  //  p L K0 pi+ pi+ pi+ pi- pi- pi0 (n L K+ pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.006, 0.009, 0.014, 0.02, 0.025},

  //  p L K0 pi+ pi+ pi- pi0 pi0 pi0 (n L K+ pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.006, 0.009, 0.014, 0.02, 0.025},

  //  p L K0 pi+ pi0 pi0 pi0 pi0 pi0 (n L K+ pi- pi0 pi0 pi0 pi0 pi0) negligible

  //  n L K+ pi+ pi+ pi+ pi- pi- pi0 (p L K0 pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.006, 0.009, 0.014, 0.02, 0.025},

  //  n L K+ pi+ pi+ pi- pi0 pi0 pi0 (p L K0 pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.006, 0.009, 0.014, 0.02, 0.025},

  //  n L K+ pi+ pi0 pi0 pi0 pi0 pi0 (p L K0 pi- pi0 pi0 pi0 pi0 pi0) negligible

  //  n L K0 pi+ pi+ pi+ pi+ pi- pi- (p L K+ pi+ pi+ pi- pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  n L K0 pi+ pi+ pi+ pi- pi0 pi0 (p L K+ pi+ pi- pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.006, 0.009, 0.014, 0.02, 0.025},

  //  n L K0 pi+ pi+ pi0 pi0 pi0 pi0 (p L K+ pi- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p S0 K+ pi+ pi+ pi+ pi- pi- pi- (n S0 K0 pi+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005},

  //  p S0 K+ pi+ pi+ pi- pi- pi0 pi0 (n S0 K0 pi+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.007, 0.011, 0.015, 0.02},

  //  p S0 K+ pi+ pi- pi0 pi0 pi0 pi0 (n S0 K0 pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005},

  //  p S0 K+ pi0 pi0 pi0 pi0 pi0 pi0 (n S0 K0 pi0 pi0 pi0 pi0 pi0 pi0) negligible

  //  p S0 K0 pi+ pi+ pi+ pi- pi- pi0 (n S0 K+ pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011},

  //  p S0 K0 pi+ pi+ pi- pi0 pi0 pi0 (n S0 K+ pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011},

  //  p S0 K0 pi+ pi0 pi0 pi0 pi0 pi0 (n S0 K+ pi- pi0 pi0 pi0 pi0 pi0) negligible

  //  n S0 K+ pi+ pi+ pi+ pi- pi- pi0 (p S0 K0 pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011},

  //  n S0 K+ pi+ pi+ pi- pi0 pi0 pi0 (p S0 K0 pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011},

  //  n S0 K+ pi+ pi0 pi0 pi0 pi0 pi0 (p S0 K0 pi- pi0 pi0 pi0 pi0 pi0) negligible

  //  n S0 K0 pi+ pi+ pi+ pi+ pi- pi- (p S0 K+ pi+ pi+ pi- pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005},

  //  n S0 K0 pi+ pi+ pi+ pi- pi0 pi0 (p S0 K+ pi+ pi- pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.006, 0.01, 0.015, 0.02, 0.025},

  //  n S0 K0 pi+ pi+ pi0 pi0 pi0 pi0 (p S0 K+ pi- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005},

  //  p S+ K+ pi+ pi+ pi- pi- pi- pi0 (n S- K0 pi+ pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.011, 0.015, 0.018, 0.021},

  //  p S+ K+ pi+ pi- pi- pi0 pi0 pi0 (n S- K0 pi+ pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.011, 0.015, 0.018, 0.021},

  //  p S+ K+ pi- pi0 pi0 pi0 pi0 pi0 (n S- K0 pi+ pi0 pi0 pi0 pi0 pi0) negligible

  //  p S+ K0 pi+ pi+ pi+ pi- pi- pi- (n S- K+ pi+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p S+ K0 pi+ pi+ pi- pi- pi0 pi0 (n S- K+ pi+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.004, 0.011, 0.016, 0.021, 0.027, 0.031},

  //  p S+ K0 pi+ pi- pi0 pi0 pi0 pi0 (n S- K+ pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011},

  //  p S+ K0 pi0 pi0 pi0 pi0 pi0 pi0 (n S- K+ pi0 pi0 pi0 pi0 pi0 pi0) negligible

  //  n S+ K+ pi+ pi+ pi+ pi- pi- pi- (p S- K0 pi+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  n S+ K+ pi+ pi+ pi- pi- pi0 pi0 (p S- K0 pi+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.004, 0.011, 0.016, 0.021, 0.027, 0.031},

  //  n S+ K+ pi+ pi- pi0 pi0 pi0 pi0 (p S- K0 pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011},

  //  n S+ K+ pi0 pi0 pi0 pi0 pi0 pi0 (p S- K0 pi0 pi0 pi0 pi0 pi0 pi0) negligible

  //  n S+ K0 pi+ pi+ pi+ pi- pi- pi0 (p S- K+ pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.011, 0.015, 0.018, 0.021},

  //  n S+ K0 pi+ pi+ pi- pi0 pi0 pi0 (p S- K+ pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.011, 0.015, 0.018, 0.021},

  //  n S+ K0 pi+ pi0 pi0 pi0 pi0 pi0 (p S- K+ pi- pi0 pi0 pi0 pi0 pi0) negligible

  //  p S- K+ pi+ pi+ pi+ pi- pi- pi0 (n S+ K0 pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.011, 0.015, 0.018, 0.021},

  //  p S- K+ pi+ pi+ pi- pi0 pi0 pi0 (n S+ K0 pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.011, 0.015, 0.018, 0.021},

  //  p S- K+ pi+ pi0 pi0 pi0 pi0 pi0 (n S+ K0 pi- pi0 pi0 pi0 pi0 pi0) negligible

  //  p S- K0 pi+ pi+ pi+ pi+ pi- pi- (n S+ K+ pi+ pi+ pi- pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p S- K0 pi+ pi+ pi+ pi- pi0 pi0 (n S+ K+ pi+ pi- pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.011, 0.015, 0.018, 0.021},

  //  p S- K0 pi+ pi+ pi0 pi0 pi0 pi0 (n S+ K+ pi- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  n S- K+ pi+ pi+ pi+ pi+ pi- pi- (p S+ K0 pi+ pi+ pi- pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  n S- K+ pi+ pi+ pi+ pi- pi0 pi0 (p S+ K0 pi+ pi- pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.011, 0.014, 0.018, 0.021},

  //  n S- K+ pi+ pi+ pi0 pi0 pi0 pi0 (p S+ K0 pi- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  n S- K0 pi+ pi+ pi+ pi+ pi- pi0 (p S+ K+ pi+ pi- pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011},

  //  n S- K0 pi+ pi+ pi+ pi0 pi0 pi0 (p S+ K+ pi- pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p p K+ K0b pi+ pi+ pi- pi- pi- (n n K- K0 pi+ pi+ pi+ pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004},

  //  p p K+ K0b pi+ pi- pi- pi0 pi0 (n n K- K0 pi+ pi+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.005, 0.009, 0.013, 0.018, 0.023},

  //  p p K+ K0b pi- pi0 pi0 pi0 pi0 (n n K- K0 pi+ pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004},

  //  p p K+ K- pi+ pi+ pi- pi- pi0 (n n K0 K0b pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.005, 0.009, 0.013, 0.018, 0.023},

  //  p p K+ K- pi+ pi- pi0 pi0 pi0 (n n K0 K0b pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.012, 0.016},

  //  p p K+ K- pi0 pi0 pi0 pi0 pi0 (n n K0 K0b pi0 pi0 pi0 pi0 pi0) negligible

  //  p p K0 K0b pi+ pi+ pi- pi- pi0 (n n K+ K- pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.005, 0.009, 0.013, 0.018, 0.023},

  //  p p K0 K0b pi+ pi- pi0 pi0 pi0 (n n K+ K- pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.012, 0.016},

  //  p p K0 K0b pi0 pi0 pi0 pi0 pi0 (n n K+ K- pi0 pi0 pi0 pi0 pi0) negligible

  //  p p K0 K- pi+ pi+ pi+ pi- pi- (n n K+ K0b pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004},

  //  p p K0 K- pi+ pi+ pi- pi0 pi0 (n n K+ K0b pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.005, 0.009, 0.013, 0.018, 0.023},

  //  p p K0 K- pi+ pi0 pi0 pi0 pi0 (n n K+ K0b pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004},

  //  p n K+ K0b pi+ pi+ pi- pi- pi0 (p n K- K0 pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.004, 0.012, 0.018, 0.026, 0.036, 0.045},

  //  p n K+ K0b pi+ pi- pi0 pi0 pi0 (p n K- K0 pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.012, 0.018, 0.024, 0.031},

  //  p n K+ K0b pi0 pi0 pi0 pi0 pi0 (p n K- K0 pi0 pi0 pi0 pi0 pi0) negligible

  //  p n K+ K- pi+ pi+ pi+ pi- pi- (p n K0 K0b pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.012, 0.016},

  //  p n K+ K- pi+ pi+ pi- pi0 pi0 (p n K0 K0b pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.004, 0.012, 0.018, 0.026, 0.036, 0.045},

  //  p n K+ K- pi+ pi0 pi0 pi0 pi0 (p n K0 K0b pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p n K0 K0b pi+ pi+ pi+ pi- pi- (p n K+ K- pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.004, 0.006, 0.009, 0.012, 0.016},

  //  p n K0 K0b pi+ pi+ pi- pi0 pi0 (p n K+ K- pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.004, 0.011, 0.018, 0.026, 0.036, 0.045},

  //  p n K0 K0b pi+ pi0 pi0 pi0 pi0 (p n K+ K- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  p n K0 K- pi+ pi+ pi+ pi- pi0 (p n K+ K0b pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.012, 0.018, 0.024, 0.031},

  //  p n K0 K- pi+ pi+ pi0 pi0 pi0 (p n K+ K0b pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.004, 0.006, 0.009, 0.012, 0.016},

  //  n n K+ K0b pi+ pi+ pi+ pi- pi- (p p K- K0 pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  n n K+ K0b pi+ pi+ pi- pi0 pi0 (p p K- K0 pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.005, 0.009, 0.013, 0.018, 0.023},

  //  n n K+ K0b pi+ pi0 pi0 pi0 pi0 (p p K- K0 pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004},

  //  n n K+ K- pi+ pi+ pi+ pi- pi0 (p p K0 K0b pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.012, 0.016},

  //  n n K+ K- pi+ pi+ pi0 pi0 pi0 (p p K0 K0b pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  n n K0 K0b pi+ pi+ pi+ pi- pi0 (p p K+ K- pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.012, 0.016},

  //  n n K0 K0b pi+ pi+ pi0 pi0 pi0 (p p K+ K- pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  n n K0 K- pi+ pi+ pi+ pi+ pi- (p p K+ K0b pi+ pi- pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004},

  //  n n K0 K- pi+ pi+ pi+ pi0 pi0 (p p K+ K0b pi- pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006}};

}

// Initialize p-p cross-section table

const G4CascadePPChannelData::data_t
G4CascadePPChannelData::data(pp2bfs, pp3bfs, pp4bfs, pp5bfs, pp6bfs, pp7bfs,
			     pp8bfs, pp9bfs, ppCrossSections, ppTotXSec,
			     pro*pro, "ProtonProton");


// Overload base class interpolator to use function for 0-10 MeV total, elastic

G4double 
G4CascadePPChannel::findCrossSection(G4double ke,
                                     const G4double (&xsec)[30]) const {
  if (ke < 0.01 &&
       (std::equal(std::cbegin(xsec), std::cend(xsec), std::cbegin(ppTotXSec))
     || std::equal(std::cbegin(xsec), std::cend(xsec), std::cbegin(ppCrossSections[0]))))
  {
    // Stepanov's function for ke < 10 MeV, up to zero-energy value
    const G4double kemin = 4.0/ppTotXSec[0];
    return (ke>0.001 ? (9.0692 - 0.0050574/ke)/ke + 6.9466 :
	    ke>kemin ? 4.0/ke : ppTotXSec[0]);
  }
  return G4PionNucSampler::findCrossSection(ke, xsec);	// Call through to base
}
