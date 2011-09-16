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
// $Id: G4CascadeXiMinusPChannel.cc,v 1.6 2010-12-15 07:40:54 gunter Exp $
//
// 20100804  M. Kelsey -- Add name string to ctor
// 20110719  M. Kelsey -- Add initial state code to ctor
// 20110725  M. Kelsey -- Instantiate cross-section object for self-registration
// 20110916  M. Kelsey -- Drop self-registration due to platform inconsistencies

#include "G4CascadeXiMinusPChannel.hh"
#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

namespace {
  // Outgoing particle types of a given multiplicity

  static const G4int xmp2bfs[6][2] =
    {{1,31}, {2,29}, {21,21}, {21,25}, {25,25}, {23,27}};

  static const G4int xmp3bfs[24][3] =
    {{1,13,21},  {1,13,25}, {1,17,27}, {1,5,29},   {1,7,31},   {2,17,21},
     {2,17,25},  {2,13,23}, {2,7,29},  {2,3,31},   {7,21,21},  {7,21,25},
     {5,21,23},  {3,21,27}, {15,21,29},{11,21,31}, {5,23,25},  {7,23,27}, 
     {15,23,31}, {7,25,25}, {3,25,27}, {15,25,29}, {11,25,31}, {11,27,29}};

  static const G4int xmp4bfs[4][4] =
    {{1,7,13,21}, {2,7,17,21}, {1,3,5,31}, {2,3,5,29}};

  static const G4int xmp5bfs[4][5] =
    {{1,3,5,13,21}, {2,3,5,17,21}, {1,3,5,7,31}, {2,3,5,7,29}};

  static const G4int xmp6bfs[4][6] =
    {{1,3,5,7,13,21}, {2,3,5,7,17,21}, {1,3,3,5,5,31}, {2,3,3,5,5,29}};

  static const G4int xmp7bfs[4][7] =
    {{1,3,3,5,5,13,21}, {2,3,3,5,5,17,21}, {1,3,3,5,5,7,31}, {2,3,3,5,5,7,29}}; 

  // Cross sections for X- p -> 2-7 body final states
  // 
  // first index:    0-5: channels for mult = 2
  //                6-29: channels for mult = 3 
  //               30-33: channels for mult = 4
  //               34-37: channels for mult = 5
  //               38-41: channels for mult = 6
  //               42-45: channels for mult = 7
  //
  // second index: kinetic energy
  // 
  static const G4double xmpCrossSections[46][31] = {
    //
    // multiplicity 2 (6 channels)
    //
    // X- p
    {22.00,20.00,18.00,16.00,15.00,14.00,13.00,12.00,11.00,10.00,
     9.00, 6.00, 5.50, 5.00, 4.50, 4.00, 3.70, 3.30, 3.00, 2.70,
     2.50, 2.20, 2.00, 1.80, 1.60, 1.40, 1.20, 1.10, 1.00, 0.90, 0.70},
 
    // X0 n
    {11.00,10.50,10.00, 9.50, 9.00, 8.50, 8.30, 8.00, 7.70, 7.50,
     7.20, 4.00, 3.50, 3.00, 2.50, 2.30, 2.00, 1.70, 1.50, 1.35,
     1.25, 1.10, 1.00, 0.90, 0.80, 0.70, 0.60, 0.55, 0.50, 0.45, 0.35},
 
    // L L
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.50, 2.00, 2.30, 2.50, 2.80, 2.50, 2.20, 2.00, 1.70, 1.50,
      1.40, 1.30, 1.20, 1.10, 1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.30},
        
    // L S0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.50, 2.00, 2.30, 2.50, 2.80, 2.50, 2.20, 2.00, 1.70, 1.50,
      1.40, 1.30, 1.20, 1.10, 1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.30},
 
    // S0 S0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.50, 2.00, 2.30, 2.50, 2.80, 2.50, 2.20, 2.00, 1.70, 1.50,
      1.40, 1.30, 1.20, 1.10, 1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.30},
 
    // S+ S-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.50, 2.00, 2.30, 2.50, 2.80, 2.50, 2.20, 2.00, 1.70, 1.50,
      1.40, 1.30, 1.20, 1.10, 1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.30},
 
    //
    //  multiplicity 3 (24 channels)
    //
    // p L K-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.10, 0.50, 0.60, 0.67, 0.73, 0.90, 0.90, 0.80, 0.70, 0.60,
      0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.12, 0.10, 0.05},
 
    // p S0 K-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.10, 0.50, 0.60, 0.67, 0.73, 0.90, 0.90, 0.80, 0.70, 0.60,
      0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.12, 0.10, 0.05},
 
    // p S- K0bar
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.10, 0.50, 0.60, 0.67, 0.73, 0.90, 0.90, 0.80, 0.70, 0.60,
      0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.12, 0.10, 0.05},
 
    // p X0 pi-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.05, 0.25, 0.27, 0.33, 0.37, 0.45, 0.45, 0.40, 0.35, 0.30,
      0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03},
 
    // p X- pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.05, 0.25, 0.27, 0.33, 0.37, 0.45, 0.45, 0.40, 0.35, 0.30,
      0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03},

    // n L K0bar
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.10, 0.50, 0.60, 0.67, 0.73, 0.90, 0.90, 0.80, 0.70, 0.60,
      0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.12, 0.10, 0.05},
 
    // n S0 K0bar
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.10, 0.50, 0.60, 0.67, 0.73, 0.90, 0.90, 0.80, 0.70, 0.60,
      0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.12, 0.10, 0.05},

    // n S+ K-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.10, 0.50, 0.60, 0.67, 0.73, 0.90, 0.90, 0.80, 0.70, 0.60,
      0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.12, 0.10, 0.05},
 
    // n X0 pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.05, 0.25, 0.27, 0.33, 0.37, 0.45, 0.45, 0.40, 0.35, 0.30,
      0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03},
 
    // n X- pi+
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.05, 0.25, 0.27, 0.33, 0.37, 0.45, 0.45, 0.40, 0.35, 0.30,
      0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03},
 
    // L L pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.05, 0.25, 0.27, 0.33, 0.37, 0.45, 0.45, 0.40, 0.35, 0.30,
      0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03},
 
    // L S0 pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.05, 0.25, 0.27, 0.33, 0.37, 0.45, 0.45, 0.40, 0.35, 0.30,
      0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03},
 
    // L S+ pi-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.05, 0.25, 0.27, 0.33, 0.37, 0.45, 0.45, 0.40, 0.35, 0.30,
      0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03},
 
    // L S- pi+
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.05, 0.25, 0.27, 0.33, 0.37, 0.45, 0.45, 0.40, 0.35, 0.30,
      0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03},
 
    // L X0 K0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01, 0.02, 0.02, 0.03,
      0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},
 
    // L X- K+
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01, 0.02, 0.02, 0.03,
      0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

    // S+ S0 pi-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.05, 0.25, 0.27, 0.33, 0.37, 0.45, 0.45, 0.40, 0.35, 0.30,
      0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03},

    // S+ S- pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.05, 0.25, 0.27, 0.33, 0.37, 0.45, 0.45, 0.40, 0.35, 0.30,
      0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03},
 
    // S+ X- K0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01, 0.02, 0.02, 0.03,
      0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},
 
    // S0 S0 pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.05, 0.25, 0.27, 0.33, 0.37, 0.45, 0.45, 0.40, 0.35, 0.30,
      0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03},
 
    // S0 S- pi+
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.05, 0.25, 0.27, 0.33, 0.37, 0.45, 0.45, 0.40, 0.35, 0.30,
      0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03},
 
    // S0 X0 K0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01, 0.02, 0.02, 0.03,
      0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},
 
    // S0 X- K+
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01, 0.02, 0.02, 0.03,
      0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},
 
    // S- X0 K+
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.01, 0.02, 0.02, 0.03,
      0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02},

    //
    //  multiplicity 4 (4 channels)
    //
    // p L K- pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.07, 0.15, 0.30, 0.35, 0.38, 0.40, 0.43, 0.45,
      0.47, 0.50, 0.53, 0.50, 0.47, 0.43, 0.37, 0.30, 0.25, 0.20, 0.05},
 
    // n L K0bar pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.07, 0.15, 0.30, 0.35, 0.38, 0.40, 0.43, 0.45,
      0.47, 0.50, 0.53, 0.50, 0.47, 0.43, 0.37, 0.30, 0.25, 0.20, 0.05},
 
    // p X- pi+ pi-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.03, 0.07, 0.13, 0.15, 0.17, 0.20, 0.22, 0.23,
      0.24, 0.25, 0.26, 0.25, 0.24, 0.23, 0.22, 0.20, 0.17, 0.15, 0.05},
 
    // n X0 pi+ pi-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.03, 0.07, 0.13, 0.15, 0.17, 0.20, 0.22, 0.23,
      0.24, 0.25, 0.26, 0.25, 0.24, 0.23, 0.22, 0.20, 0.17, 0.15, 0.05},

    //
    //  multiplicity 5 (4 channels)
    // 
    // p L K- pi+ pi-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.05, 0.10, 0.15, 0.17,
      0.19, 0.21, 0.23, 0.25, 0.27, 0.27, 0.25, 0.23, 0.21, 0.19, 0.15},
 
    // n L K0bar pi+ pi-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.05, 0.10, 0.15, 0.17,
      0.19, 0.21, 0.23, 0.25, 0.27, 0.27, 0.25, 0.23, 0.21, 0.19, 0.15},
 
    // p X- pi+ pi- pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.05, 0.07, 0.08,
      0.09, 0.10, 0.11, 0.12, 0.13, 0.13, 0.12, 0.11, 0.10, 0.09, 0.07},
 
    // n X0 pi+ pi- pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.05, 0.07, 0.08,
      0.09, 0.10, 0.11, 0.12, 0.13, 0.13, 0.12, 0.11, 0.10, 0.09, 0.07},

    //
    //  multiplicity 6 (4 channels)
    // 
    // p L K- pi+ pi- pi0     
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 
      0.05, 0.10, 0.15, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.22, 0.21},
 
    // n L K0bar pi+ pi- pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 
      0.05, 0.10, 0.15, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.22, 0.21},
 
    // p X- 2pi+ 2pi-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 
      0.02, 0.05, 0.07, 0.08, 0.09, 0.09, 0.10, 0.10, 0.11, 0.11, 0.10},
 
    // n X0 2pi+ 2pi-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 
      0.02, 0.05, 0.07, 0.08, 0.09, 0.09, 0.10, 0.10, 0.11, 0.11, 0.10},
 
    //
    //  multiplicity 7 (4 channels)
    // 
    // p L K- 2pi+ 2pi-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.00, 0.0,  0.0,  0.0,
      0.01, 0.01, 0.01, 0.02, 0.02, 0.03, 0.03, 0.04, 0.04, 0.05, 0.05},
 
    // n L K0bar 2pi+ 2pi-
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.00, 0.0,  0.0,  0.0,
      0.01, 0.01, 0.01, 0.02, 0.02, 0.03, 0.03, 0.04, 0.04, 0.05, 0.05},
 
    // p X- 2pi+ 2pi- pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.00, 0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02},
 
    // n X0 2pi+ 2pi- pi0
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.00, 0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02}};

}

G4CascadeXiMinusPChannelData::data_t
G4CascadeXiMinusPChannelData::data(xmp2bfs, xmp3bfs, xmp4bfs,
				   xmp5bfs, xmp6bfs, xmp7bfs,
				   xmpCrossSections, xim*pro, "XiMinusP");
