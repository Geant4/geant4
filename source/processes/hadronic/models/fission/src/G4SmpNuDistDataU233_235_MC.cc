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
//
// This software was developed by Lawrence Livermore National Laboratory.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// 3. The name of the author may not be used to endorse or promote products
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Copyright (c) 2006 The Regents of the University of California.
// All rights reserved.
// UCRL-CODE-224807
//
//
// $Id: G4SmpNuDistDataU233_235_MC.cc 68799 2013-04-05 13:29:46Z gcosmo $
//

#include <cmath>
#include "G4fissionEvent.hh"

G4int G4fissionEvent::G4SmpNuDistDataU233_235_MC(G4double nubar) {

/*
  Description
    Sample Number of Neutrons from fission in U-233 and U-235 using 
    Zucker and Holden's tabulated data for U-235
    The 11 P(nu) distributions are given as a function of nubar, 
    the average number of neutrons from induced fission for the 
    11 different energies (0 to 10 MeV), based on the U-235 data 
    from Zucker and Holden.
*/

/*
  Input
    nubar    - average number of neutrons per fission
  Output
    G4SmpNuDistDataU233_235_MC  - sampled multiplicity
    
*/

  static G4double U235nu [11] [8] = {
     {.0317223, .1717071, .3361991, .3039695, .1269459, .0266793, .0026322, .0001449},
     {.0237898, .1555525, .3216515, .3150433, .1444732, .0356013, .0034339, .0004546},
     {.0183989, .1384891, .3062123, .3217566, .1628673, .0455972, .0055694, .0011093},
     {.0141460, .1194839, .2883075, .3266568, .1836014, .0569113, .0089426, .0019504},
     {.0115208, .1032624, .2716849, .3283426, .2021206, .0674456, .0128924, .0027307},
     {.0078498, .0802010, .2456595, .3308175, .2291646, .0836912, .0187016, .0039148},
     {.0046272, .0563321, .2132296, .3290407, .2599806, .1045974, .0265604, .0056322},
     {.0024659, .0360957, .1788634, .3210507, .2892537, .1282576, .0360887, .0079244},
     {.0012702, .0216090, .1472227, .3083032, .3123950, .1522540, .0462449, .0107009},
     {.0007288, .0134879, .1231200, .2949390, .3258251, .1731879, .0551737, .0135376},
     {.0004373, .0080115, .1002329, .2779283, .3342611, .1966100, .0650090, .0175099}
    };
  static G4double U235nubar [11] = {
      2.4140000, 
      2.5236700, 
      2.6368200, 
      2.7623400, 
      2.8738400, 
      3.0386999, 
      3.2316099, 
      3.4272800, 
      3.6041900, 
      3.7395900, 
      3.8749800
    };
  G4double fraction, r, cum;
  G4int engind, nu;

/* 
  Check if nubar is within the range of experimental values
*/
  if(nubar >= U235nubar[0] && nubar <= U235nubar[10]) {
/*
     Use Zucker and Holden Data
*/
     engind = 1;
     while (nubar > U235nubar[engind]){ engind++;}
     // Loop checking, 11.03.2015, T. Koi
     fraction = (nubar-U235nubar[engind-1])/(U235nubar[engind]-U235nubar[engind-1]);
     if(fisslibrng() > fraction) engind--;

     r = fisslibrng();
     nu = 0;
     cum = U235nu[engind][0];
     while (r > cum && nu < 7){ 
     // Loop checking, 11.03.2015, T. Koi
       nu++;
       cum += U235nu[engind][nu];
     }
     return nu;
  } else {
/*
     Use Terrell's formula
*/
     return (G4int) G4SmpTerrell(nubar);
  }
}
