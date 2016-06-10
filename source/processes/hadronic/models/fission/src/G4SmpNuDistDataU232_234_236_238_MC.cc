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
// $Id: G4SmpNuDistDataU232_234_236_238_MC.cc 68799 2013-04-05 13:29:46Z gcosmo $
//

#include <cmath>
#include "G4fissionEvent.hh"

G4int G4fissionEvent::G4SmpNuDistDataU232_234_236_238_MC(G4double nubar) {

/*
  Description
    Sample Number of Neutrons from fission in U-232, U-234, U-236, 
    and U-238 using Zucker and Holden's tabulated data for U-238
    The 11 P(nu) distributions are given as a function of nubar, 
    the average number of neutrons from induced fission for the 
    11 different energies (0 to 10 MeV), based on the U-238 data 
    from Zucker and Holden.
*/

/*
  Input
    nubar    - average number of neutrons per fission
  Output
    G4SmpNuDistDataU232_234_236_238_MC  - sampled multiplicity
    
*/

  static G4double U238nu [11] [9] = {
     {.0396484, .2529541, .2939544, .2644470, .1111758, .0312261, .0059347, .0005436, .0001158},
     {.0299076, .2043215, .2995886, .2914889, .1301480, .0363119, .0073638, .0006947, .0001751},
     {.0226651, .1624020, .2957263, .3119098, .1528786, .0434233, .0097473, .0009318, .0003159},
     {.0170253, .1272992, .2840540, .3260192, .1779579, .0526575, .0130997, .0013467, .0005405},
     {.0124932, .0984797, .2661875, .3344938, .2040116, .0640468, .0173837, .0020308, .0008730},
     {.0088167, .0751744, .2436570, .3379711, .2297901, .0775971, .0225619, .0030689, .0013626},
     {.0058736, .0565985, .2179252, .3368863, .2541575, .0933127, .0286200, .0045431, .0031316},
     {.0035997, .0420460, .1904095, .3314575, .2760413, .1112075, .0355683, .0065387, .0031316},
     {.0019495, .0309087, .1625055, .3217392, .2943792, .1313074, .0434347, .0091474, .0046284},
     {.0008767, .0226587, .1356058, .3076919, .3080816, .1536446, .0522549, .0124682, .0067176},
     {.0003271, .0168184, .1111114, .2892434, .3160166, .1782484, .0620617, .0166066, .0095665}
    };
  static G4double U238nubar [11] = {
      2.2753781,
      2.4305631,
      2.5857481,
      2.7409331,
      2.8961181,
      3.0513031,
      3.2064881,
      3.3616731,
      3.5168581,
      3.6720432,
      3.8272281
    };
  G4double fraction, r, cum;
  G4int engind, nu;

/* 
  Check if nubar is within the range of experimental values
*/
  if(nubar >= U238nubar[0] && nubar <= U238nubar[10]) {
/*
     Use Zucker and Holden Data
*/
     engind = 1;
     while (nubar > U238nubar[engind]){ engind++;}
     // Loop checking, 11.03.2015, T. Koi
     fraction = (nubar-U238nubar[engind-1])/(U238nubar[engind]-U238nubar[engind-1]);
     if(fisslibrng() > fraction) engind--;

     r = fisslibrng();
     nu = 0;
     cum = U238nu[engind][0];
     while (r > cum && nu < 8){ 
     // Loop checking, 11.03.2015, T. Koi
       nu++;
       cum += U238nu[engind][nu];
     }
     return nu;
  } else {
/*
     Use Terrell's formula
*/
     return (G4int) G4SmpTerrell(nubar);
  }
}
