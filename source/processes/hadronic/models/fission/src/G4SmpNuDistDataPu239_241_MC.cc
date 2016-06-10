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
// $Id: G4SmpNuDistDataPu239_241_MC.cc 68799 2013-04-05 13:29:46Z gcosmo $
//

#include <cmath>
#include "G4fissionEvent.hh"

G4int G4fissionEvent::G4SmpNuDistDataPu239_241_MC(G4double nubar) {

/*
  Description
    Sample Number of Neutrons from fission in Pu-239 and Pu-241 using 
    Zucker and Holden's tabulated data for Pu-239
    The 11 P(nu) distributions are given as a function of nubar, 
    the average number of neutrons from induced fission for the 
    11 different energies (0 to 10 MeV), based on the Pu-239 data 
    from Zucker and Holden.
*/

/*
  Input
    nubar    - average number of neutrons per fission
  Output
    G4SmpNuDistDataPu239_241_MC  - sampled multiplicity
    
*/

  static G4double Pu239nu [11] [9] = {
     {.0108826, .0994916, .2748898, .3269196, .2046061, .0726834, .0097282, .0006301, .0001685},
     {.0084842, .0790030, .2536175, .3289870, .2328111, .0800161, .0155581, .0011760, .0003469},
     {.0062555, .0611921, .2265608, .3260637, .2588354, .0956070, .0224705, .0025946, .0005205},
     {.0045860, .0477879, .1983002, .3184667, .2792811, .1158950, .0301128, .0048471, .0007233},
     {.0032908, .0374390, .1704196, .3071862, .2948565, .1392594, .0386738, .0078701, .0010046},
     {.0022750, .0291416, .1437645, .2928006, .3063902, .1641647, .0484343, .0116151, .0014149},
     {.0014893, .0222369, .1190439, .2756297, .3144908, .1892897, .0597353, .0160828, .0029917},
     {.0009061, .0163528, .0968110, .2558524, .3194566, .2134888, .0729739, .0213339, .0020017},
     {.0004647, .0113283, .0775201, .2335926, .3213289, .2356614, .0886183, .0274895, .0039531},
     {.0002800, .0071460, .0615577, .2089810, .3200121, .2545846, .1072344, .0347255, .0054786},
     {.0002064, .0038856, .0492548, .1822078, .3154159, .2687282, .1295143, .0432654, .0075217}
    };
  static G4double Pu239nubar [11] = {
      2.8760000,
      3.0088800,
      3.1628300,
      3.3167800,
      3.4707300,
      3.6246800,
      3.7786300,
      3.9325800,
      4.0865300,
      4.2404900,
      4.3944400
    };
  G4double fraction, r, cum;
  G4int engind, nu;

/* 
  Check if nubar is within the range of experimental values
*/
  if(nubar >= Pu239nubar[0] && nubar <= Pu239nubar[10]) {
/*
     Use Zucker and Holden Data
*/
     engind = 1;
     while (nubar > Pu239nubar[engind]){ engind++;}
     // Loop checking, 11.03.2015, T. Koi
     fraction = (nubar-Pu239nubar[engind-1])/(Pu239nubar[engind]-Pu239nubar[engind-1]);
     if(fisslibrng() > fraction) engind--;

     r = fisslibrng();
     nu = 0;
     cum = Pu239nu[engind][0];
     while (r > cum && nu < 8){ 
     // Loop checking, 11.03.2015, T. Koi
       nu++;
       cum += Pu239nu[engind][nu];
     }
     return nu;
  } else {
/*
     Use Terrell's formula
*/
     return (G4int) G4SmpTerrell(nubar);
  }
}
