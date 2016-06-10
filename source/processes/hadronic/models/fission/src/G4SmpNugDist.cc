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
// $Id: G4SmpNugDist.cc 68799 2013-04-05 13:29:46Z gcosmo $
//

#include <cmath>
#include "G4fissionEvent.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

#define nfissg 40
#define alphanegbin 26

G4int G4fissionEvent::G4SmpNugDist(G4int isotope, G4double nubar) {

/*
  Description
    Sample Number of Photons from neutron induced fission in 
    all isotopes using Tim Valentine's model (negative binomial
    distribution, using nubar as a model parameter)
*/

/*
  Input
    iso          - isotope
  Output
    G4SmpNugDist - sampled multiplicity
*/
 
  static G4double logcoeff[nfissg+1] = {
     0.00000000000000e+00,
     3.25809653802149e+00,
     5.86078622346587e+00,
     8.09437844497297e+00,
     1.00753799138395e+01,
     1.18671393830676e+01,
     1.35093671183247e+01,
     1.50291928720691e+01,
     1.64462588918558e+01,
     1.77753948391357e+01,
     1.90281578076311e+01,
     2.02137814732888e+01,
     2.13397927361450e+01,
     2.24124295384099e+01,
     2.34369338549243e+01,
     2.44177631079360e+01,
     2.53587464524005e+01,
     2.62632027266277e+01,
     2.71340310844251e+01,
     2.79737817391769e+01,
     2.87847119553932e+01,
     2.95688309141589e+01,
     3.03279360625106e+01,
     3.10636428574894e+01,
     3.17774093252521e+01,
     3.24705565058120e+01,
     3.31442856005149e+01,
     3.37996924530920e+01,
     3.44377798564689e+01,
     3.50594680730467e+01,
     3.56656038766170e+01,
     3.62569683628670e+01,
     3.68342837279018e+01,
     3.73982191769817e+01,
     3.79493960962713e+01,
     3.84883925970040e+01,
     3.90157475227212e+01,
     3.95319639951220e+01,
     4.00375125617872e+01,
     4.05328339990172e+01,
     4.10183418147990e+01
  };
  G4int i, A, Z;
  G4double cpi[nfissg+1];
  G4double p, q, nubarg;
  G4double r;

/* 
  No data is available for induced fission gamma number
  distributions. Sample the negative binomial cumulative 
  probability distribution.
*/
  A = (G4int) (isotope-1000*((G4int)(isotope/1000)));
  Z = (G4int) ((isotope-A)/1000);
  G4Pow* Pow = G4Pow::GetInstance();
  nubarg = ((2.51-1.13e-5*Pow->powA(G4double(Z),2.)*std::sqrt(G4double(A)))*nubar+4.0)
           /(-1.33+119.6*Pow->A13(G4double(Z))/G4double(A));
  p = 1.*alphanegbin/(alphanegbin+nubarg);
  q = 1.-p;
  cpi[0] = G4Exp(logcoeff[0]+26.*G4Log(p));
  for (i=1; i<=nfissg; i++) cpi[i] = cpi[i-1] + G4Exp(logcoeff[i]+26.*G4Log(p)+i*G4Log(q));
  for (i=0; i<=nfissg; i++) cpi[i] = cpi[i]/cpi[nfissg-1];

  r=fisslibrng();

  for(i=0; i<=nfissg; i++) if (r <= cpi[i]) return i;

  //
  // Fall through
  //

  G4cout << " SmpNugDist: random number " << r << " out of range " << G4endl;
  return -1;

}
