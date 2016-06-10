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
// $Id: G4SmpTerrell.cc 67966 2013-03-13 09:38:38Z gcosmo $
//

#include <cmath>
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4fissionEvent.hh"

#define TWOPI 6.283185307
#define SQRT2 1.414213562
#define BSHIFT -0.43287
#define WIDTH 1.079

G4double G4fissionEvent::G4SmpTerrell(G4double nubar) {
/*
  Description
    Sample Fission Number from Terrell's modified Gaussian distribution

    method uses Red Cullen's algoritm UCRL-TR-222526
*/

/*
  Input
    nubar    - average number of neutrons per fission
  Output
    G4SmpTerrell  - sampled multiplicity
    
*/

  G4double width;
  G4double temp1, temp2, expo, cshift;
  G4double rw, theta, sampleg;


  if (nubar < WIDTH) {
    std::ostringstream o;
    o << nubar;
    std::string errMsg = "fission nubar out of range, nubar=" + o.str();
    G4fissionerr(6, "SmpTerrell", errMsg);
  }

  width = SQRT2 * WIDTH;
  temp1 = nubar + 0.5;
  temp2 = temp1/width;
  temp2 *= temp2;
  expo = G4Exp(-temp2);
  cshift = temp1 + BSHIFT * WIDTH * expo/(1. - expo);

  G4int icounter = 0;
  G4int icounter_max = 1024;
  do {
    rw = std::sqrt(-G4Log(fisslibrng()));
    theta = TWOPI * fisslibrng();
    sampleg = width * rw * std::cos(theta) + cshift;
    icounter++;
    if ( icounter > icounter_max ) { 
      G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
      break;
    }
  } while (sampleg < 0.0);
  // Loop checking, 11.03.2015, T. Koi

  return std::floor(sampleg);
}
