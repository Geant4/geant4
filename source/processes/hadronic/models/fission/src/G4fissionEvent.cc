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
// $Id: G4fissionEvent.cc 68799 2013-04-05 13:29:46Z gcosmo $
//

#include "G4fissionEvent.hh"

G4int G4fissionEvent::delayoption=0;
G4int G4fissionEvent::correlationoption=0;
G4int G4fissionEvent::nudistoption=3;
G4int G4fissionEvent::Cf252ndistoption=0;
G4int G4fissionEvent::Cf252nengoption=0;

G4fissionEvent::G4fissionEvent(G4int isotope, G4double time,
                               G4double nubar, G4double eng)
 :neutronNu(0), neutronEnergies(0), neutronVelocities(0), neutronDircosu(0),
  neutronDircosv(), neutronDircosw(), neutronAges(0),
  photonNu(0), photonEnergies(0), photonVelocities(0), photonDircosu(0),
  photonDircosv(0), photonDircosw(0), photonAges(0)
{
   /*
    * Constructs a fission event with neutronNu neutrons and photonNu
    * photons.
    */
   G4int i;

   if (nubar == -1.) {
      /* spontaneous fission */
      neutronNu = G4SmpSpNuDistData(isotope, Cf252ndistoption);
      photonNu = G4SmpSpNugDistData(isotope);
   } else {
      /* induced fission */
      if (nudistoption == 0 || nudistoption == 1) {
         switch (isotope) {
            case 92235:
               neutronNu = G4SmpNuDistDataU235(eng,nudistoption);
               break;
            case 92238:
               neutronNu = G4SmpNuDistDataU238(eng);
               break;
            case 94239:
               neutronNu = G4SmpNuDistDataPu239(eng);
               break;
            default:
               neutronNu = (G4int) G4SmpTerrell(nubar);
               break;
         } 
      } else if (nudistoption == 2) {
         switch (isotope) {
            case 92232:
            case 92234:
            case 92236:
            case 92238:
               neutronNu = G4SmpNuDistDataU232_234_236_238(nubar);
               break;
            case 92233:
            case 92235:
               neutronNu = (G4int) G4SmpNuDistDataU233_235(nubar);
               break;
            case 94239:
            case 94241:
               neutronNu = G4SmpNuDistDataPu239_241(nubar);
               break;
            default:
               neutronNu = (G4int) G4SmpTerrell(nubar);
               break;
         }
      } else if (nudistoption == 3) {
         switch (isotope) {
            case 92232:
            case 92234:
            case 92236:
            case 92238:
               neutronNu = G4SmpNuDistDataU232_234_236_238_MC(nubar);
               break;
            case 92233:
            case 92235:
               neutronNu = (G4int) G4SmpNuDistDataU233_235_MC(nubar);
               break;
            case 94239:
            case 94241:
               neutronNu = G4SmpNuDistDataPu239_241_MC(nubar);
               break;
            default:
               neutronNu = (G4int) G4SmpTerrell(nubar);
               break;
         } 
      }
      photonNu = G4SmpNugDist(isotope, nubar);
   }
   if (neutronNu > 0) {
      neutronEnergies = new G4double[ neutronNu ];
      neutronVelocities = new G4double[ neutronNu ];
      neutronDircosu = new G4double[ neutronNu ];
      neutronDircosv = new G4double[ neutronNu ];
      neutronDircosw = new G4double[ neutronNu ];
      neutronAges = new G4double[neutronNu];
      for (i=0; i<neutronNu; i++) {
         if (isotope == 98252) neutronEnergies[i] = G4SmpNEngCf252(Cf252nengoption);
         else neutronEnergies[i] = G4SmpWatt(eng, isotope);
         neutronVelocities[i] = G4SmpNVel(
                 neutronEnergies[i],
                 &(neutronDircosu[i]),
                 &(neutronDircosv[i]),
                 &(neutronDircosw[i])
                );
         neutronAges[i] = time;
      }
   }
   if (photonNu > 0) {
      photonEnergies = new G4double[photonNu];
      photonVelocities = new G4double[photonNu];
      photonDircosu = new G4double[photonNu];
      photonDircosv = new G4double[photonNu];
      photonDircosw = new G4double[photonNu];
      photonAges = new G4double[photonNu];
      for (i=0; i<photonNu; i++) {
         photonEnergies[i] = G4SmpGEng();
         photonVelocities[i] = G4SmpPVel(
                 photonEnergies[i],
                 &(photonDircosu[i]),
                 &(photonDircosv[i]),
                 &(photonDircosw[i])
                );
         photonAges[i] = time;
      }
   }
}

G4fissionEvent::~G4fissionEvent() {
   if (neutronNu > 0) {
      delete [] neutronEnergies;
      delete [] neutronVelocities;
      delete [] neutronDircosu;
      delete [] neutronDircosv;
      delete [] neutronDircosw;
      delete [] neutronAges;
   }

   if (photonNu > 0) {
      delete [] photonEnergies;
      delete [] photonVelocities;
      delete [] photonDircosu;
      delete [] photonDircosv;
      delete [] photonDircosw;
      delete [] photonAges;
   }
}
