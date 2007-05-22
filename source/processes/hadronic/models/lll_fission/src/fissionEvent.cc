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
// $Id: fissionEvent.cc,v 1.1 2007-05-22 00:53:11 dennis Exp $
//

#include "fissionEvent.hh"
#include <stdio.h>
#include <stdlib.h>

int fissionEvent::delayoption=0;
int fissionEvent::correlationoption=0;
int fissionEvent::nudistoption=3;
int fissionEvent::Cf252ndistoption=0;
int fissionEvent::Cf252nengoption=0;

fissionEvent::fissionEvent(int isotope, double time, double nubar, double eng) {
   /*
    * Constructs a fission event with neutronNu neutrons and photonNu
    * photons.
    */
   int i;

   neutronNu = 0;
   photonNu = 0;
   if (nubar == -1.) {
      /* spontaneous fission */
      neutronNu = SmpSpNuDistData(isotope, Cf252ndistoption);
      photonNu = SmpSpNugDistData(isotope);
   } else {
      /* induced fission */
      if (nudistoption == 0 || nudistoption == 1) {
         switch (isotope) {
            case 92235:
               neutronNu = SmpNuDistDataU235(eng,nudistoption);
               break;
            case 92238:
               neutronNu = SmpNuDistDataU238(eng);
               break;
            case 94239:
               neutronNu = SmpNuDistDataPu239(eng);
               break;
            default:
               neutronNu = (int) SmpTerrell(nubar);
               break;
         } 
      } else if (nudistoption == 2) {
         switch (isotope) {
            case 92232:
            case 92234:
            case 92236:
            case 92238:
               neutronNu = SmpNuDistDataU232_234_236_238(nubar);
               break;
            case 92233:
            case 92235:
               neutronNu = (int) SmpNuDistDataU233_235(nubar);
               break;
            case 94239:
            case 94241:
               neutronNu = SmpNuDistDataPu239_241(nubar);
               break;
            default:
               neutronNu = (int) SmpTerrell(nubar);
               break;
         }
      } else if (nudistoption == 3) {
         switch (isotope) {
            case 92232:
            case 92234:
            case 92236:
            case 92238:
               neutronNu = SmpNuDistDataU232_234_236_238_MC(nubar);
               break;
            case 92233:
            case 92235:
               neutronNu = (int) SmpNuDistDataU233_235_MC(nubar);
               break;
            case 94239:
            case 94241:
               neutronNu = SmpNuDistDataPu239_241_MC(nubar);
               break;
            default:
               neutronNu = (int) SmpTerrell(nubar);
               break;
         } 
      }
      photonNu = SmpNugDist(isotope, nubar);
   }
   if (neutronNu > 0) {
      neutronEnergies = new double[ neutronNu ];
      neutronVelocities = new double[ neutronNu ];
      neutronDircosu = new double[ neutronNu ];
      neutronDircosv = new double[ neutronNu ];
      neutronDircosw = new double[ neutronNu ];
      neutronAges = new double[neutronNu];
      for (i=0; i<neutronNu; i++) {
         if (isotope == 98252) neutronEnergies[i] = SmpNEngCf252(Cf252nengoption);
         else neutronEnergies[i] = SmpWatt(eng, isotope);
         neutronVelocities[i] = SmpNVel(
                 neutronEnergies[i],
                 &(neutronDircosu[i]),
                 &(neutronDircosv[i]),
                 &(neutronDircosw[i])
                );
         neutronAges[i] = time;
      }
   }
   if (photonNu > 0) {
      photonEnergies = new double[photonNu];
      photonVelocities = new double[photonNu];
      photonDircosu = new double[photonNu];
      photonDircosv = new double[photonNu];
      photonDircosw = new double[photonNu];
      photonAges = new double[photonNu];
      for (i=0; i<photonNu; i++) {
         photonEnergies[i] = SmpGEng();
         photonVelocities[i] = SmpPVel(
                 photonEnergies[i],
                 &(photonDircosu[i]),
                 &(photonDircosv[i]),
                 &(photonDircosw[i])
                );
         photonAges[i] = time;
      }
   }
};

fissionEvent::~fissionEvent() {
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
};
