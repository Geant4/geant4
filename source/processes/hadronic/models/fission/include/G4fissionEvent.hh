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
// $Id: G4fissionEvent.hh 68799 2013-04-05 13:29:46Z gcosmo $
//
#ifndef G4fissionEvent_hh
#define G4fissionEvent_hh

#include "globals.hh"
#include <string>

class G4fissionEvent {
   private:
      G4int neutronNu; // number of neutrons in this fission event
      G4double* neutronEnergies; 
      G4double* neutronVelocities; 
      G4double* neutronDircosu; 
      G4double* neutronDircosv; 
      G4double* neutronDircosw; 
      G4double* neutronAges; 

      G4int photonNu; // number of photons in this fission event
      G4double* photonEnergies; 
      G4double* photonVelocities; 
      G4double* photonDircosu; 
      G4double* photonDircosv; 
      G4double* photonDircosw; 
      G4double* photonAges; 

      // options
      static G4int delayoption;
      static G4int correlationoption;
      static G4int nudistoption;
      static G4int Cf252ndistoption;
      static G4int Cf252nengoption;
      static G4double (*rngdptr)(void);
      static float (*rngfptr)(void);

   public:
      // These are all the methods of this class accessible to the caller of the object 
      G4fissionEvent(G4int isotope, G4double time, G4double nubar, G4double eng);
      ~G4fissionEvent();
      G4int getNeutronNu() {
         return neutronNu;
      }
      G4int getPhotonNu() {
         return photonNu;
      }
      G4double getNeutronEnergy(G4int index) {
         if (index >= 0 && index < neutronNu) return neutronEnergies[index];
         else return -1;
      }
      G4double getNeutronVelocity(G4int index) {
         if (index >= 0 && index < neutronNu) return neutronVelocities[index];
         else return -1;
      }
      G4double getNeutronDircosu(G4int index) {
         if (index >= 0 && index < neutronNu) return neutronDircosu[index];
         else return -1;
      }
      G4double getNeutronDircosv(G4int index) {
         if (index >= 0 && index < neutronNu) return neutronDircosv[index];
         else return -1;
      }
      G4double getNeutronDircosw(G4int index) {
         if (index >= 0 && index < neutronNu) return neutronDircosw[index];
         else return -1;
      }
      G4double getPhotonEnergy(G4int index) {
         if (index >= 0 && index < photonNu) return photonEnergies[index];
         else return -1;
      }
      G4double getPhotonVelocity(G4int index) {
         if (index >= 0 && index < photonNu) return photonVelocities[index];
         else return -1;
      }
      G4double getPhotonDircosu(G4int index) {
         if (index >= 0 && index < photonNu) return photonDircosu[index];
         else return -1;
      }
      G4double getPhotonDircosv(G4int index) {
         if (index >= 0 && index < photonNu) return photonDircosv[index];
         else return -1;
      }
      G4double getPhotonDircosw(G4int index) {
         if (index >= 0 && index < photonNu) return photonDircosw[index];
         else return -1;
      }
      G4double getNeutronAge(G4int index) {
         if (index >= 0 && index < neutronNu) return neutronAges[index];
         else return -1;
      }
      G4double getPhotonAge(G4int index) {
         if (index >= 0 && index < photonNu) return photonAges[index];
         else return -1;
      }
      static void setDelayOption(G4int delay) {
         delayoption = delay;
      };
      static void setCorrelationOption(G4int correlation) {
         correlationoption = correlation;
      };
      static void setNudistOption(G4int nudist) {
         nudistoption = nudist;
      };
      static void setCf252Option(G4int ndist, G4int neng) {
         Cf252ndistoption = ndist;
         Cf252nengoption = neng;
      };
      static void setRNGf(float (*funcptr) (void)) {
         rngfptr = funcptr;
         rngdptr = rngf2d;
      }
      static void setRNGd(G4double (*funcptr) (void)) {
         rngdptr = funcptr;
      }


   private:
      G4int G4SmpNuDistDataU232_234_236_238(G4double nubar);
      G4int G4SmpNuDistDataU232_234_236_238_MC(G4double nubar);
      G4int G4SmpNuDistDataU233_235(G4double nubar);
      G4int G4SmpNuDistDataU233_235_MC(G4double nubar);
      G4int G4SmpNuDistDataU235(G4double erg, G4int option);
      G4int G4SmpNuDistDataPu239(G4double erg);
      G4double G4SmpNVel(G4double eng, G4double* cosdiru, G4double* cosdirv, G4double* cosdirw);
      G4double G4SmpNEngCf252(G4int option);
      void G4SmpIsoDir(G4double* cosdiru, G4double* cosdirv, G4double* cosdirw);
      G4double G4SmpGEng();
      G4int G4SmpNuDistDataPu239_241(G4double nubar);
      G4int G4SmpNuDistDataPu239_241_MC(G4double nubar);
      G4int G4SmpNuDistDataU238(G4double erg);
      G4int G4SmpNugDist(G4int isotope, G4double nubar);
      G4double G4SmpPVel(G4double eng, G4double* cosdiru, G4double* cosdirv, G4double* cosdirw);
      G4int G4SmpSpNuDistData(G4int isotope, G4int Cf252option);
      G4double G4SmpSpNubarData(G4int isotope);
      G4int G4SmpSpNugDistData(G4int isotope);
      G4double G4SmpTerrell(G4double nubar);
      G4double G4SmpWatt(G4double ePart, G4int iso);
      void G4fissionerr(G4int iSever, std::string chSubNam, std::string chMsg);
      static G4double fisslibrng(void);
      static G4double rngf2d(void);
};

#endif
