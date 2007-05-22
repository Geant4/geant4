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
// $Id: fissionEvent.hh,v 1.1 2007-05-22 00:49:32 dennis Exp $
//

#include <string>

class fissionEvent {
   private:
      int neutronNu; // number of neutrons in this fission event
      double* neutronEnergies; 
      double* neutronVelocities; 
      double* neutronDircosu; 
      double* neutronDircosv; 
      double* neutronDircosw; 
      double* neutronAges; 

      int photonNu; // number of photons in this fission event
      double* photonEnergies; 
      double* photonVelocities; 
      double* photonDircosu; 
      double* photonDircosv; 
      double* photonDircosw; 
      double* photonAges; 

      // options
      static int delayoption;
      static int correlationoption;
      static int nudistoption;
      static int Cf252ndistoption;
      static int Cf252nengoption;
      static double (*rngdptr)(void);
      static float (*rngfptr)(void);

   public:
      // These are all the methods of this class accessible to the caller of the object 
      fissionEvent(int isotope, double time, double nubar, double eng);
      ~fissionEvent();
      int getNeutronNu() {
         return neutronNu;
      }
      int getPhotonNu() {
         return photonNu;
      }
      double getNeutronEnergy(int index) {
         if (index >= 0 && index < neutronNu) return neutronEnergies[index];
         else return -1;
      }
      double getNeutronVelocity(int index) {
         if (index >= 0 && index < neutronNu) return neutronVelocities[index];
         else return -1;
      }
      double getNeutronDircosu(int index) {
         if (index >= 0 && index < neutronNu) return neutronDircosu[index];
         else return -1;
      }
      double getNeutronDircosv(int index) {
         if (index >= 0 && index < neutronNu) return neutronDircosv[index];
         else return -1;
      }
      double getNeutronDircosw(int index) {
         if (index >= 0 && index < neutronNu) return neutronDircosw[index];
         else return -1;
      }
      double getPhotonEnergy(int index) {
         if (index >= 0 && index < photonNu) return photonEnergies[index];
         else return -1;
      }
      double getPhotonVelocity(int index) {
         if (index >= 0 && index < photonNu) return photonVelocities[index];
         else return -1;
      }
      double getPhotonDircosu(int index) {
         if (index >= 0 && index < photonNu) return photonDircosu[index];
         else return -1;
      }
      double getPhotonDircosv(int index) {
         if (index >= 0 && index < photonNu) return photonDircosv[index];
         else return -1;
      }
      double getPhotonDircosw(int index) {
         if (index >= 0 && index < photonNu) return photonDircosw[index];
         else return -1;
      }
      double getNeutronAge(int index) {
         if (index >= 0 && index < neutronNu) return neutronAges[index];
         else return -1;
      }
      double getPhotonAge(int index) {
         if (index >= 0 && index < photonNu) return photonAges[index];
         else return -1;
      }
      static void setDelayOption(int delay) {
         delayoption = delay;
      };
      static void setCorrelationOption(int correlation) {
         correlationoption = correlation;
      };
      static void setNudistOption(int nudist) {
         nudistoption = nudist;
      };
      static void setCf252Option(int ndist, int neng) {
         Cf252ndistoption = ndist;
         Cf252nengoption = neng;
      };
      static void setRNGf(float (*funcptr) (void)) {
         rngfptr = funcptr;
         rngdptr = rngf2d;
      }
      static void setRNGd(double (*funcptr) (void)) {
         rngdptr = funcptr;
      }


   private:
      int SmpNuDistDataU232_234_236_238(double nubar);
      int SmpNuDistDataU232_234_236_238_MC(double nubar);
      int SmpNuDistDataU233_235(double nubar);
      int SmpNuDistDataU233_235_MC(double nubar);
      int SmpNuDistDataU235(double erg, int option);
      int SmpNuDistDataPu239(double erg);
      double SmpNVel(double eng, double* cosdiru, double* cosdirv, double* cosdirw);
      double SmpNEngCf252(int option);
      void SmpIsoDir(double* cosdiru, double* cosdirv, double* cosdirw);
      double SmpGEng();
      int SmpNuDistDataPu239_241(double nubar);
      int SmpNuDistDataPu239_241_MC(double nubar);
      int SmpNuDistDataU238(double erg);
      int SmpNugDist(int isotope, double nubar);
      double SmpPVel(double eng, double* cosdiru, double* cosdirv, double* cosdirw);
      int SmpSpNuDistData(int isotope, int Cf252option);
      double SmpSpNubarData(int isotope);
      int SmpSpNugDistData(int isotope);
      double SmpTerrell(double nubar);
      double SmpWatt(double ePart, int iso);
      void fissionerr(int iSever, std::string chSubNam, std::string chMsg);
      static double fisslibrng(void);
      static double rngf2d(void);
};
