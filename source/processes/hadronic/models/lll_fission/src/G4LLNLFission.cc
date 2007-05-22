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
// $Id: G4LLNLFission.cc,v 1.1 2007-05-22 00:52:06 dennis Exp $
//
// This class is a copy of Fission.cc, made for use with Geant4.
//

#include "fissionEvent.hh"
#include <stdio.h>
#include <stdlib.h>

fissionEvent* fe;

extern "C" {
   extern float (*rngfptr) (void);

   extern double (*rngdptr) (void);

   extern double rngf2d(void);

   void genspfissevt_(int *isotope, double *time) {
      if (fe != 0) delete fe;
      fe = new fissionEvent(*isotope, *time, -1., 0.);
   };

   void genfissevt_(int *isotope, double *time, double *nubar, double *eng) {
      if (fe != 0) delete fe;
      fe = new fissionEvent(*isotope, *time, *nubar, *eng);
   };

   int getnnu_() {
      return (*fe).getNeutronNu();
   };

   int getpnu_() {
      return (*fe).getPhotonNu();
   };

   double getneng_(int *index) {
      return (*fe).getNeutronEnergy(*index);
   };

   double getnvel_(int *index) {
      return (*fe).getNeutronVelocity(*index);
   };

   double getndircosu_(int *index) {
      return (*fe).getNeutronDircosu(*index);
   };

   double getndircosv_(int *index) {
      return (*fe).getNeutronDircosv(*index);
   };

   double getndircosw_(int *index) {
      return (*fe).getNeutronDircosw(*index);
   };

   double getpeng_(int *index) {
      return (*fe).getPhotonEnergy(*index);
   };

   double getpvel_(int *index) {
      return (*fe).getPhotonVelocity(*index);
   };

   double getpdircosu_(int *index) {
      return (*fe).getPhotonDircosu(*index);
   };

   double getpdircosv_(int *index) {
      return (*fe).getPhotonDircosv(*index);
   };

   double getpdircosw_(int *index) {
      return (*fe).getPhotonDircosw(*index);
   };

   double getnage_(int *index) {
      return (*fe).getNeutronAge(*index);
   };

   double getpage_(int *index) {
      return (*fe).getPhotonAge(*index);
   };

   void setdelay_(int *delay) {
      (*fe).setDelayOption(*delay);
   };

   void setcorrel_(int *correlation) {
      (*fe).setCorrelationOption(*correlation);
   };

   void setnudist_(int *nudist) {
/*
      where the argument *nudist affects induced fissions only, it
      is set to
         0 for sampling Zucker and Holden probability distributions 
           for U-235,238 and Pu-239. Terrell for other isotopes.
         1 same as above, but using Gwin, Spencer and Ingle 
           tabulated distributions for thermal energies for U-235.
           Terrell for other isotopes.
         2 for sampling fission-induced neutron multiplicity in 
           (a) U-232, U-234, U-236 and U-238 using Zucker and 
               Holden's tabulated data for U-238
           (b) U-233 and U-235 using Zucker and Holden's tabulated 
               data for U-235
           (c) Pu-239 and Pu-241 using Zucker and Holden's tabulated 
               data for Pu-239
           The P(nu) distributions for *nudist=2 are given as a 
           function of the average number of neutrons from fission, 
           based on interpolation of the data from Zucker and Holden.
           Terrell for other isotopes.
         3 for sampling fission-induced neutron multiplicity in 
           (a) U-232, U-234, U-236 and U-238 using Zucker and 
               Holden's tabulated data for U-238
           (b) U-233 and U-235 using Zucker and Holden's tabulated 
               data for U-235
           (c) Pu-239 and Pu-241 using Zucker and Holden's tabulated 
               data for Pu-239
           The Z&H tables have P(nu) distributions for 11 energies 
           (0 MeV through 10 MeV), along with their nubars. For 
           *nudist=3, we select the P(nu) distribution that has
           a nubar closest either from above, or from below, to the 
           to the nubar entered for the induced fission, based on a 
           random number and fractional distances to the end of the 
           nubar interval thus formed.
           Terrell for other isotopes.
*/

      (*fe).setNudistOption(*nudist);
   };

   void setcf252_(int *ndist, int *neng) {
/*
      where the argument
      *ndist is set to 
         0 to sample the spontaneous fission neutron multiplicity 
           using tabulated data from Spencer
         1 to sample the spontaneous fission neutron multiplicity
           using tabulated data from Boldeman
      *neng is set to 
         0 to sample the Mannhart corrected Maxwellian spectrum
         1 to sample the Madland-Nix theoretical spectrum
         2 to sample the Froehner Watt spectrum
*/
      (*fe).setCf252Option(*ndist, *neng);
   }; 

   void setrngf_(float (*funcptr) (void)) {
      fissionEvent::setRNGf(funcptr);
   }

   void setrngd_(double (*funcptr) (void)) {
      fissionEvent::setRNGd(funcptr);
   }
}
