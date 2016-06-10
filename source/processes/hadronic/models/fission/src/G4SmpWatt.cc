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
// $Id: G4SmpWatt.cc 68799 2013-04-05 13:29:46Z gcosmo $
//

#define nZAfis 39    /* 38 fissionable isotopes in ENDL + U-232 from ENDF7 */
#define WATTEMIN 1.0e-6
#define WATTEMAX 20.0

#include <cmath>
#include "G4Log.hh"
#include "G4fissionEvent.hh"

G4double G4fissionEvent::G4SmpWatt(G4double ePart, G4int iso) {

/*
  Description
    Sample Watt Spectrum as in TART (Kalos algorithm)
*/

/*
  Input
    ePart     - energy of incoming particle
    iso       - isotope
  Output
              - energy of incoming particle
*/

  static G4int nZA [nZAfis]= {
                 90231, 90232, 90233,
                 91233,
                 92232, 92233, 92234, 92235, 92236, 92237, 92238, 92239, 92240,
                 93235, 93236, 93237, 93238,
                 94237, 94238, 94239, 94240, 94241, 94242, 94243,
                 95241, 95242, 95243,
                 96242, 96243, 96244, 96245, 96246, 96247, 96248,
                 97249,
                 98249, 98250, 98251, 98252};

  static G4double Watta [nZAfis][3] = {
                      {6.00949285e-05, -8.36695381e-03,  9.50939496e-01},
                      {6.54348443e-05, -8.86574327e-03,  9.55404490e-01},
                      {7.08173682e-05, -9.22676286e-03,  9.50088329e-01},
                      {6.35839062e-05, -8.63645973e-03,  9.24583535e-01},
                      {8.21929628e-05,  4.01922936e-03,  1.152121164e00},
                      {6.21335718e-05, -8.45651858e-03,  9.14717276e-01},
                      {6.81386135e-05, -8.99142394e-03,  9.21954824e-01},
                      {7.32627297e-05, -9.36908697e-03,  9.20107976e-01},
                      {8.06505279e-05, -9.95416671e-03,  9.27890410e-01},
                      {8.33208285e-05, -1.01073057e-02,  9.17691654e-01},
                      {8.96944680e-05, -1.06491070e-02,  9.25496030e-01},
                      {9.44608097e-05, -1.08940419e-02,  9.17795511e-01},
                      {1.01395704e-04, -1.15098159e-02,  9.29395462e-01},
                      {6.81110009e-05, -8.91619352e-03,  9.00047566e-01},
                      {7.21126359e-05, -9.20179363e-03,  8.95722889e-01},
                      {7.82371142e-05, -9.67050621e-03,  8.99574933e-01},
                      {8.27256297e-05, -9.99353009e-03,  8.97461897e-01},
                      {7.29458059e-05, -9.22415170e-03,  8.80996165e-01},
                      {8.02383914e-05, -9.78291439e-03,  8.88964070e-01},
                      {8.50641730e-05, -1.01099145e-02,  8.87304833e-01},
                      {9.10537157e-05, -1.05303084e-02,  8.89438514e-01},
                      {9.43014320e-05, -1.07133543e-02,  8.82632055e-01},
                      {1.02655616e-04, -1.13154691e-02,  8.91617174e-01},
                      {1.06118094e-04, -1.14971777e-02,  8.85181637e-01},
                      {9.08474473e-05, -1.04296303e-02,  8.71942958e-01},
                      {9.35633054e-05, -1.05612167e-02,  8.63930371e-01},
                      {1.01940441e-04, -1.11573929e-02,  8.73153437e-01},
                      {9.19501202e-05, -1.04229157e-02,  8.58681822e-01},
                      {9.42991674e-05, -1.05098872e-02,  8.49103546e-01},
                      {1.02747171e-04, -1.11371417e-02,  8.60434431e-01},
                      {1.05024967e-04, -1.12138980e-02,  8.51101942e-01},
                      {1.14130011e-04, -1.18692049e-02,  8.62838259e-01},
                      {1.15163673e-04, -1.18553822e-02,  8.51306646e-01},
                      {1.27169055e-04, -1.27033210e-02,  8.68623539e-01},
                      {1.24195213e-04, -1.24047085e-02,  8.48974077e-01},
                      {1.12616150e-04, -1.15135023e-02,  8.19708800e-01},
                      {1.23637465e-04, -1.22869889e-02,  8.35392018e-01},
                      {1.22724317e-04, -1.21677963e-02,  8.22569523e-01},
                      {1.33891595e-04, -1.29267762e-02,  8.37122909e-01} };

   G4double a;  /* Watt Parameters */
   G4double b=1.0;

   G4double rand1,rand2;
   G4double x,y,z;
   G4double eSmp;
   G4int i;


/*
   Find Watt parameters for isotope
*/
   G4int isoindex=-1;
   for (i=0; isoindex == -1 && i<nZAfis; i++) {
      if (iso == nZA[i]) isoindex = i;
   }
   if (isoindex == -1) {
      std::ostringstream o;
      o << iso;
      std::string errMsg = "No Watt spectrum available for iso " + o.str();
      G4fissionerr(6, "SmpWatt", errMsg);
   }
   
   a= Watta[isoindex][2] + ePart*(Watta[isoindex][1] + ePart*Watta[isoindex][0]);

   x= 1. + (b/(8.*a));
   y= (x + std::sqrt(x*x-1.))/a;
   z= a*y - 1.;

   G4int icounter = 0;
   G4int icounter_max = 1024;
   do {

      rand1= -G4Log(fisslibrng());
      rand2= -G4Log(fisslibrng());
      eSmp= y*rand1;

      icounter++;
      if ( icounter > icounter_max ) { 
	 G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
         break;
      }

   } while ((rand2-z*(rand1+1.))*(rand2-z*(rand1+1.)) > b*y*rand1 ||
             eSmp < WATTEMIN || eSmp > WATTEMAX);
   // Loop checking, 11.03.2015, T. Koi
   
   return eSmp;
}
