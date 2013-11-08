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
// $Id: G4SmpSpNugDistData.cc 68799 2013-04-05 13:29:46Z gcosmo $
//

#include "G4fissionEvent.hh"

#define nSPfissg 30

G4int G4fissionEvent::G4SmpSpNugDistData(G4int isotope) {

/*
  Description
    Sample Number of Photons from spontaneous fission in 
    (a) Cf-252 using the double Poisson model from Brunson;
    (b) Th-232, 
        U-232, U-233, U-234, U-235, U-236, U-238*,
        Np-237, 
        Pu-238*, Pu-239, Pu-240*, Pu-241, Pu-242*,
        Am-241, 
        Cm-242*, Cm-244*, 
        Bk-249
        using the negative binomial distribution based on the
        spontaneous fission neutron nubar from Ensslin's 
        tabulated data or Holden and Zucker's tabulated data 
        (for isotopes denoted with asterix *).
*/

/*
  Input
    iso          - isotope
  Output
    G4SmpSpNugDistData - sampled multiplicity
                       -1 if there is no multiplicity data for that isotope
*/
 
  G4int i;
  G4double sum, nubar;
  G4double r;

  static G4double Cf252spdist [nSPfissg] = { 
         5.162699e-4,3.742057e-3,1.360482e-2,3.312786e-2,6.090540e-2,
         9.043537e-2,1.133984e-1,1.240985e-1,1.216759e-1,1.092255e-1,
         9.137106e-2,7.219960e-2,5.438060e-2,3.923091e-2,2.714690e-2,
         1.800781e-2,1.143520e-2,6.942099e-3,4.025720e-3,2.229510e-3,
         1.179602e-3,5.966936e-4,2.888766e-4,1.340137e-4,5.965291e-5,
         2.551191e-5,1.049692e-5,4.160575e-6,1.590596e-6,0.000000e+0
      };

/*
  sample the spontaneous fission photon number distribution 
*/
  nubar=0.;
  if (isotope == 98252) {
//  Cf-252 using the G4double Poisson model from Brunson;
    r=fisslibrng();

    sum = 0.;
    for (i = 0; i < nSPfissg-1; i++) {
      sum = sum + Cf252spdist[i];
      if (r <= sum || Cf252spdist[i+1] == 0.) return i;
    }
  } else if (isotope == 92238) {
/*
    using the spontaneous fission nubar from
    Holden and Zucker's tabulated data 
*/
    nubar = 1.9900002;
  } else if (isotope == 94240) {
    nubar = 2.1540006;
  } else if (isotope == 94242) {
    nubar = 2.1489998;
  } else if (isotope == 96242) {
    nubar = 2.54;
  } else if (isotope == 96244) {
    nubar = 2.7200005;
  } else if (isotope == 94238) {
    nubar = 2.2100301;
  }

  if (nubar != 0.) {
    return G4SmpNugDist(isotope, nubar);
  } else {
/*
    using the spontaneous fission nubar from
    N. Ensslin, et.al., "Application Guide to Neutron
    Multiplicity Counting," LA-13422-M (November 1998)
*/
    nubar = G4SmpSpNubarData(isotope);
    if (nubar != -1.) {
      return G4SmpNugDist(isotope, nubar);
    } else {
// There is no nubar information for that isotope, return -1,
// meaning no data available for that isotope
      return -1;
    }
  }
}
