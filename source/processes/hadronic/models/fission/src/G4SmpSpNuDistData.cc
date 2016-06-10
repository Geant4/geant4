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
// $Id: G4SmpSpNuDistData.cc 68799 2013-04-05 13:29:46Z gcosmo $
//

#include "G4fissionEvent.hh"

#define nSPfissIso 8
#define nSPfissn 11

G4int G4fissionEvent::G4SmpSpNuDistData(G4int isotope, G4int Cf252option) {

/*
  Description
    Sample Number of Neutrons from spontaneous fission 
    (a) from the neutron multiplicity data for 
        U-238, Pu-238, Pu-240, Pu-242, Cm-242, Cm-244
           using Holden and Zucker's tabulated data
        Cf-252 using either Spencer's tabulated data or 
           Boldeman's data
    (b) from Terrell's approximation using nubar for
        Th-232, 
        U-232, U-233, U-234, U-235, U-236,
        Np-237, 
        Pu-239, Pu-241, 
        Am-241, 
        Bk-249
           using Ensslin's data.
*/

/*
  Input
    iso          - isotope
    Cf252option  - 0 to use Spencer's tabulated data
                   1 to use Boldeman's data
  Output
    G4SmpSpNuDistData - sampled multiplicity
                      -1 is the isotope has 
                         no multiplicity data,
                         nor any nubar data
*/
 
  G4int i, index;
  G4double sum, nubar;
  G4double r;

  static G4double sfnu [nSPfissIso][nSPfissn] = { 
    {0.0481677,0.2485215,0.4253044,0.2284094,0.0423438,0.0072533,
     0.0000000,0.0000000,0.0000000,0.0000000,0.0000000},

    {0.0631852,0.2319644,0.3333230,0.2528207,0.0986461,0.0180199,
     0.0020407,0.0000000,0.0000000,0.0000000,0.0000000},

    {0.0679423,0.2293159,0.3341228,0.2475507,0.0996922,0.0182398,
     0.0031364,0.0000000,0.0000000,0.0000000,0.0000000},

    {0.0212550,0.1467407,0.3267531,0.3268277,0.1375090,0.0373815,
     0.0025912,0.0007551,0.0001867,0.0000000,0.0000000},

    {0.0150050,0.1161725,0.2998427,0.3331614,0.1837748,0.0429780,
     0.0087914,0.0002744,0.0000000,0.0000000,0.0000000},

    {0.0540647,0.2053880,0.3802279,0.2248483,0.1078646,0.0276366,
     0.0000000,0.0000000,0.0000000,0.0000000,0.0000000},

    {0.0021100,0.0246700,0.1229000,0.2714400,0.3076300,0.1877000,
     0.0677000,0.0140600,0.0016700,0.0001000,0.0000000},

    {0.0020900,0.0262100,0.1262000,0.2752000,0.3018000,0.1846000,
     0.0668000,0.0150000,0.0021000,0.0000000,0.0000000} };

/*
  sample the spontaneous fission neutron number distribution
*/
  index = -1;

  if (isotope == 92238) index = 0;
  else if (isotope == 94240) index = 1;
  else if (isotope == 94242) index = 2;
  else if (isotope == 96242) index = 3;
  else if (isotope == 96244) index = 4;
  else if (isotope == 94238) index = 5;
  else if (isotope == 98252 && Cf252option == 0) index = 6;
  else if (isotope == 98252 && Cf252option == 1) index = 7;

  if (index != -1) { 
    r=fisslibrng();

    sum = 0.;
    for (i = 0; i < nSPfissn-1; i++) {
      sum = sum + sfnu[index][i];
      if (r <= sum || sfnu[index][i+1] == 0.) return i;
    }
    //
    // Fall through
    //
    G4cout << " Random number out of range in SmpSpNuDistData " << G4endl;
    return -1;

  } else {
// There is no full multiplicity distribution data available
// for that isotope, let's try to find a nubar for it in
// N. Ensslin, et.al., "Application Guide to Neutron
// Multiplicity Counting," LA-13422-M (November 1998)
// and use Terrell's approximation
    nubar = G4SmpSpNubarData(isotope);
    if (nubar != -1.) {
      return (G4int) G4SmpTerrell(nubar);
    } else {
// There is no nubar information for that isotope, return -1,
// meaning no data available for that isotope
      return -1;
    }
  }
}
