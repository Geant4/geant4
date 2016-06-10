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
// $Id: G4SmpSpNubarData.cc 68799 2013-04-05 13:29:46Z gcosmo $
//

#include "G4fissionEvent.hh"

#define nSPfissNubarIso 18

G4double G4fissionEvent::G4SmpSpNubarData(G4int isotope) {

/*
  Description
    Determine average number of neutrons from spontaneous fission for
        Th-232, 
        U-232, U-233, U-234, U-235, U-236, U-238
        Np-237, 
        Pu-239, Pu-240, Pu-241,  Pu-242
        Am-241, 
        Cm-242, Cm-244, 
        Bk-249,
        Cf-252
    Based on Ensslin's data.
    N. Ensslin, et.al., "Application Guide to Neutron Multiplicity Counting," 
    LA-13422-M (November 1998)
*/

/*
  Input
    iso          - isotope
  Output
    G4SmpSpNubarData - average number of neutrons
                     -1. is the isotope has 
                         no nubar data
*/
 
  G4int i;

  static G4int spzaid [nSPfissNubarIso] = {
      90232, 92232, 92233, 92234, 92235,
      92236, 92238, 93237, 94238, 94239,
      94240, 94241, 94242, 95241, 96242,
      96244, 97249, 98252 };
  static G4double spnubar [nSPfissNubarIso] = {
      2.14,  1.71, 1.76,  1.81, 1.86,
      1.91,  2.01, 2.05,  2.21, 2.16,
      2.156, 2.25, 2.145, 3.22, 2.54,
      2.72,  3.40, 3.757
      };

// Find nubar
  for (i=0; i<nSPfissNubarIso; i++) {
    if (isotope == spzaid[i]) {
      return spnubar[i];
    }
  }
// no nubar available for that isotope
  return -1.;
}
