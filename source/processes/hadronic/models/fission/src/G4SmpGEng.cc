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
// $Id: G4SmpGEng.cc 67966 2013-03-13 09:38:38Z gcosmo $
//

#include <cmath>
#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4fissionEvent.hh"

G4double G4fissionEvent::G4SmpGEng() {

/*
  Description
    Sample energy spectrum for photons emitted by fission-induced reactions. 
    The energy spectrum of the prompt fission gamma rays is obtained from
    Maienschein's measurements.
*/

/*
  Input
  Return
    energy of photon emitted by neutron-induced fission
*/

   G4double r;
/*
  Calculate the energy of photons emitted from fission
*/
   G4Pow* Pow=G4Pow::GetInstance();
   r=fisslibrng();

   if (r == 0.0) return 0.085;

   if (r <= 0.0001) return 0.0855+0.01692*(r/0.0001)-0.02401*Pow->powA(r/0.0001,2.)+0.01274*Pow->powA(r/0.0001,3.);

   if (r > 0.0001 && r <= 0.01) return 0.09141 +  0.23846*((r-0.0001)/0.0099) 
				               -  1.75947*Pow->powA((r-.0001)/0.0099,2.)
                                               + 10.98611*Pow->powA((r-0.0001)/0.0099,3.)
                                               - 43.19181*Pow->powA((r-.0001)/0.0099,4.)
                                               +105.70005*Pow->powA((r-.0001)/.0099,5.)
                                               -160.72894*Pow->powA((r-.0001)/.0099,6.)
                                               +147.43399*Pow->powA((r-.0001)/.0099,7.)
                                               - 74.60043*Pow->powA((r-.0001)/0.0099,8.)
                                               + 15.97547*Pow->powA((r-.0001)/0.0099,9.);

   if (r > 0.01 && r <= 0.1537) return 0.14486 + 0.40914*((r-.01)/.1437)
				               - 1.28150*Pow->powA((r-0.01)/0.1437,2.)
                                               + 5.07377*Pow->powA((r-0.01)/0.1437,3.)
                                               -15.42031*Pow->powA((r-0.01)/0.1437,4.)
                                               +31.96346*Pow->powA((r-0.01)/0.1437,5.)
                                               -43.12605*Pow->powA((r-0.01)/0.1437,6.)
                                               +36.02908*Pow->powA((r-0.01)/0.1437,7.)
                                               -16.87185*Pow->powA((r-0.01)/0.1437,8.)
                                               + 3.37941*Pow->powA((r-0.01)/0.1437,9.);

   if (r > 0.1537 && r <= 0.7114) return (-1./2.3)*G4Log(0.71956*(0.1537-r)+0.50158);

   if (r > 0.7114 && r <= 1.0) return (-1./1.1)*G4Log(1.15292*(0.7114-r)+0.33287);
   //
   // Fall through 
   //
   G4cout << " Random number out of range in SmpGEng " << G4endl;
   return -1.0;
}
