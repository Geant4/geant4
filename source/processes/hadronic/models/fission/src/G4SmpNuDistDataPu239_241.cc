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
// $Id: G4SmpNuDistDataPu239_241.cc 67966 2013-03-13 09:38:38Z gcosmo $
//

#include <cmath>
#include "G4Pow.hh"
#include "G4fissionEvent.hh"

G4int G4fissionEvent::G4SmpNuDistDataPu239_241(G4double nubar) {

/*
  Description
    Sample Number of Neutrons from fission in Pu-239 and Pu-241 using 
    Zucker and Holden's tabulated data for Pu-239
    The P(nu) distribution is given as a function of the average 
    number of neutrons from fission, based on interpolation of the 
    Pu-239 data from Zucker and Holden.
*/

/*
  Input
    nubar    - average number of neutrons per fission
  Output
    G4SmpNuDistDataPu239_241  - sampled multiplicity
    
*/

  G4double pnu[9], cpnu, sum;
  G4double r;

/* 
  Check if nubar is within the range of experimental values
*/
  if(nubar >= 2.85 && nubar <= 4.25) {
/*
     Use Zucker and Holden Data
*/
     G4Pow* Pow=G4Pow::GetInstance();
     pnu[0]=-2.412937e-3*Pow->powN(nubar,3)+3.210687e-2*Pow->powN(nubar,2)-1.434037e-1*nubar+2.150733e-1;
     pnu[1]=-2.650615e-2*Pow->powN(nubar,3)+3.290389e-1*Pow->powN(nubar,2)-1.389007*nubar+2.002327;
     pnu[2]=3.232028e-2*Pow->powN(nubar,3)-3.176093e-1*Pow->powN(nubar,2)+8.605098e-1*nubar-3.411191e-1;
     pnu[3]=1.623289e-2*Pow->powN(nubar,3)-2.414705e-1*Pow->powN(nubar,2)+1.007282*nubar-9.583769e-1;
     pnu[4]=1.932275e-2*Pow->powN(nubar,3)-2.923666e-1*Pow->powN(nubar,2)+1.421383*nubar-1.924025;
     pnu[5]=-6.185679e-2*Pow->powN(nubar,3)+6.82888e-1*Pow->powN(nubar,2)-2.347653*nubar+2.647049;
     pnu[6]=1.79773e-2*Pow->powN(nubar,3)-1.60516e-1*Pow->powN(nubar,2)+5.228077e-1*nubar-5.939556e-1;
     pnu[7]=3.530038e-3*Pow->powN(nubar,4)-4.925425e-2*Pow->powN(nubar,3)+2.726784e-1*Pow->powN(nubar,2)-6.81281e-1*nubar+6.347577e-1;
     pnu[8]=2.837523e-3*Pow->powN(nubar,3)-2.678644e-2*Pow->powN(nubar,2)+8.545638e-2*nubar-9.156078e-2;

     sum=pnu[0]+pnu[1]+pnu[2]+pnu[3]+pnu[4]+pnu[5]+pnu[6]+pnu[7]+pnu[8];

     pnu[0]=pnu[0]/sum;
     pnu[1]=pnu[1]/sum;
     pnu[2]=pnu[2]/sum;
     pnu[3]=pnu[3]/sum;
     pnu[4]=pnu[4]/sum;
     pnu[5]=pnu[5]/sum;
     pnu[6]=pnu[6]/sum;
     pnu[7]=pnu[7]/sum;
     pnu[8]=pnu[8]/sum;

     r=fisslibrng();

     if(r <= pnu[0]) return 0;

     cpnu=pnu[0]+pnu[1];
     if(r <= cpnu) return 1;

     cpnu=cpnu+pnu[2];
     if(r <= cpnu) return 2;

     cpnu=cpnu+pnu[3];
     if(r <= cpnu) return 3;

     cpnu=cpnu+pnu[4];
     if(r <= cpnu) return 4;

     cpnu=cpnu+pnu[5];
     if(r <= cpnu) return 5;

     cpnu=cpnu+pnu[6];
     if(r <= cpnu) return 6;

     cpnu=cpnu+pnu[7];
     if(r <= cpnu) return 7;
     else return 8;

  } else {
/*
     Use Terrell's formula
*/
     return (G4int) G4SmpTerrell(nubar);
  }
}
