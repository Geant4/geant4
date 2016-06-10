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
// $Id: G4SmpNuDistDataU233_235.cc 67966 2013-03-13 09:38:38Z gcosmo $
//

#include <cmath>
#include "G4fissionEvent.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"

int G4fissionEvent::G4SmpNuDistDataU233_235(G4double nubar) {

  G4Pow* Pow = G4Pow::GetInstance();
/*
  Description
    Sample Number of Neutrons from fission in U-233 and U-235 using 
    Zucker and Holden's tabulated data for U-235
    The P(nu) distribution is given as a function of the average 
    number of neutrons from fission, based on interpolation of the 
    U-235 data from Zucker and Holden.
*/

/*
  Input
    nubar    - average number of neutrons per fission
  Output
    G4SmpNuDistDataU233_235  - sampled multiplicity
    
*/

  G4double pnu[8], cpnu, sum;
  G4double r;

/* 
  Check if nubar is within the range of experimental values
*/
  if(nubar >= 2.25 && nubar <= 4.0) {
/*
     Use Zucker and Holden Data
*/
     if(nubar <= 2.8738) pnu[0]=-9.279554e-02*Pow->powN(nubar,3)+8.036687e-01*Pow->powN(nubar,2)-2.342684*nubar+2.309035;
     else if(nubar > 2.8738 && nubar <= 3.4272) pnu[0]=1.50072e-2*Pow->powN(nubar,2)-1.109109e-1*nubar+2.063133e-1;
     else pnu[0]=1.498897e+3*G4Exp(-3.883864*nubar);

     if(nubar <= 3.2316) pnu[1]=3.531126e-2*Pow->powN(nubar,3)-2.787213e-1*Pow->powN(nubar,2)+5.824072e-1*nubar-1.067136e-1;
     else pnu[1]=6.574492e-2*Pow->powN(nubar,2)-5.425741e-1*nubar+1.123199;

     pnu[2]=1.274643e-2*Pow->powN(nubar,3)-1.387954e-1*Pow->powN(nubar,2)+3.264669e-1*nubar+1.77148e-1;

     pnu[3]=5.473738e-2*Pow->powN(nubar,5)-8.835826e-1*Pow->powN(nubar,4)+5.657201*Pow->powN(nubar,3)-1.802669e+1*Pow->powN(nubar,2)+2.867937e+1*nubar-1.794296e+1;

     pnu[4]=-3.591076e-2*Pow->powN(nubar,3)+3.092624e-1*Pow->powN(nubar,2)-7.184805e-1*nubar+5.649400e-1;

     if(nubar <= 2.8738) pnu[5]=1.699374e-2*Pow->powN(nubar,2)-1.069558e-3*nubar-6.981430e-2;
     else pnu[5]=2.100175e-2*Pow->powN(nubar,3)-1.705788e-1*Pow->powN(nubar,2)+5.575467e-1*nubar-6.245873e-1;

     if(nubar <= 3.0387) pnu[6]=9.431919e-7*Pow->powA(nubar,8.958848);
     else pnu[6]=4.322428e-3*Pow->powN(nubar,3)-2.094790e-2*Pow->powN(nubar,2)+4.449671e-2*nubar-4.435987e-2;

     pnu[7]=5.689084e-3*Pow->powN(nubar,4)-6.591895e-2*Pow->powN(nubar,3)+2.886861e-1*Pow->powN(nubar,2)-5.588146e-1*nubar+4.009166e-1;

     sum=pnu[0]+pnu[1]+pnu[2]+pnu[3]+pnu[4]+pnu[5]+pnu[6]+pnu[7];

     pnu[0]=pnu[0]/sum;
     pnu[1]=pnu[1]/sum;
     pnu[2]=pnu[2]/sum;
     pnu[3]=pnu[3]/sum;
     pnu[4]=pnu[4]/sum;
     pnu[5]=pnu[5]/sum;
     pnu[6]=pnu[6]/sum;
     pnu[7]=pnu[7]/sum;

     r=fisslibrng();

     if(r <= pnu[0]) return (int) 0;

     cpnu=pnu[0]+pnu[1];
     if(r <= cpnu) return (int) 1;

     cpnu=cpnu+pnu[2];
     if(r <= cpnu) return (int) 2;

     cpnu=cpnu+pnu[3];
     if(r <= cpnu) return (int) 3;

     cpnu=cpnu+pnu[4];
     if(r <= cpnu) return (int) 4;

     cpnu=cpnu+pnu[5];
     if(r <= cpnu) return (int) 5;

     cpnu=cpnu+pnu[6];
     if(r <= cpnu) return (int) 6;
     else return (int) 7;
  } else {
/*
     Use Terrell's formula
*/
     return (int) G4SmpTerrell(nubar);
  }
}
