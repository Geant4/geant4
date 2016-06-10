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
// $Id: G4SmpNuDistDataU235.cc 67966 2013-03-13 09:38:38Z gcosmo $
//

#include <cmath>
#include "G4Pow.hh"
#include "G4fissionEvent.hh"

G4int G4fissionEvent::G4SmpNuDistDataU235(G4double erg, G4int option) {

/*
  Description
    Sample Number of Neutrons from fission in U-235 using probability
    distributions based on either
      (option 0) Zucker and Holden's tabulated data for U-235
      (option 1) Zucker and Holden's tabulated data for U-235 and
                 Gwin, Spencer and Ingle tabulated data for U-235 
                 at thermal energies
*/

/*
  Input
    erg      - incident neutron energy
    option   - 0 for sampling Zucker and Holden probability distributions
               1 for sampling probability distributions based on Zucker 
                 and Holden tabulated distributions as well as Gwin, 
                 Spencer and Ingle tabulated distributions at thermal 
                 energies
  Output
    G4SmpNuDistDataU235  - sampled multiplicity
    
*/
 
  G4double cpnu;
  G4double pnu[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  G4double eng;
  G4double r;

/* 
  Check if energy is within the range of experimental values
*/
  if (erg > 10) eng=10.;
  else eng=erg;

  r=fisslibrng();
/*
  U-235 nu distribution
*/
  G4Pow* Pow=G4Pow::GetInstance();
  if (option == 0) {
     if (eng <= 3.0) pnu[0]=0.0317223e0-9.67117e-3*eng+1.9726e-3*Pow->powN(eng,2)-2.33933e-4*Pow->powN(eng,3);
     if (eng > 3 && eng <= 7) pnu[0]=-1.24147e-2+2.52982e-2*eng-7.88108e-3*Pow->powN(eng,2)+9.10008e-4*Pow->powN(eng,3)-3.67208e-5*Pow->powN(eng,4);
     if (eng > 7 && eng <= 10) pnu[0]=6.31258e-2-1.89764e-2*eng+1.94475e-3*Pow->powN(eng,2)-6.74e-5*Pow->powN(eng,3);
     if (r <= pnu[0]) return 0;

     if (eng <= 4.0) pnu[1]=0.171707e0-0.0178305e0*eng+3.42286e-3*Pow->powN(eng,2)-2.1168e-3*Pow->powN(eng,3)+3.84226e-4*Pow->powN(eng,4)-1.44289e-5*Pow->powN(eng,5);
     if (eng > 4 && eng <= 7) pnu[1]=9.8633e-2+3.53323e-2*eng-1.15037e-2*Pow->powN(eng,2)+7.4e-4*Pow->powN(eng,3);
     if (eng > 7 && eng <= 10) pnu[1]=0.628295-0.180677*eng+1.80664e-2*Pow->powN(eng,2)-6.2015e-4*Pow->powN(eng,3);
     cpnu=pnu[0]+pnu[1];
     if (r <= cpnu) return 1;

     if (eng <= 4.0) pnu[2]=0.336199e0-1.59569e-2*eng+2.78036e-3*Pow->powN(eng,2)-1.59278e-3*Pow->powN(eng,3)+2.21742e-4*Pow->powN(eng,4);
     if (eng > 4 && eng <= 8) pnu[2]=0.229153e0+5.27561e-2*eng-1.29288e-2*Pow->powN(eng,2)+5.67233e-4*Pow->powN(eng,3)+8.06667e-6*Pow->powN(eng,4);
     if (eng > 8 && eng <= 10) pnu[2]=-0.395206e0+0.227399e0*eng-2.86051e-2*Pow->powN(eng,2)+1.08196e-3*Pow->powN(eng,3);
     cpnu=cpnu+pnu[2];
     if (r <= cpnu) return 2;

     if (eng <= 5.0) pnu[3]=0.30395461e0+0.01348261e0*eng-0.00262298e0*Pow->powN(eng,2)+1.99482407e-4*Pow->powN(eng,3);
     if (eng > 5 && eng <= 10) pnu[3]=0.10992355e0+0.09246839e0*eng-0.00885344e0*Pow->powN(eng,2)-7.60589252e-4*Pow->powN(eng,3)+1.50973591e-4*Pow->powN(eng,4)-6.20436503e-6*Pow->powN(eng,5);
     cpnu=cpnu+pnu[3];
     if (r <= cpnu) return 3;

     if (eng <= 4.0) pnu[4]=0.126946e0+1.64489e-2*eng+2.44029e-3*Pow->powN(eng,2)-2.1019e-3*Pow->powN(eng,3)+8.50104e-4*Pow->powN(eng,4)-1.10127e-4*Pow->powN(eng,5);
     if (eng > 4 && eng <= 8) pnu[4]=0.263373e0-7.47799e-2*eng+2.0588e-2*Pow->powN(eng,2)-1.55132e-3*Pow->powN(eng,3)+3.025e-5*Pow->powN(eng,4);
     if (eng > 8 && eng <= 10) pnu[4]=-0.277491e0+0.157606e0*eng-1.38467e-2*Pow->powN(eng,2)+4.20357e-4*Pow->powN(eng,3);
     cpnu=cpnu+pnu[4];
     if (r <= cpnu) return 4;

     if (eng <= 4.0) pnu[5]=0.0266793e0+9.05206e-3*eng-6.58754e-4*Pow->powN(eng,2)+6.26292e-4*Pow->powN(eng,3)-9.75958e-5*Pow->powN(eng,4);
     if (eng > 4 && eng <= 8) pnu[5]=0.0693092e0-1.46524e-2*eng+3.2841e-3*Pow->powN(eng,2)+1.50833e-4*Pow->powN(eng,3)-2.13e-5*Pow->powN(eng,4);
     if (eng > 8 && eng <= 10) pnu[5]=0.881442e0-0.271486e0*eng+3.15097e-2*Pow->powN(eng,2)-1.12095e-3*Pow->powN(eng,3);
     cpnu=cpnu+pnu[5];
     if (r <= cpnu) return 5;


     if (eng <= 4.0) pnu[6]=0.0026322e0+2.44017e-4*eng+4.55992e-4*Pow->powN(eng,2)+1.25233e-4*Pow->powN(eng,3)-2.35417e-5*Pow->powN(eng,4);
     if (eng > 4 && eng <= 8) pnu[6]=-5.3989e-3+9.48298e-3*eng-2.95864e-3*Pow->powN(eng,2)+5.43025e-4*Pow->powN(eng,3)-2.75625e-5*Pow->powN(eng,4);
     if (eng > 8 && eng <= 10) pnu[6]=0.177058-5.57839e-2*eng+6.81359e-3*Pow->powN(eng,2)-2.35568e-4*Pow->powN(eng,3);
     cpnu=cpnu+pnu[6];
     if (r <= cpnu) return 6;
     else return 7;

  } else if (option == 1) {
     if (eng <= 3.0) pnu[0]=0.0291000e0-4.836167e-3*eng-6.72500e-4*Pow->powN(eng,2)+2.076667e-4*Pow->powN(eng,3);
     if (eng > 3 && eng <= 7) pnu[0]=-1.23950e-2+2.52790e-2*eng-7.874333e-3*Pow->powN(eng,2)+9.09000e-4*Pow->powN(eng,3)-3.666667e-5*Pow->powN(eng,4);
     if (eng > 7 && eng <= 10) pnu[0]=6.328200e-2-1.903283e-2*eng+1.951500e-3*Pow->powN(eng,2)-6.766667e-5*Pow->powN(eng,3);
     if (r <= pnu[0]) return 0;

     if (eng <= 4.0) pnu[1]=0.166000e0-0.005591833e0*eng-5.624500e-3*Pow->powN(eng,2)+7.673333e-4*Pow->powN(eng,3)-2.00000e-6*Pow->powN(eng,4);
     if (eng > 4 && eng <= 7) pnu[1]=9.860600e-2+3.534733e-2*eng-1.150650e-2*Pow->powN(eng,2)+7.401667e-4*Pow->powN(eng,3);
     if (eng > 7 && eng <= 10) pnu[1]=0.628401e0-0.1807157e0*eng+1.807100e-2*Pow->powN(eng,2)-6.203333e-4*Pow->powN(eng,3);
     cpnu=pnu[0]+pnu[1];
     if (r <= cpnu) return 1;

     if (eng <= 4.0) pnu[2]=0.336200e0-1.596058e-2*eng+2.783625e-3*Pow->powN(eng,2)-1.593917e-3*Pow->powN(eng,3)+2.21875e-4*Pow->powN(eng,4);
     if (eng > 4 && eng <= 8) pnu[2]=0.2292350e0+5.26925e-2*eng-1.291067e-2*Pow->powN(eng,2)+5.650000e-4*Pow->powN(eng,3)+8.166667e-6*Pow->powN(eng,4);
     if (eng > 8 && eng <= 10) pnu[2]=0.3838230e0-3.4439e-2*eng+6.0800e-4*Pow->powN(eng,2);
     cpnu=cpnu+pnu[2];
     if (r <= cpnu) return 2;

     if (eng <= 4.0) pnu[3]=0.3074000e0+0.00794125e0*eng-0.0002580417e0*Pow->powN(eng,2)-1.875000e-5*Pow->powN(eng,3)-2.145833e-5*Pow->powN(eng,4);
     if (eng > 4 && eng <= 7) pnu[3]=0.3152270e0-2.623667e-3*eng+2.785000e-3*Pow->powN(eng,2)-3.273333e-4*Pow->powN(eng,3);
     if (eng > 7 && eng <= 10) pnu[3]=0.6476430e0-0.1046148e0*eng+1.181600e-2*Pow->powN(eng,2)-5.051667e-4*Pow->powN(eng,3);
     cpnu=cpnu+pnu[3];
     if (r <= cpnu) return 3;

     if (eng <= 4.0) pnu[4]=0.133300e0+5.853750e-3*eng+6.200875e-3*Pow->powN(eng,2)-8.95250e-4*Pow->powN(eng,3)+1.36250e-5*Pow->powN(eng,4);
     if (eng > 4 && eng <= 7) pnu[4]=0.2379650e0-5.548167e-2*eng+1.517350e-2*Pow->powN(eng,2)-8.858333e-4*Pow->powN(eng,3);
     if (eng > 7 && eng <= 10) pnu[4]=-0.5408690e0+0.2461313e0*eng-2.372350e-2*Pow->powN(eng,2)+7.861667e-4*Pow->powN(eng,3);
     cpnu=cpnu+pnu[4];
     if (r <= cpnu) return 4;

     if (eng <= 4.0) pnu[5]=0.025900e0+1.067450e-2*eng-1.794000e-3*Pow->powN(eng,2)+9.50500e-4*Pow->powN(eng,3)-1.3000e-4*Pow->powN(eng,4);
     if (eng > 4 && eng <= 7) pnu[5]=0.0871960e0-2.823683e-2*eng+7.0955e-3*Pow->powN(eng,2)-3.176667e-4*Pow->powN(eng,3);
     if (eng > 7 && eng <= 10) pnu[5]=-0.591650e0+0.2236360e0*eng-2.373100e-2*Pow->powN(eng,2)+9.25000e-4*Pow->powN(eng,3);
     cpnu=cpnu+pnu[5];
     if (r <= cpnu) return 5;

     if (eng <= 4.0) pnu[6]=0.002100e0+1.35500e-3*eng-3.235833e-4*Pow->powN(eng,2)+3.48500e-4*Pow->powN(eng,3)-4.591667e-5*Pow->powN(eng,4);
     if (eng > 4 && eng <= 8) pnu[6]=1.767200e-2-8.055667e-3*eng+1.96650e-3*Pow->powN(eng,2)-6.283333e-5*Pow->powN(eng,3);
     if (eng > 8 && eng <= 10) pnu[6]=-0.2485310e0+8.72590e-2*eng-9.14550e-3*Pow->powN(eng,2)+3.555000e-4*Pow->powN(eng,3);
     cpnu=cpnu+pnu[6];
     if (r <= cpnu) return 6;
     else return 7;
  }
  //
  // Fall through
  //

  G4cout << " SmpNuDistDataU235: unrecognized option = " << option << G4endl;
  return -1;

}
