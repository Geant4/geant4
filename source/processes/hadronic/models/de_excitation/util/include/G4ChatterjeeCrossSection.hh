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
// $Id: G4ChatterjeeCrossSection.hh 66241 2012-12-13 18:34:42Z gunter $
//
// V.Ivanchenko 13.04.2015
// 
// J.M. Quesada 22.04.2015 several fixes

#ifndef G4ChatterjeeCrossSection_h
#define G4ChatterjeeCrossSection_h 1

#include <CLHEP/Units/SystemOfUnits.h>
#include "globals.hh"
#include "G4Pow.hh"

static const G4double emax = 50*CLHEP::MeV;

// index: 0-neutron, 1-proton, 2-deuteron, 3-triton, 4-He3, 5-He4
// parameters: p0, p1, p2, landa0, landa1, mu0, mu1, nu0, nu1, nu2, ra

static const G4double paramC[6][11] = {
// n from Becchetti and Greenlees 
  { 0.,   0.,      0.,  18.57,  -22.93, 381.7, 24.31, 0.172,  -15.39, 804.8, 0.0},
// p from Bechetti and Greenlees 
  {15.72, 9.65, -449.,  0.00437,-16.58, 244.7, 0.503, 273.1, -182.4, -1.872, 0.0},
// d from Lohr and Haeberli 
  {-38.21,922.6,-2804.,-0.0323, -5.48,  336.1, 0.48,  524.3, -371.8, -5.924, 1.2},
// t from Becchetti and Greenlees 
  {-11.04,619.1,-2147., 0.0426, -10.33, 601.9, 0.37,  583.0, -546.2,  1.718, 1.2},
// 3he from  Becchetti and Greenlees 
  {-3.06, 278.5,-1389.,-0.00535,-11.16, 555.5, 0.4,   687.4, -476.3,  0.509, 1.2},
//  alpha from huizenga and igo 
  {10.95, -85.2, 1146., 0.0643, -13.96, 781.2, 0.29, -304.7, -470.0,  -8.58, 1.2}
};

class G4ChatterjeeCrossSection
{
public:

static G4double ComputePowerParameter(G4int resA, G4int idx)
  {
    return G4Pow::GetInstance()->powZ(resA, paramC[idx][6]);
  }

  static G4double ComputeCrossSection(G4double K, G4double cb, G4double resA13, G4double amu1, 
				      G4int idx, G4int Z, G4int resA)
  {
    G4double sig;
    G4double Kc  = std::min(K, emax);

    // parameterisation for neutron
    if(0 == Z) {
      G4double landa = paramC[idx][3]/resA13 + paramC[idx][4];
      G4double mu = (paramC[idx][5] + paramC[idx][6]*resA13)*resA13;
      G4double nu = std::abs((paramC[idx][7]*resA + paramC[idx][8]*resA13)*resA13 
			     + paramC[idx][9]);
      sig = landa*Kc + mu + nu/Kc;

      // parameterisation for charged 
    } else {
      //JMQ 20.04.2015 1.5 F
      G4double ec = cb;
      //G4double ec = 1.44 * Z * resZ / (1.5*resA13 + paramC[idx][10]);
      G4double ecsq = ec*ec;
      G4double p = paramC[idx][0] + paramC[idx][1]/ec + paramC[idx][2]/(ecsq); 
      G4double landa = paramC[idx][3]*resA + paramC[idx][4];
      G4double mu = paramC[idx][5]*amu1;
      G4double nu = amu1* (paramC[idx][7] + paramC[idx][8]*ec + paramC[idx][9]*ecsq);
      G4double q = landa - nu/ecsq - 2*p*ec;
      G4double r = mu + 2*nu/ec + p*ecsq;
      G4double ji= std::max(Kc, ec);
      if(Kc < ec) { sig = p*Kc*Kc + q*Kc + r; }
      else { sig = p*(Kc - ji)*(Kc - ji) + landa*Kc + mu + nu*(2 - Kc/ji)/ji; }
    }
    sig = std::max(sig, 0.0);
    return sig;
  }
};

#endif
