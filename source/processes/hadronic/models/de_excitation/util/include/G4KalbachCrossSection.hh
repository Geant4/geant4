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
// $Id: G4KalbachCrossSection.hh 66241 2012-12-13 18:34:42Z gunter $
//
// V.Ivanchenko 13.04.2015
// 
// J.M. Quesada 22.04.2015 several fixes

#ifndef G4KalbachCrossSection_h
#define G4KalbachCrossSection_h 1

#include "globals.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"

//from subroutine sigpar of PRECO-2000  by Constance Kalbach Walker
//     Calculate optical model reaction cross sections
//     using the empirical parameterization
//     of Narasimha Murthy, Chaterjee, and Gupta
//     going over to the geometrical limit at high energy.
//
//           Proton cross sections scaled down with signor for a<100
//           (appropriate for becchetti-greenlees potential).
//           p2 reduced and global red'n factor introduced below Bc
//           Neutron cross sections scaled down with signor for a<40
//           Scaled up for A>210 (added June '98 to conform with
//           my published papers)
//           (appropriate for Mani et al potential)
//

// index: 0-neutron, 1-proton, 2-deuteron, 3-triton, 4-He3, 5-He4
// parameters: p0, p1, p2, lambda0, lambda1, mu0, mu1, nu0, nu1, nu2, ra

static const G4double paramK[6][11] = {
// n from mani, melkanoff and iori
  {-312., 0.,      0.,  12.10,  -11.27, 234.1, 38.26, 1.55,  -106.1, 1280.8, 0.0},
// p from  becchetti and greenlees (but modified with sub-barrier
// correction function and p2 changed from -449)
  {15.72, 9.65, -300.,  0.00437,-16.58, 244.7, 0.503, 273.1, -182.4, -1.872, 0.0},
// d from o.m. of perey and perey
  {0.798, 420.3,-1651., 0.00619, -7.54, 583.5, 0.337, 421.8, -474.5, -3.592, 0.8},
// t from o.m. of hafele, flynn et al
  {-21.45,484.7,-1608., 0.0186, -8.9,   686.3, 0.325, 368.9, -522.2, -4.998, 0.8},
// 3he from o.m. of gibson et al
  {-2.88,205.6, -1487.,0.00459,-8.93,   611.2, 0.35 , 473.8, -468.2, -2.225, 0.8},
// alpha from huizenga and igo 
  { 10.95,-85.2, 1146.,  0.0643,-13.96, 781.2,  0.29, -304.7,-470.0, -8.580, 1.2}
};

class G4KalbachCrossSection
{
public:

static G4double ComputePowerParameter(G4int resA, G4int idx)
  {
    return G4Pow::GetInstance()->powZ(resA, paramK[idx][6]);
  }

  static G4double ComputeCrossSection(G4double K, G4double cb,  
				      G4double resA13, G4double amu1, 
				      G4int idx, G4int Z, G4int A, 
				      G4int resA)
  {    
    G4double sig = 0.0;
    G4double signor = 1.0; 
    G4double lambda, mu, nu;
    G4double ec = 0.5;
    if(0 < Z) { ec = cb; }
      //JMQ 13.02.2009 tuning for improving cluster emission ddxs 
      //               (spallation benchmark) 
    /*
      G4double xx = 1.7;
      if(1 == A) { xx = 1.5; }
      ec = 1.44 * Z * resZ / (xx*resA13 + paramK[idx][10]);
    }
    */
    G4double ecsq = ec*ec;
    G4double elab = K * (A + resA) / G4double(resA);    
    
    if(idx == 0) { // parameterization for neutron

      if(resA < 40)       { signor =0.7 + resA*0.0075; }
      else if(resA > 210) { signor = 1. + (resA-210)*0.004; }
      lambda = paramK[idx][3]/resA13 + paramK[idx][4];
      mu = (paramK[idx][5] + paramK[idx][6]*resA13)*resA13;
      // JMQ 20.11.2008 very low energy behaviour corrected 
      //                (problem for A (apprx.)>60) fix for avoiding  
      //                neutron xs going to zero at very low energies
      nu = std::abs((paramK[idx][7]*resA + paramK[idx][8]*resA13)*resA13 
		    + paramK[idx][9]);

    } else { // parameterization for charged 
      // proton correction
      if(idx == 1) {
	if (resA <= 60)      { signor = 0.92; }
	else if (resA < 100) { signor = 0.8 + resA*0.002; }
      }
      lambda = paramK[idx][3]*resA + paramK[idx][4];
      mu = paramK[idx][5]*amu1;
      nu = amu1* (paramK[idx][7] + paramK[idx][8]*ec + paramK[idx][9]*ecsq);
    }
    /*
    G4cout << "## idx= " << idx << " K= " << K << " elab= " << elab << "  ec= " << ec 
	   << " lambda= " << lambda << " mu= " << mu << "  nu= " << nu << G4endl; 
    */
    // threashold cross section
    if(elab < ec) {
      G4double p = paramK[idx][0];
      if(0 < Z) { p += paramK[idx][1]/ec + paramK[idx][2]/ecsq; }
      G4double a = -2*p*ec + lambda - nu/ecsq;
      G4double b = p*ecsq + mu + 2*nu/ec;
      G4double ecut;
      G4double det = a*a - 4*p*b;
      if (det > 0.0) { ecut = (std::sqrt(det) - a)/(2*p); }
      else           { ecut = -a/(2*p); }

      //G4cout << "  elab= " << elab << " ecut= " << ecut << " sig= " << sig
      //	     << "  sig1= " << (p*elab*elab + a*elab + b)*signor << G4endl;
      // If ecut>0, sig=0 at elab=ecut
      if(0 == idx) {
      	sig = (lambda*ec + mu + nu/ec)*signor*std::sqrt(elab/ec);
      } else if(elab >= ecut) { 
	sig = (p*elab*elab + a*elab + b)*signor; 

	// extra proton correction
        if(1 == idx) {
	  // c and w are for global correction factor for 
	  // they are scaled down for light targets where ec is low.
	  G4double cc = std::min(3.15, ec*0.5);
	  G4double signor2 = (ec - elab - cc) *3.15/ (0.7*cc);
	  sig /= (1. + G4Exp(signor2));
	}
      }
      //G4cout << "       ecut= " << ecut << " a= " << a << "  b= " << b 
      //     <<  " signor= " << signor << " sig= " << sig << G4endl;  

      // high energy cross section
    } else {
      // etest is the energy above which the rxn cross section is
      // compared with the geometrical limit and the max taken.

      // neutron parameters
      G4double etest  = 32.;
      G4double xnulam = 1.0;

      // parameters for charged
      static const G4double flow = 1.e-18;
      static const G4double spill= 1.e+18;
      if(0 < Z) {
        etest = 0.0;
	xnulam = nu / lambda;
	xnulam = std::min(xnulam, spill); 
	if (xnulam >= flow) { 
	  if(1 == idx) { etest = std::sqrt(xnulam) + 7.; }
	  else         { etest = 1.2 *std::sqrt(xnulam); }
	}
      }
      // ** For xnulam.gt.0, sig reaches a maximum at sqrt(xnulam).
      sig = (lambda*elab + mu + nu/elab)*signor;
      if (xnulam >= flow && elab >= etest) {
	G4double geom = std::sqrt(A*K);
	geom = 1.23*resA13 + paramK[idx][10] + 4.573/geom;
	geom = 31.416 * geom * geom;
	sig = std::max(sig, geom);
      }
    }
    sig = std::max(sig, 0.0);
    //G4cout << "  ---- sig= " << sig << G4endl;
    return sig;
  }
};

#endif
