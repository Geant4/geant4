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
// $Id: GFlashSamplingShowerTuning.hh 68057 2013-03-13 14:46:00Z gcosmo $
//
//
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  GFlashSamplingShowerTuning
//
//  Class description:
//
//  Tuning class for GFlash homogeneous shower parameterisation.
//  Definitions:
//    <t>: shower center of gravity
//      T: Depth at shower maximum
//     Ec: Critical energy
//     X0: Radiation length
//     y = E/Ec
//
//  Please, see hep-ex/0001020 for details.

// Author: Joanna Weng - 11.2005
//---------------------------------------------------------------
#ifndef GFlashSamplingShowerTuning_hh
#define GFlashSamplingShowerTuning_hh

#include "GVFlashHomoShowerTuning.hh"

class GFlashSamplingShowerTuning : public GVFlashHomoShowerTuning
{
  public:

    GFlashSamplingShowerTuning() {}
    virtual ~GFlashSamplingShowerTuning() {}
  

  public: // with description

  G4double ParsAveT1(){ return -0.55;} // t1
  G4double ParsAveT2(){ return -0.69;} // t2
    // T_sam =  log(exp( log T_hom) + t1*Fs-1 + t2*(1-ehat))

  G4double ParsAveA1(){ return -0.476;  } // a1
    // alpha_sam = log(exp(log alphah_hom) +(a1*Fs-1))

  G4double ParsSigLogT1(){ return -2.5;} // t1
  G4double ParsSigLogT2(){ return 1.25;} // t2
    // std::sqrt(var(ln(T_sam))) = 1/(t+t2*ln(y))

  G4double ParsSigLogA1(){ return -0.82;} // a1
  G4double ParsSigLogA2(){ return 0.79; } // a2
    // std::sqrt(var(ln(alpha_sam))) = 1/(a1+a2*ln(y))

  G4double ParsRho1(){ return 0.784; } // r1
  G4double ParsRho2(){ return -0.023;} // r2
    // Correlation(ln(T),ln(alpha))=r1+r2*ln(y)

  // Radial profiles
  // f(r) := (1/dE(t))(dE(t,r)/dr)
  // Ansatz:
  // f(r) = p(2*r*Rc**2)/(r**2+Rc**2)**2+(1-p)*(2*r*Rt**2)/(r**2+Rt**2)**2,
  //        0<p<1

  G4double ParsRC1(){ return -0.0203;   } // c1
  G4double ParsRC2(){ return 0.0397;  }   // c2
    // Rc_sam = Rc_hom + c1 * (1-ehat) + c2 *Fs-1*exp (-tau)

  G4double ParsRT1(){ return -0.14;  }   // t1
  G4double ParsRT2(){ return -0.495; }   // t2
    // Rt_sam = Rc_hom + t1 * (1-ehat) + t2 *Fs-1*exp (-tau)

  G4double ParsWC1(){ return 0.348;   } // c1
  G4double ParsWC2(){ return -0.642;} // c2
    // W_sam = W_hom + (1-ehat)*(c1 + c2 *Fs-1 * exp (- (tau -1 )**2))

  // Fluctuations on radial profiles through number of spots
  // The total number of spots needed for a shower is

  G4double ParsSpotN1(){ return 10.3; } // n1
  G4double ParsSpotN2(){ return 0.959;} // n2
    // Ns = n1*ln(Z)(E/GeV)**n2

  // The number of spots per longitudinal interval is:
  // (1/Ns)(dNs(t)/dt) = f(t)
  //  = (beta*t)**(alpha-1)*beta*std::exp(-beta*t)/Gamma(alpha)
  // <t> = alpha_s/beta_s
  // Ts = (alpha_s-1)/beta_s
  // and
  // Ts = T*(t1+t2*Z)
  // alpha_s = alpha*(a1+a2*Z)

  G4double ParsSpotT1(){ return 0.813; } // t1
  G4double ParsSpotT2(){ return 0.0019;} // t2

  G4double ParsSpotA1(){ return 0.844; } //a1
  G4double ParsSpotA2(){ return 0.0026;} //a2

  // Resolution

  G4double ConstantResolution(){ return 0.00;  }  
  G4double NoiseResolution()   { return 0.00;  } // not used    
  G4double SamplingResolution(){ return 0.11;  } // not used

};

#endif
