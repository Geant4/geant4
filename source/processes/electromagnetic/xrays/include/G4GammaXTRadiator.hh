//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4GammaXTRadiator.hh,v 1.1 2002-01-22 15:22:53 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// Rough process describing a radiator of X-ray transition radiation.  
// Thicknesses of plates and gas gaps are distributed according to gamma 
// distribution. x are thicknesses of plates or gas gaps:
//
// p(x) = (alpha/<x>)^alpha * x^(alpha-1) * exp(-alpha*x/<x>) / G(alpha)
//
// G(alpha) is Euler's gamma function.
// Plates have mean <x> = fPlateThick > 0 and power alpha = fAlphaPlate > 0 :
// Gas gaps have mean <x> = fGasThick > 0 and power alpha = fAlphaGas > 0 :
// We suppose that:
// formation zone ~ mean thickness << absorption length
// for each material and in the range 1-100 keV. This allows us to simplify
// interference effects in radiator stack (GetStackFactor method).
// 
// 
// History:
// 21.01.02 V. Grichine, first version 
//


#ifndef G4GammaXTRadiator_h
#define G4GammaXTRadiator_h 1

#include "G4VXTRenergyLoss.hh"

class G4GammaXTRadiator : public G4VXTRenergyLoss
{
public:

   G4GammaXTRadiator (G4LogicalVolume *anEnvelope,
                           G4double,G4double,
                           G4Material*,G4Material*,
                        G4double,G4double,G4int,
                        const G4String & processName = "XTRgammaRadiator");
  ~G4GammaXTRadiator ();

  // Pure virtual function from base class

  G4double GetStackFactor( G4double energy, G4double gamma, G4double varAngle);

private:

  G4double fAlphaPlate, fAlphaGas ;
};

#endif






