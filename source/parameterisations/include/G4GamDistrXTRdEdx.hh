// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GamDistrXTRdEdx.hh,v 1.1 2001-02-27 15:22:49 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// Rough model describing a radiator of X-ray transition radiation.  
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
// 27.02.01 V. Grichine, first version 
//


#ifndef G4GamDistrXTRdEdx_h
#define G4GamDistrXTRdEdx_h 1

#include "G4VFastSimulationModel.hh"
// #include "G4ForwardXrayTR.hh"

#include "G4VXTRdEdx.hh"

class G4GamDistrXTRdEdx : public G4VXTRdEdx
{
public:

   G4GamDistrXTRdEdx (G4LogicalVolume *anEnvelope,
                           G4double,G4double,
                           G4double,G4double     );
  ~G4GamDistrXTRdEdx ();

  // Pure virtual function from base class

  G4double GetStackFactor( G4double energy, G4double gamma, G4double varAngle);

private:

  G4double fAlphaPlate, fAlphaGas ;
};

#endif
