// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8GamDistrXrayTRmodel.hh,v 1.1 2000-03-03 09:15:41 grichine Exp $
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
// 11.02.00 V. Grichine, first version 
//


#ifndef Em8GamDistrXrayTRmodel_h
#define Em8GamDistrXrayTRmodel_h 1

#include "G4VFastSimulationModel.hh"
#include "G4ForwardXrayTR.hh"

#include "Em8XrayTRmodel.hh"

class Em8GamDistrXrayTRmodel : public Em8XrayTRmodel
{
public:

   Em8GamDistrXrayTRmodel (G4LogicalVolume *anEnvelope,
                           G4double,G4double,
                           G4double,G4double     );
  ~Em8GamDistrXrayTRmodel ();

  // Pure virtual function from base class

  G4double GetStackFactor( G4double energy, G4double gamma, G4double varAngle);

private:

  G4double fAlphaPlate, fAlphaGas ;
};

#endif
