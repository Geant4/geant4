// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FoamXTRdEdx.hh,v 1.1 2001-02-27 15:22:49 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// Rough model describing a radiator of X-ray transition radiation.  
// Thicknesses of plates and gas gaps are exponentially distributed.
// We suppose that:
// formation zone ~ mean thickness << absorption length
// for each material and in the range 1-100 keV. This allows us to simplify
// interference effects in radiator stack (GetStackFactor method).
// 
// 
// History:
// 27.02.01 V. Grichine, first version 
//


#ifndef G4FoamXTRdEdx_h
#define G4FoamXTRdEdx_h 1

#include "G4VFastSimulationModel.hh"
// #include "G4ForwardXrayTR.hh"

#include "G4VXTRdEdx.hh"

class G4FoamXTRdEdx : public G4VXTRdEdx
{
public:

   G4FoamXTRdEdx (G4LogicalVolume *anEnvelope,G4double,G4double);
  ~G4FoamXTRdEdx ();

  // Pure virtual function from base class

  G4double GetStackFactor( G4double energy, G4double gamma, G4double varAngle);
};

#endif
