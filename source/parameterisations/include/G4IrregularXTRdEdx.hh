// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IrregularXTRdEdx.hh,v 1.1 2001-02-27 15:22:49 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// Very rough model describing a radiator of X-ray transition radiation.  
// Thicknesses of plates and gas gaps are exponentially distributed.
// We suppose that:
// formation zone << mean thickness << absorption length
// for each material and in the range 1-100 keV. This allows us to simplify
// essentially interference effects in radiator stack (GetStackFactor method).
// The price is decreasing of X-ray TR photon yield. 
// 
// History:
// 27.02.01 V. Grichine first version 
//


#ifndef G4IrregularXTRdEdx_h
#define G4IrregularXTRdEdx_h 1

#include "G4VFastSimulationModel.hh"
#include "G4VXTRdEdx.hh"

class G4IrregularXTRdEdx : public G4VXTRdEdx
{
public:

   G4IrregularXTRdEdx (G4LogicalVolume *anEnvelope,G4double,G4double);
  ~G4IrregularXTRdEdx();

  // Pure virtual function from base class

  //  void DoIt(const G4FastTrack&, G4FastStep&);
  G4double GetStackFactor( G4double energy, G4double gamma, G4double varAngle);
};

#endif
