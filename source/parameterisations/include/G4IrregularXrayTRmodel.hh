// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IrregularXrayTRmodel.hh,v 1.2 2001-02-27 07:37:18 grichine Exp $
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
// 23.01.00 V. Grichine first version based on ExN05PiModel class
// 08.02.00 V. Grichine, DoIt was placed in base class
//


#ifndef G4IrregularXrayTRmodel_h
#define G4IrregularXrayTRmodel_h 1

#include "G4VFastSimulationModel.hh"
// #include "G4ForwardXrayTR.hh"

#include "G4VXrayTRadModel.hh"

class G4IrregularXrayTRmodel : public G4VXrayTRadModel
{
public:

   G4IrregularXrayTRmodel (G4LogicalVolume *anEnvelope,G4double,G4double);
  ~G4IrregularXrayTRmodel ();

  // Pure virtual function from base class

  //  void DoIt(const G4FastTrack&, G4FastStep&);
  G4double GetStackFactor( G4double energy, G4double gamma, G4double varAngle);
};

#endif
