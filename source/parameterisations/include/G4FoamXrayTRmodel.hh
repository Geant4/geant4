// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FoamXrayTRmodel.hh,v 1.2 2001-02-27 07:36:47 grichine Exp $
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
// 09.02.00 V. Grichine, first version 
//


#ifndef G4FoamXrayTRmodel_h
#define G4FoamXrayTRmodel_h 1

#include "G4VFastSimulationModel.hh"
// #include "G4ForwardXrayTR.hh"

#include "G4VXrayTRadModel.hh"

class G4FoamXrayTRmodel : public G4VXrayTRadModel
{
public:

   G4FoamXrayTRmodel (G4LogicalVolume *anEnvelope,G4double,G4double);
  ~G4FoamXrayTRmodel ();

  // Pure virtual function from base class

  G4double GetStackFactor( G4double energy, G4double gamma, G4double varAngle);
};

#endif
