// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8FoamXrayTRmodel.hh,v 1.1 2000-03-03 09:14:55 grichine Exp $
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


#ifndef Em8FoamXrayTRmodel_h
#define Em8FoamXrayTRmodel_h 1

#include "G4VFastSimulationModel.hh"
#include "G4ForwardXrayTR.hh"

#include "Em8XrayTRmodel.hh"

class Em8FoamXrayTRmodel : public Em8XrayTRmodel
{
public:

   Em8FoamXrayTRmodel (G4LogicalVolume *anEnvelope,G4double,G4double);
  ~Em8FoamXrayTRmodel ();

  // Pure virtual function from base class

  G4double GetStackFactor( G4double energy, G4double gamma, G4double varAngle);
};

#endif
