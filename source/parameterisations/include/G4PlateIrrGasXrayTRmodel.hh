// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PlateIrrGasXrayTRmodel.hh,v 1.2 2001-02-27 07:37:28 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// Model describing a radiator of X-ray transition radiation.  
// Thicknesses of plates is fixed while gas gaps are fully irregular.
// We suppose that:
// formation zone ~ mean thickness << absorption length
// for each material and in the range 1-100 keV. This allows us to simplify
// interference effects in radiator stack (GetStackFactor method).
// 
// 
// History:
// 10.02.00 V. Grichine, first version 
//


#ifndef G4PlateIrrGasXrayTRmodel_h
#define G4PlateIrrGasXrayTRmodel_h 1

#include "G4VFastSimulationModel.hh"
// #include "G4ForwardXrayTR.hh"

#include "G4VXrayTRadModel.hh"

class G4PlateIrrGasXrayTRmodel : public G4VXrayTRadModel
{
public:

   G4PlateIrrGasXrayTRmodel (G4LogicalVolume *anEnvelope,G4double,G4double);
  ~G4PlateIrrGasXrayTRmodel ();

  // Pure virtual function from base class

  G4double GetStackFactor( G4double energy, G4double gamma, G4double varAngle);
};

#endif
