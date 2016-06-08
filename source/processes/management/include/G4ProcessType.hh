// This code implementation is the intellectual property of
// the  GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcessType.hh,v 1.3 1999/11/07 17:11:45 kurasige Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
//
//---------------------------------------------------------------
//
// G4ProcessType.hh
//
// Class Description:
//   This is an enumerator to define process type
//
//
//---------------------------------------------------------------

#ifndef G4ProcessType_h
#define G4ProcessType_h 1

enum G4ProcessType
{
  fNotDefined,
  fTransportation,
  fElectromagnetic,
  fOptical,             
  fHadronic,
  fPhotolepton_hadron,
  fDecay,
  fGeneral,
  fParameterisation,
  fUserDefined
};
#endif


