// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05PiModel.hh,v 1.1 1999-01-07 16:06:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//----------------------------------------
// Simple example of a "parameterisation".
//----------------------------------------
#ifndef ExN05PiModel_h
#define ExN05PiModel_h 1

#include "G4VFastSimulationModel.hh"

class ExN05PiModel : public G4VFastSimulationModel
{
public:
  //-------------------------
  // Constructor, destructor
  //-------------------------
  ExN05PiModel (G4LogicalVolume *anEnvelope);
  ~ExN05PiModel ();

  //------------------------------
  // Virtual methods of the base
  // class to be coded by the user
  //------------------------------

  // -- IsApplicable
  G4bool IsApplicable(const G4ParticleDefinition&);
  // -- ModelTrigger
  G4bool ModelTrigger(const G4FastTrack &);
  // -- User method DoIt
  void DoIt(const G4FastTrack&, G4FastStep&);

};
#endif
