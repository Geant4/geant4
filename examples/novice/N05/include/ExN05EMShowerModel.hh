// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05EMShowerModel.hh,v 1.1 1999-01-07 16:06:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//----------------------------------------------
// Parameterisation of e+/e-/gamma producing hits
// The hits are the same as defined in the detailed
// simulation.
//----------------------------------------------
#ifndef ExN05EMShowerModel_h
#define ExN05EMShowerModel_h 1

#include "ExN05EnergySpot.hh"

#include "G4VFastSimulationModel.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include <rw/tvordvec.h>

class ExN05EMShowerModel : public G4VFastSimulationModel
{
public:
  //-------------------------
  // Constructor, destructor
  //-------------------------
  ExN05EMShowerModel (G4String, G4LogicalVolume*);
  ExN05EMShowerModel (G4String);
  ~ExN05EMShowerModel ();

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

private:
  void AssignSpotAndCallHit(const ExN05EnergySpot &eSpot);
  void FillFakeStep(const ExN05EnergySpot &eSpot);
  void Explode(const G4FastTrack&);
  void BuildDetectorResponse();
  G4Step      *fFakeStep;
  G4StepPoint *fFakePreStepPoint, *fFakePostStepPoint;
  G4VTouchable*fpTouchable;
  G4Navigator *fpNavigator;
  G4bool       fNaviSetup;

  RWTValOrderedVector<ExN05EnergySpot> feSpotList;

};
#endif
