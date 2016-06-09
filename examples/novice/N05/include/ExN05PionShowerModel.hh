//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExN05PionShowerModel.hh,v 1.9 2005/11/16 07:39:03 mverderi Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
//----------------------------------------------
// Parameterisation for pi+/pi- producing hits.
//----------------------------------------------
#ifndef ExN05PionShowerModel_h
#define ExN05PionShowerModel_h 1

#include "ExN05EnergySpot.hh"

#include "G4VFastSimulationModel.hh"
#include "G4Step.hh"
#include "G4TouchableHandle.hh"
#include <vector>

class ExN05PionShowerModel : public G4VFastSimulationModel
{
public:
  //-------------------------
  // Constructor, destructor
  //-------------------------
  ExN05PionShowerModel (G4String, G4Region*);
  ExN05PionShowerModel (G4String);
  ~ExN05PionShowerModel ();

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
  
private:
  G4Step                         *fFakeStep;
  G4StepPoint                    *fFakePreStepPoint, *fFakePostStepPoint;
  G4TouchableHandle              fTouchableHandle;
  G4Navigator                    *fpNavigator;
  G4bool                         fNaviSetup;

  std::vector<ExN05EnergySpot> feSpotList;

};
#endif
