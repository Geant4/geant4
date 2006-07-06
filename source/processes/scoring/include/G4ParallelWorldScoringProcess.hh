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
// $Id: G4ParallelWorldScoringProcess.hh,v 1.1 2006-07-06 15:28:53 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//---------------------------------------------------------------
//
//  G4ParallelWorldScoringProcess.hh
//
//  Description:
//    This procss takes a parallel world and limits a step
//   on the boundaries of volumes in the parallel world.
//    It invokes sensitive detectors assigned in the parallel
//   world.
//
//---------------------------------------------------------------


#ifndef G4ParallelWorldScoringProcess_h
#define G4ParallelWorldScoringProcess_h 1

#include "globals.hh"
class G4Step;
class G4Navigator;
class G4TransportationManager;
class G4PathFinder;
class G4VTouchable;
class G4VPhysicalVolume;
class G4ParticleChange;
#include "G4VProcess.hh"
#include "G4FieldTrack.hh"
#include "G4TouchableHandle.hh"

//------------------------------------------
//
//        G4ParallelWorldScoringProcess class
//
//------------------------------------------


// Class Description:

class G4ParallelWorldScoringProcess : public G4VProcess
{
public: // with description

  //------------------------
  // Constructor/Destructor
  //------------------------
  
  G4ParallelWorldScoringProcess(const G4String& processName = "ParaWorldScore",
				 G4ProcessType theType = fParameterisation);
  virtual ~G4ParallelWorldScoringProcess();
  
  //--------------------------------------------------------------
  // Set Paralle World
  //--------------------------------------------------------------

  void SetParallelWorld(G4String parallelWorldName);
  void SetParallelWorld(G4VPhysicalVolume* parallelWorld);

  //--------------------------------------------------------------
  //     Process interface
  //--------------------------------------------------------------

  void StartTracking(G4Track*);
  
  //------------------------------------------------------------------------
  // GetPhysicalInteractionLength() and DoIt() methods for AtRest 
  //------------------------------------------------------------------------
  
  G4double AtRestGetPhysicalInteractionLength(
					      const G4Track& ,
					      G4ForceCondition* 
					      );

  G4VParticleChange* AtRestDoIt(
			       const G4Track& ,
			       const G4Step&
			       );

  //------------------------------------------------------------------------
  // GetPhysicalInteractionLength() and DoIt() methods for AlongStep 
  //------------------------------------------------------------------------
  
  G4double AlongStepGetPhysicalInteractionLength(
						 const G4Track&,
						 G4double  ,
						 G4double  ,
						 G4double&,
						 G4GPILSelection*
						 );

  G4VParticleChange* AlongStepDoIt(
				  const G4Track& ,
				  const G4Step& 
				  );

  //-----------------------------------------------------------------------
  // GetPhysicalInteractionLength() and DoIt() methods for PostStep
  //-----------------------------------------------------------------------
  
  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
						G4double   previousStepSize,
						G4ForceCondition* condition);
  
  G4VParticleChange* PostStepDoIt(const G4Track& ,const G4Step& );

private:
  void CopyStep(const G4Step & step);

  G4Step * fGhostStep;
  G4StepPoint * fGhostPreStepPoint;
  G4StepPoint * fGhostPostStepPoint;

  G4VParticleChange aDummyParticleChange;
  G4ParticleChange xParticleChange;

  G4TransportationManager* fTransportationManager;
  G4PathFinder*        fPathFinder;

  // -------------------------------
  // Navigation in the Ghost World:
  // -------------------------------
  G4String             fGhostWorldName;
  G4VPhysicalVolume*   fGhostWorld;
  G4Navigator*         fGhostNavigator;
  G4int                fNavigatorID;
  G4TouchableHandle    fOldGhostTouchable;
  G4TouchableHandle    fNewGhostTouchable;
  G4FieldTrack         fFieldTrack;
  G4double             fGhostSafety;
  G4bool               fOnBoundary;

  // ******************************************************
  // ******************************************************
  //
  //  For TESTS:
  //
  // ******************************************************
  // ******************************************************
public:
  void Verbose(const G4Step&) const;
};

#endif
